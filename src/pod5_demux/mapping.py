# src/pod5_demux/mapping.py
import os
import re
import gzip
import io
import glob
from pathlib import Path
from functools import partial
from multiprocessing import Pool
from typing import Tuple, Dict, Optional
import struct
import zlib

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

from Bio import SeqIO
from pod5_demux.utils import detect_format

BARCODE_RE = re.compile(r"(barcode\d+)")

# BGZF meta-header according to SAMv1 Section 4.1
# ID1 ID2 CM FLG(4B) + MTIME(4B) + XFL OS(2B) + XLEN(2B) + SI1 SI2 SLEN(2B)
_BGZF_HEADER = struct.Struct('<4BI2BH2BH')
 
# CRC32 + ISIZE at the end of each BGZF block (SAMv1 Section 4.1)
_BGZF_TAIL = struct.Struct('<2I')
 
# Fixed-size portion of a BAM alignment record (SAMv1 Section 4.2)
# refID(4) pos(4) l_read_name(1) mapq(1) bin(2) n_cigar_op(2) flag(2) l_seq(4)
# next_refID(4) next_pos(4) tlen(4) = Total 32 bytes
# We only use: l_read_name (byte 8), n_cigar_op (bytes 12-13), l_seq (bytes 16-19)
_BAM_RECORD_HEADER = struct.Struct('<4s4sBBHHHI4s4s4s')  # 32 bytes
 
 
# ---------------------------------------------------------------------------
# BGZF Decompression — SAMv1 Section 4.1
# ---------------------------------------------------------------------------
 
def _read_bgzf_blocks(file_path: str) -> bytes:
    """
    Decompresses a BGZF file into a single bytearray.
 
    According to SAMv1 Section 4.1: A BGZF file is a series of gzip blocks, 
    each containing an extra field with BSIZE (total block size - 1).
    Data is compressed using raw deflate (zlib without gzip header).
    EOF is signaled by an empty block with ISIZE == 0.
    """
    result = bytearray()
    hsize = _BGZF_HEADER.size  # 18 bytes
 
    with open(file_path, 'rb') as f:
        while True:
            raw = f.read(hsize)
            if len(raw) < hsize:
                break
 
            ID1, ID2, CM, FLG, MTIME, XFL, OS, XLEN, SI1, SI2, SLEN = _BGZF_HEADER.unpack(raw)
 
            # BGZF signature verification (SAMv1 Section 4.1)
            if not (ID1 == 31 and ID2 == 139 and CM == 8 and FLG == 4
                    and SI1 == 66 and SI2 == 67 and SLEN == 2):
                raise ValueError('Invalid BGZF header.')
 
            # BSIZE: total block size - 1
            bsize = struct.unpack('<H', f.read(2))[0]
 
            # Compressed data size: BSIZE+1 - 18(header) - 2(BSIZE field) - 8(CRC32+ISIZE)
            # = BSIZE - XLEN - 19 (SAMv1 Section 4.1: CDATA = BSIZE-XLEN-19)
            cdata = f.read(bsize - XLEN - 19)
            crc32, isize = _BGZF_TAIL.unpack(f.read(8))
 
            # EOF marker: empty block with ISIZE == 0 (SAMv1 Section 4.1.2)
            if isize == 0:
                break
 
            # Raw deflate decompression: zlib.decompressobj(-15) bypasses gzip header
            # This is the same approach used by htslib and bamnostic
            data = zlib.decompressobj(-15).decompress(cdata)
 
            # Integrity check
            if zlib.crc32(data) & 0xFFFFFFFF != crc32 or len(data) != isize:
                raise ValueError('BGZF: Integrity error (CRC32 or ISIZE mismatch).')
 
            result.extend(data)
 
    return bytes(result)
 
 
# ---------------------------------------------------------------------------
# Custom BAM Parser — reads only read_id and BC/RG tags
# ---------------------------------------------------------------------------
 
def _parse_bam_without_pysam(file_path: str, is_fmap: bool, clean_filename: str, known_bc: Optional[str]) -> Tuple[Dict, Dict]:
    """
    Reads a BAM file without the pysam dependency. 
    Extracts only the UUID (read_name) and the BC or RG optional tags.
 
    BAM record structure (SAMv1 Section 4.2):
      32B             Fixed header (refID, pos, l_read_name, mapq, bin,
                                    n_cigar_op, flag, l_seq, next_refID, next_pos, tlen)
      l_read_name B   read_name (NUL-terminated)
      n_cigar_op×4 B  cigar     ← skipped
      (l_seq+1)//2 B  seq       ← skipped
      l_seq B         qual      ← skipped
      optional tags   BC, RG    ← extracted
    """
    mapping = {}
    filename_map = {}
 
    try:
        buf = io.BytesIO(_read_bgzf_blocks(file_path))
 
        # --- BAM Header ---
        if buf.read(4) != b'BAM\x01':
            raise ValueError('Invalid BAM magic bytes.')
 
        l_text = struct.unpack('<I', buf.read(4))[0]
        buf.read(l_text)  # Plain-text SAM header
 
        n_ref = struct.unpack('<I', buf.read(4))[0]
        for _ in range(n_ref):
            l_name = struct.unpack('<I', buf.read(4))[0]
            buf.read(l_name + 4)  # Reference name + length
 
        # --- Alignment Records ---
        while True:
            bs = buf.read(4)
            if len(bs) < 4:
                break
            block_size = struct.unpack('<I', bs)[0]
            block = buf.read(block_size)
            if len(block) < block_size:
                break
 
            # Unpack the 32-byte fixed portion
            _, _, l_read_name, _, _, n_cigar_op, _, l_seq, _, _, _ = _BAM_RECORD_HEADER.unpack(block[:32])
 
            # read_name (UUID) — without the trailing NUL byte
            read_id = block[32: 32 + l_read_name - 1].decode('ascii', errors='replace')

            if known_bc:
                found_bc = known_bc
            else:
                # Offset where optional tags begin:
                # 32 (fixed) + l_read_name + n_cigar_op×4 + (l_seq+1)//2 + l_seq
                tag_offset = 32 + l_read_name + n_cigar_op * 4 + (l_seq + 1) // 2 + l_seq
    
                # Parse optional tags (SAMv1 Section 4.2.4)
                found_bc = _extract_bc_tag(block, tag_offset)
 
            if found_bc:
                mapping[read_id] = found_bc
                if is_fmap:
                    filename_map[read_id] = clean_filename
 
    except Exception as e:
        print(f'[Error] Failed to read BAM (custom parser) {file_path}: {e}')
 
    return mapping, filename_map
 
 
def _extract_bc_tag(block: bytes, offset: int) -> str|None:
    """
    Iterates through the optional tags of a BAM record and returns the 
    value of the BC (Barcode) or RG (Read Group) tag.
 
    Optional tag structure (SAMv1 Section 4.2.4):
      2B  tag name (e.g., b'BC')
      1B  val_type ('A','c','C','s','S','i','I','f','Z','H','B')
      xB  value (length depends on val_type)
    """
    # Byte sizes for fixed-length numeric tag types 
    TAG_VALUE_SIZES = {b'A': 1, b'c': 1, b'C': 1, b's': 2, b'S': 2,
                       b'i': 4, b'I': 4, b'f': 4}
 
    rg_val = None
    i = offset
 
    while i < len(block) - 2:
        tag      = block[i:i+2]
        val_type = block[i+2:i+3]
        i += 3
 
        if val_type in TAG_VALUE_SIZES:
            size = TAG_VALUE_SIZES[val_type]
            val = block[i:i+size]
            i += size
 
            if tag in (b'BC', b'RG'):
                # Numeric types (unlikely for BC/RG, but handled safely)
                decoded = val.decode('ascii', errors='replace').rstrip('\x00')
                if tag == b'BC':
                    return decoded
                elif tag == b'RG':
                    rg_val = decoded
 
        elif val_type in (b'Z', b'H'):
            # NUL-terminated string
            end = block.index(b'\x00', i)
            val = block[i:end].decode('ascii', errors='replace')
            i = end + 1
 
            if tag == b'BC':
                return val
            elif tag == b'RG' and rg_val is None:
                match = BARCODE_RE.search(val)
                rg_val = match.group() if match else val
 
        elif val_type == b'B':
            # Array of values: sub-type(1B) + count(4B) + count×size elements
            sub_type = block[i:i+1]
            count    = struct.unpack('<I', block[i+1:i+5])[0]
            size     = TAG_VALUE_SIZES.get(sub_type, 1)
            i += 5 + count * size
 
        else:
            break  # Unknown type encountered — stop parsing
 
    return rg_val  # BC not found, return RG (or None)


def parse_single_mapping_file(file_path: str, file_type: str, is_fmap: bool, known_bc: Optional[str]) -> Tuple[Dict, Dict, Dict]:
    """
    Parses a single sequence mapping file (BAM, SAM, or FASTQ) and extracts 
    read-to-barcode mappings. Automatically falls back to custom parsers 
    if pysam is unavailable.
    """
    mapping = {}
    filename_map = {}
    map_bc_count = {}
    
    clean_filename = Path(file_path).name
    for ext in [".fastq.gz", ".fastq", f".{file_type}"]:
        clean_filename = clean_filename.replace(ext, "")

    if file_type in ["sam", "bam"]:
        if not PYSAM_AVAILABLE and file_type == "bam":
            mapping, filename_map = _parse_bam_without_pysam(file_path, is_fmap, clean_filename, known_bc)
            
        elif not PYSAM_AVAILABLE and file_type == "sam":
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('@'): 
                        continue
                    parts = line.split('\t')
                    if len(parts) < 11: 
                        continue
                    
                    read_id = parts[0]
                    found_bc = None
                    
                    if known_bc:
                        found_bc = known_bc
                    else:
                        # Search optional fields for BC or RG tags
                        for field in parts[11:]:
                            if field.startswith("BC:Z:"):
                                found_bc = field.split(":")[2].strip()
                                break
                            elif field.startswith("RG:Z:"):
                                match = BARCODE_RE.search(field)
                                found_bc = match.group() if match else field
                    if found_bc:
                        mapping[read_id] = found_bc
                        if is_fmap: 
                            filename_map[read_id] = clean_filename
        else:
            # Parse using pysam if available
            mode = "rb" if file_type == "bam" else "r"
            try:
                with pysam.AlignmentFile(file_path, mode, check_sq=False) as samfile:
                    for read in samfile.fetch(until_eof=True):
                        read_id = read.query_name
                        if known_bc:
                            bc = known_bc
                        else:
                            try:
                                bc = read.get_tag("BC")
                            except KeyError:
                                try:
                                    rg = read.get_tag("RG")
                                    match = BARCODE_RE.search(str(rg))
                                    bc = match.group() if match else str(rg)
                                except KeyError:
                                    bc = None
                        if bc:
                            mapping[read_id] = bc
                            if is_fmap: 
                                filename_map[read_id] = clean_filename
            except Exception as e:
                print(f"[Error] Failed to read BAM/SAM with pysam {file_path}: {e}")

    elif file_type == "fastq":
        # Handle both gzipped and plain FASTQ files
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = "rt" if file_path.endswith('.gz') else "r"
        try:
            with open_func(file_path, mode) as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    if known_bc:
                        mapping[record.id] = known_bc
                    else:
                        match = BARCODE_RE.search(record.description)
                        if match:
                            mapping[record.id] = match.group()
                    if is_fmap: 
                        filename_map[record.id] = clean_filename
        except Exception as e:
            print(f"[Error] Failed to read FASTQ {file_path}: {e}")

    for bc in mapping.values():
        map_bc_count[bc] = map_bc_count.get(bc, 0) + 1
    return mapping, filename_map, map_bc_count

def load_barcode_map_parallel(input_path: str, output_format: str, known_bc: Optional[str], n_cores: int) -> Tuple[Dict, Dict, Dict]:
    """
    Locates mapping files in the input path and constructs a large 
    lookup dictionary in parallel using multiple CPU cores.
    """
    is_fmap = (output_format == "folder")
    input_format, path_type = detect_format(input_path)

    if not input_format:
        raise ValueError(f"Cannot detect input format for {input_path}")

    # Discover files based on detected format
    if input_format == "fastq":
        files = glob.glob(os.path.join(input_path, "**", f"*.{input_format}*"), recursive=True) if path_type == "dir" else [input_path]
    else:
        files = glob.glob(os.path.join(input_path, "**", f"*.{input_format}"), recursive=True) if path_type == "dir" else [input_path]

    print(f"-> Found {len(files)} mapping files ({input_format.upper()}). Loading in parallel...")

    barcode_map = {}
    filename_map = {}
    map_bc_count = {}
    
    worker = partial(parse_single_mapping_file, file_type=input_format, is_fmap=is_fmap, known_bc=known_bc)
    
    with Pool(processes=n_cores) as pool:
        results = pool.map(worker, files)
        for m, fm, bc_count in results:
            barcode_map.update(m)
            filename_map.update(fm)
            for bc, count in bc_count.items():
                map_bc_count[bc] = map_bc_count.get(bc, 0) + count
            
    return barcode_map, filename_map, map_bc_count