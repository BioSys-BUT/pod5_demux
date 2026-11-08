# src/pod5_demux/mapping.py
import os
import re
import gzip
import glob
from pathlib import Path
from functools import partial
from multiprocessing import Pool
from typing import Tuple, Dict, Optional
import pysam


from Bio import SeqIO
from pod5_demux.utils import detect_format

BARCODE_RE = re.compile(r"(barcode\d+)")

def parse_single_mapping_file(
    file_path: str,
    file_type: str,
    is_fmap: bool,
    path_bc: Optional[str],
    known_bc: Optional[str],
) -> Tuple[Dict, Dict, Dict]:
    """
    Parses a single sequence mapping file (BAM, SAM, or FASTQ) and extracts
    read-to-barcode mappings. 
    Barcode assignment priority (applied per read):
      1. Annotation inside the file (BC tag, or RG tag only if it matches
         barcodeXX; barcode in FASTQ header).
         If a BC/RG tag exists but does not match barcodeXX (e.g. "unclassified"),
         the read is treated as explicitly unclassified and known_bc is NOT applied.
      2. Barcode inferred from path (file or parent directory named barcodeXX)
      3. Explicit user fallback (--bc / known_bc) — only when no tag was present
    """
    mapping = {}
    filename_map = {}
    map_bc_count = {}
    
    clean_filename = Path(file_path).name
    for ext in [".fastq.gz", ".fastq", f".{file_type}"]:
        clean_filename = clean_filename.replace(ext, "")

    if file_type in ["sam", "bam"]:
        # Parse using pysam
        mode = "rb" if file_type == "bam" else "r"
        try:
            with pysam.AlignmentFile(file_path, mode, check_sq=False) as samfile:
                for read in samfile.fetch(until_eof=True):
                    read_id = read.query_name
                    found_bc = None
                    tag_present = False

                    try:
                        bc_val = read.get_tag("BC")
                        tag_present = True
                        match = BARCODE_RE.search(str(bc_val))
                        found_bc = match.group() if match else None
                    except KeyError:
                        try:
                            rg = read.get_tag("RG")
                            tag_present = True
                            match = BARCODE_RE.search(str(rg))
                            found_bc = match.group() if match else None
                        except KeyError:
                            pass

                    if not found_bc:
                        found_bc = path_bc
                    if not found_bc and not tag_present:
                        found_bc = known_bc

                    if found_bc:
                        mapping[read_id] = found_bc
                        if is_fmap:
                            filename_map[read_id] = clean_filename

        except Exception as e:
            print(f"[Error] Failed to read BAM/SAM with pysam {file_path}: {e}")
            raise

    elif file_type == "fastq":
        open_func = gzip.open if file_path.endswith('.gz') else open
        mode = "rt" if file_path.endswith('.gz') else "r"
        try:
            with open_func(file_path, mode) as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    found_bc = None
                    tag_present = False

                    # For FASTQ: check if header contains any barcode-like field.
                    # "barcode=unclassified" or similar counts as tag_present.
                    if "barcode" in record.description.lower():
                        tag_present = True
                        match = BARCODE_RE.search(record.description)
                        found_bc = match.group() if match else None

                    if not found_bc:
                        found_bc = path_bc
                    if not found_bc and not tag_present:
                        found_bc = known_bc

                    if found_bc:
                        mapping[record.id] = found_bc
                        if is_fmap:
                            filename_map[record.id] = clean_filename

        except Exception as e:
            print(f"[Error] Failed to read FASTQ {file_path}: {e}")
            raise

    for bc in mapping.values():
        map_bc_count[bc] = map_bc_count.get(bc, 0) + 1

    return mapping, filename_map, map_bc_count


def _detect_bc_from_path(input_path: str) -> Optional[str]:
    """
    Tries to detect a barcode name from the input path.
    
    Handles two cases:
      1. The path is a directory named barcodeXX   → e.g. /data/fastq/barcode05/
      2. The path is a file named barcodeXX.fastq  → e.g. /data/fastq/barcode05.fastq
    
    Returns the detected barcode string (e.g. 'barcode05') or None.
    """
    p = Path(input_path)
    name = p.name if p.is_dir() else p.stem
    name = name.replace(".fastq", "")
    match = BARCODE_RE.fullmatch(name)
    return match.group() if match else None


def load_barcode_map_parallel(
    input_path: str,
    output_format: str,
    known_bc: Optional[str],
    n_cores: int,
) -> Tuple[Dict, Dict, Dict]:
    """
    Locates mapping files in the input path and constructs a large
    lookup dictionary in parallel using multiple CPU cores.

    Barcode assignment priority (applied per read, inside each worker):
      1. Annotation inside the file (BC/RG tag in BAM/SAM; barcode in FASTQ header)
         Reads with a BC/RG tag that does not match barcodeXX (e.g. "unclassified")
         are treated as explicitly unclassified — known_bc is NOT applied to them.
      2. Barcode inferred from path (file or parent directory named barcodeXX)
      3. Explicit user fallback (--bc / known_bc) — only when no tag was present
    """
    is_fmap = (output_format == "folder")
    input_format, path_type = detect_format(input_path)

    if not input_format:
        raise ValueError(f"Cannot detect input format for {input_path}")

    if input_format == "fastq":
        files = glob.glob(os.path.join(input_path, "**", f"*.{input_format}*"), recursive=True) if path_type == "dir" else [input_path]
        files = [f for f in files if f.lower().endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))]
    else:
        files = glob.glob(os.path.join(input_path, "**", f"*.{input_format}"), recursive=True) if path_type == "dir" else [input_path]
        files = [f for f in files if f.lower().endswith(f".{input_format}")]

    print(f"-> Found {len(files)} mapping files ({input_format.upper()}). Loading in parallel...")

    def make_worker_args(file_path: str) -> tuple:
        path_bc = _detect_bc_from_path(file_path)
        if not path_bc:
            path_bc = _detect_bc_from_path(str(Path(file_path).parent))
        return (file_path, input_format, is_fmap, path_bc, known_bc)

    worker_args = [make_worker_args(f) for f in files]

    barcode_map = {}
    filename_map = {}
    map_bc_count = {}
    duplicates = 0
    conflicts = 0

    with Pool(processes=n_cores) as pool:
        results = pool.starmap(parse_single_mapping_file, worker_args)
        for m, fm, bc_count in results:
            for read_id, bc in m.items():
                if read_id in barcode_map:
                    duplicates += 1
                    if barcode_map[read_id] != bc:
                        conflicts += 1
                else:
                    barcode_map[read_id] = bc

            filename_map.update(fm)

    # Count barcodes in the final merged mapping
    for bc in barcode_map.values():
        map_bc_count[bc] = map_bc_count.get(bc, 0) + 1

    if duplicates > 0:
        if conflicts > 0:
            print(f"[Warning] Found {duplicates} duplicate reads in mapping files ({conflicts} with conflicting barcodes)!")
        else:
            print(f"[Warning] Found {duplicates} duplicate reads in mapping files!")

    return barcode_map, filename_map, map_bc_count