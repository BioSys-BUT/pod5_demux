# src/pod5_demux/pod5_core.py
import os
from pathlib import Path
from typing import Optional, Tuple

import pod5

# Barcode lookup table shared by workers (set once per worker process)
_BC_MAP = None


def _init_worker(bc_map):
    global _BC_MAP
    _BC_MAP = bc_map


def split_one_pod5_by_barcode(
    pod5_file: str,
    output_dir: str,
    output_mode: str,
    known_bc: Optional[str],
) -> dict:
    """
    Reads a single POD5 file and writes its reads, grouped by barcode.

    folder mode:      <output_dir>/<barcode>/<source>.pod5 (final files)
    single_file mode: <output_dir>/temp_parts/<barcode>__<source>.pod5 (merged later)

    The barcode lookup table is taken from the worker-global _BC_MAP.
    """
    mapping = _BC_MAP
    source_name = Path(pod5_file).stem
    temp_dir = os.path.join(output_dir, "temp_parts")

    writers = {}
    stats = {"processed": 0, "sorted": 0, "unclassified": 0, "bc_counts": {}, "n_files": 0}

    try:
        with pod5.Reader(pod5_file) as reader:
            for read_record in reader:
                stats["processed"] += 1
                read_id = str(read_record.read_id)
                bc_name = mapping.get(read_id, "unclassified")

                # If a specific barcode is enforced, skip all other reads
                if known_bc and bc_name != known_bc:
                    stats["unclassified"] += 1
                    continue

                # Open a writer the first time a barcode is seen in this file
                if bc_name not in writers:
                    if output_mode == "folder":
                        bc_dir = os.path.join(output_dir, bc_name)
                        os.makedirs(bc_dir, exist_ok=True)
                        out_path = os.path.join(bc_dir, f"{source_name}.pod5")
                    else:
                        out_path = os.path.join(temp_dir, f"{bc_name}__{source_name}.pod5")
                    writers[bc_name] = pod5.Writer(out_path)

                try:
                    writers[bc_name].add_read(read_record.to_read())
                    if bc_name != "unclassified":
                        stats["sorted"] += 1
                    else:
                        stats["unclassified"] += 1
                    stats["bc_counts"][bc_name] = stats["bc_counts"].get(bc_name, 0) + 1
                except Exception as e:
                    print(f"[Error] Failed to write ID {read_id}: {e}")
    finally:
        stats["n_files"] = len(writers)
        # Close all writers to avoid corrupted files
        for w in writers.values():
            w.close()

    return stats


def merge_pod5_files(args: Tuple[list, str]):
    """
    Merges multiple temporary POD5 files into one final file using pod5 Python API.
    """
    input_files, output_file = args
    if os.path.exists(output_file):
        os.remove(output_file)
    # Check that inputs exist
    for file in input_files:
        if not os.path.exists(file):
            raise FileNotFoundError(f"POD5 file not found: {file}")
    with pod5.Writer(output_file) as writer:
        for file in input_files:
            try:
                with pod5.Reader(file) as reader:
                    for read_rec in reader.reads():
                        writer.add_read(read_rec.to_read())
            except Exception as e:
                print(f"[Error] Failed to read POD5 file {file}: {e}")