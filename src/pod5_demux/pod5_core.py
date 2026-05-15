# src/pod5_demux/pod5_core.py
import os
import sys
import subprocess
from pathlib import Path
from typing import Optional, Tuple, Dict

import pod5

def split_one_pod5_by_barcode(
    pod5_file: str, 
    output_dir: str, 
    mapping: Dict[str, str], 
    filename_map: Dict[str, str], 
    output_mode: str, 
    known_bc: Optional[str]
) -> dict:
    """
    Reads a single POD5 file and splits its reads into temporary files 
    based on the provided barcode mapping.
    
    Args:
        pod5_file: Path to the input POD5 file.
        output_dir: Base directory for output.
        mapping: Dictionary mapping read UUIDs to barcode names.
        filename_map: Dictionary mapping read UUIDs to their original FASTQ/BAM filenames.
        output_mode: 'folder' or 'single_file'.
        known_bc: Optional specific barcode to filter by.
        
    Returns:
        A dictionary containing processing statistics (processed, sorted, unclassified).
    """
    source_name = Path(pod5_file).stem
    writers = {}
    stats = {"processed": 0, "sorted": 0, "unclassified": 0, "bc_counts": {}}
    temp_dir = os.path.join(output_dir, "temp_parts")

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

                # Determine the grouping key for the writer
                if output_mode == "folder" and filename_map:
                    origin_name = filename_map.get(read_id, "unknown")
                    writer_key = f"{bc_name}__{origin_name}"
                else:
                    writer_key = bc_name

                # Initialize a new writer if this is the first time seeing this key
                if writer_key not in writers:
                    out_name = f"{writer_key}__{source_name}.pod5"
                    out_path = os.path.join(temp_dir, out_name)
                    writers[writer_key] = pod5.Writer(out_path)
                    
                    # Note: We add the run_info for the first read. 
                    writers[writer_key]._add_run_info(read_record.run_info)

                try:
                    writers[writer_key].add_read(read_record.to_read())
                    if bc_name != "unclassified": 
                        stats["sorted"] += 1
                    else: 
                        stats["unclassified"] += 1
                    stats["bc_counts"][bc_name] = stats["bc_counts"].get(bc_name, 0) + 1
                except Exception as e:
                    print(f"[Error] Failed to write ID {read_id}: {e}")
    finally:
        # Ensure all writers are properly closed to avoid data corruption
        for w in writers.values():
            w.close()
            
    return stats

def merge_pod5_files(args: Tuple[list, str]):
    """
    Merges multiple temporary POD5 files into a single final file.
    Attempts a fast subprocess merge first. If that fails (e.g., missing CLI tool 
    or OS argument length limits), it falls back to a slower Python API merge.
    """
    input_files, output_file = args
    try:
        # Attempt 1: Call pod5 merge directly via CLI
        # Note: If `input_files` has thousands of items, this can hit the OS 
        # "Argument list too long" limit.
        subprocess.run(
            ["pod5", "merge"] + input_files + ["-o", output_file],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        try:
            # Attempt 2: Call via the current Python executable
            python_bin_dir = os.path.dirname(sys.executable)
            pod5_exe = os.path.join(python_bin_dir, "pod5")
            subprocess.run(
                [pod5_exe, "-m", "pod5", "merge"] + input_files + ["-o", output_file],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            # Fallback: Merge using the pod5 Python API (slowest method)
            print(f"  [Info] Subprocess merge failed, using Python API fallback for {Path(output_file).name}")
            
            with pod5.Writer(output_file) as writer:
                for file in input_files:
                    with pod5.Reader(file) as reader:
                        for i, read_rec in enumerate(reader.reads()):
                            if i == 0:
                                writer._add_run_info(read_rec.run_info)
                            writer.add_read(read_rec.to_read())
