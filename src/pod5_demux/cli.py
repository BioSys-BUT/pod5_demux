# src/pod5_demux/cli.py
import os
import glob
import time
import shutil
from functools import partial
from multiprocessing import Pool
from typing import Annotated, Optional
import re
from enum import Enum

import typer

# Import internal functions
from pod5_demux.utils import ensure_unique_dir
from pod5_demux.mapping import load_barcode_map_parallel
from pod5_demux.pod5_core import split_one_pod5_by_barcode, merge_pod5_files

try:
    # Change color for [required] parameters in Typer help
    from typer import rich_utils as _ru
    _ru.STYLE_REQUIRED_SHORT = "bright_red"
    _ru.STYLE_REQUIRED_LONG  = "bright_red"
except (ImportError, AttributeError):
    pass

app = typer.Typer(
    help="Optimized tool for demultiplexing POD5 files using a BAM/SAM/FASTQ map.",
    add_completion=False,
)
BARCODE_RE = re.compile(r"(barcode\d+)")

class OutputMode(str, Enum):
    """Allowed values for output_mode to ensure valid user input."""
    single_file = "single_file"
    folder = "folder"


def run_demultiplexing(
    input_seq: str, 
    input_pod5: str, 
    output_dir: str, 
    known_bc: Optional[str], 
    output_mode: OutputMode, 
    n_cores: int
):
    """
    Core function that orchestrates the entire demultiplexing workflow.
    It loads the mapping, splits the POD5 files by barcode, and then merges them.
    """
    start_time = time.perf_counter()

    output_dir = ensure_unique_dir(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print("*" * 50)
    print(" POD5 DEMULTIPLEXER")
    print("*" * 50)
    print(f"Input map:      {input_seq}")
    print(f"Input POD5:     {input_pod5}")
    print(f"Output folder:  {output_dir}")
    print(f"Output mode:    {output_mode}")
    print(f"Thread count:   {n_cores}")
    print("=" * 50)

    # LOAD MAP
    try:
        bc_map, f_map = load_barcode_map_parallel(input_seq, output_mode, known_bc, n_cores)
    except ValueError as e:
        print(f"Critical error: {e}")
        return

    if not bc_map:
        print("!!! Process terminated: No barcodes found in mapping files.")
        return

    map_time = time.perf_counter()
    print(f"-> Map loaded: {len(bc_map)} IDs (time: {map_time - start_time:.2f} s)")

    # FIND POD5 FILES
    pod5_files = glob.glob(os.path.join(input_pod5, "**", "*.pod5"), recursive=True) if os.path.isdir(input_pod5) else [input_pod5]
    if not pod5_files:
        print("!!! Process terminated: No source POD5 files found.")
        return

    # SPLIT PHASE
    # Create a temporary directory for split parts
    temp_dir = os.path.join(output_dir, "temp_parts")
    os.makedirs(temp_dir, exist_ok=True)

    print(f"-> Splitting {len(pod5_files)} POD5 files...")
    worker_pod5 = partial(
        split_one_pod5_by_barcode, 
        mapping=bc_map, 
        output_dir=output_dir, 
        output_mode=output_mode, 
        filename_map=f_map, 
        known_bc=known_bc
    )

    with Pool(n_cores) as p:
        results = p.map(worker_pod5, pod5_files)

    total_reads = sum(r["processed"] for r in results)
    process_time = time.perf_counter()
    print(f"-> Splitting complete: Processed {total_reads} reads (time: {process_time - map_time:.2f} s)")

    # MERGE PHASE
    print("-> Merging split files...")
    temp_files = glob.glob(os.path.join(temp_dir, "*.pod5"))

    groups = {}
    for f in temp_files:
        fname = os.path.basename(f)
        parts = fname.split("__")
        bc = parts[0]

        if output_mode == "folder":
            origin = parts[1] if len(parts) > 1 else "unknown"
            final_folder = os.path.join(output_dir, bc)
            final_name = f"{origin}.pod5"
        else:
            final_folder = output_dir
            final_name = f"{bc}.pod5"

        target_key = (final_folder, final_name)
        if target_key not in groups:
            groups[target_key] = []
        groups[target_key].append(f)

    # Create output folders and prepare arguments for parallel merge
    merge_args = []
    single_moves = []
    for (folder, name), parts in groups.items():
        os.makedirs(folder, exist_ok=True)
        final_path = os.path.join(folder, name)
        if len(parts) > 1:
            merge_args.append((parts, final_path))
        else:
            single_moves.append((parts[0], final_path))

    # Move single-part files without merging
    for src, dst in single_moves:
        shutil.move(src, dst)

    # Parallel merge of multi-part groups
    if merge_args:
        with Pool(n_cores) as p:
            p.map(merge_pod5_files, merge_args)

    shutil.rmtree(temp_dir, ignore_errors=True)

    end_time = time.perf_counter()

    # REPORT
    print("\n" + "="*50)
    print(" PROCESSING COMPLETE ")
    print("="*50)
    print(f"Total time:         {end_time - start_time:.2f} seconds")
    print(f"Split time:         {process_time - map_time:.2f} seconds")
    print(f"Merge time:         {end_time - process_time:.2f} seconds")

    with open(os.path.join(output_dir, "demux_stats.txt"), "w") as stats_file:
        stats_file.write(f"Processed reads: {total_reads}\n")
        stats_file.write(f"Total time: {end_time - start_time:.2f} s\n")
        stats_file.write(f"Mode: {output_mode.value}\n")
        stats_file.write(f"Unique barcodes: {len(set(bc_map.values()))}\n")
        stats_file.write(f"Created files: {len(groups)}\n")


@app.command()
def main(
    seq: Annotated[str, typer.Option("-s", "--seq",
                                     exists=True,
                                     resolve_path=True,
                                     readable=True,
                                     help="Path to the reference map (BAM/SAM/FASTQ).")],
    pod5_path: Annotated[str, typer.Option("-p", "--pod5",
                                     exists=True,
                                     resolve_path=True,
                                     readable=True,
                                     help="Path to the POD5 data.")],
    output: Annotated[str, typer.Option("-o", "--output", 
                                     help="Output directory. (Default: 'pod5_demux_output')")] = "pod5_demux_output",
    known_bc: Annotated[Optional[str], typer.Option("-b", "--bc", 
                                     help="Barcode name to assign to all mapped reads.")] = None,
    mode: Annotated[OutputMode, typer.Option("-m", "--mode", 
                                     help="Mode: 'single_file' or 'folder'.")] = OutputMode.folder,
    threads: Annotated[Optional[int], typer.Option("-t", "--threads", 
                                     help="Number of CPU threads used for processing.")] = 8
):
    """
    Starts the demultiplexing process.
    """
    run_demultiplexing(seq, pod5_path, output, known_bc, mode, threads)