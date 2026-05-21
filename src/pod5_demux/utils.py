# src/pod5_demux/utils.py
import os
import shutil
from typing import Tuple

def ensure_unique_dir(path: str, overwrite: bool = False) -> str:
    """
    Checks whether the output directory is safe to use.
    - If it doesn't exist, it's safe.
    - If it exists but is empty, it's safe.
    - If it's not empty and overwrite=True, it wipes it clean.
    - If it's not empty and overwrite=False, raises FileExistsError.
    """
    base = os.path.abspath(path.rstrip("/\\"))
    
    if not os.path.exists(base):
        return base

    # Check if it's a directory and empty   
    if os.path.isdir(base) and not os.listdir(base):
        return base

    # Exists and is not empty
    # If overwrite is True, remove the existing directory and create a new one
    if overwrite:
        shutil.rmtree(base, ignore_errors=True)
        return base
    # If overwrite is False, raise an error
    raise FileExistsError(
        f"Output directory already exists and is NOT empty: '{base}'\n"
        f"  Use --overwrite to clear it automatically, or choose a different --output path."
    )

def detect_format(input_path: str) -> Tuple[str, str]:
    """
    Detects the format of the input data (BAM, SAM, FASTQ) 
    and whether the path is a directory or a single file.
    
    Args:
        input_path: Path to the input file or directory.
        
    Returns:
        A tuple (format_type, path_type). 
        Example: ('bam', 'dir') or ('fastq', 'file').
        Returns ('', '') if the path is invalid or format is unknown.
    """
    if not os.path.exists(input_path):
        return "", ""
        
    if os.path.isdir(input_path):
        for _, _, files in os.walk(input_path):
            for f_name in files:
                name_lower = f_name.lower()
                if name_lower.endswith(".bam"): 
                    return "bam", "dir"
                elif name_lower.endswith(".sam"): 
                    return "sam", "dir"
                elif name_lower.endswith((".fastq", ".fastq.gz")): 
                    return "fastq", "dir"
                    
    elif os.path.isfile(input_path):
        name_lower = input_path.lower()
        if name_lower.endswith(".bam"): 
            return "bam", "file"
        elif name_lower.endswith(".sam"): 
            return "sam", "file"
        elif name_lower.endswith((".fastq", ".fastq.gz")): 
            return "fastq", "file"
            
    return "", ""