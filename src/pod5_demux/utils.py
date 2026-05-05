# src/pod5_demux/utils.py
import os
from typing import Tuple

def ensure_unique_dir(path: str) -> str:
    """
    Ensures a unique output directory name by appending a number 
    if the directory already exists.
    
    Args:
        path: The desired directory path.
        
    Returns:
        A unique directory path (e.g., 'dir', 'dir(1)', 'dir(2)').
    """
    base = path.rstrip("/\\")
    if not os.path.exists(base):
        return base
        
    i = 1
    while True:
        candidate = f"{base}({i})"
        if not os.path.exists(candidate):
            return candidate
        i += 1

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