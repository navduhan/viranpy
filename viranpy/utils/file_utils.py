"""
File utility functions for ViRAnPy pipeline.
"""

import os
import re
import shutil
import subprocess
import glob
from typing import List, Iterator, Any
from pathlib import Path


def batch_iterator(iterator: Iterator[Any], batch_size: int) -> Iterator[List[Any]]:
    """
    Create batches from an iterator.
    
    Args:
        iterator: Input iterator
        batch_size: Size of each batch
        
    Yields:
        List of items in each batch
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def cat_all(input_files: List[str], output_file: str) -> None:
    """
    Concatenate all input files into a single output file.
    
    Args:
        input_files: List of input file paths
        output_file: Output file path
    """
    with open(output_file, 'wb') as outfile:
        for fname in input_files:
            with open(fname, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)


def cmd_exists(cmd: str) -> bool:
    """
    Check if a command exists in the system PATH.
    
    Args:
        cmd: Command name to check
        
    Returns:
        True if command exists, False otherwise
    """
    return subprocess.call(
        ["which", cmd], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE
    ) == 0


def string_split_by_numbers(x: str) -> List[Any]:
    """
    Split a string numerically and not by alphabetic order.
    
    Args:
        x: String to split
        
    Returns:
        List of mixed strings and integers
    """
    r = re.compile('(\d+)')
    l = r.split(x)
    return [int(y) if y.isdigit() else y for y in l]


def remove_files(file_patterns: List[str], additional_patterns: List[str] = None) -> None:
    """
    Remove files matching specified patterns.
    
    Args:
        file_patterns: List of file patterns to remove
        additional_patterns: Additional glob patterns to remove
    """
    files_to_remove = file_patterns.copy()
    
    if additional_patterns:
        for pattern in additional_patterns:
            files_to_remove.extend(glob.glob(pattern))
    
    for file_path in files_to_remove:
        if os.path.isfile(file_path):
            os.remove(file_path)


def ensure_directory_exists(directory: str) -> None:
    """
    Ensure a directory exists, create if it doesn't.
    
    Args:
        directory: Directory path to ensure exists
    """
    Path(directory).mkdir(parents=True, exist_ok=True)


def get_file_extension(file_path: str) -> str:
    """
    Get the file extension from a file path.
    
    Args:
        file_path: Path to the file
        
    Returns:
        File extension (including the dot)
    """
    return Path(file_path).suffix


def get_file_basename(file_path: str) -> str:
    """
    Get the basename (filename without extension) from a file path.
    
    Args:
        file_path: Path to the file
        
    Returns:
        Basename without extension
    """
    return Path(file_path).stem


def is_file_empty(file_path: str) -> bool:
    """
    Check if a file is empty.
    
    Args:
        file_path: Path to the file
        
    Returns:
        True if file is empty, False otherwise
    """
    return os.path.getsize(file_path) == 0


def count_lines_in_file(file_path: str) -> int:
    """
    Count the number of lines in a file.
    
    Args:
        file_path: Path to the file
        
    Returns:
        Number of lines in the file
    """
    with open(file_path, 'r') as f:
        return sum(1 for _ in f) 