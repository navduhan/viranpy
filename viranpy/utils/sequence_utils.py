#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Sequence utility functions for ViRAnPy pipeline.
"""

import os
import re
from typing import Dict, List, Tuple, Optional
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Try-Except Block for Compatibility
try:
    from Bio.Alphabet import IUPAC
except ImportError:
    IUPAC = None


def get_sequence_length_from_fasta(contig_id: str, fasta_file: str) -> int:
    """
    Get the sequence length for a specific contig from a FASTA file.
    
    Args:
        contig_id: Contig identifier
        fasta_file: Path to FASTA file
        
    Returns:
        Length of the sequence
        
    Raises:
        ValueError: If contig not found
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == contig_id:
            return len(record.seq)
    raise ValueError(f"Contig {contig_id} not found in FASTA file.")


def rename_sequences(args, record_iter) -> Dict[str, str]:
    """
    Rename sequences and create temporary files.
    
    Args:
        args: Pipeline arguments containing locus and other settings
        record_iter: Iterator of sequence records
        
    Returns:
        Dictionary mapping new names to original names
    """
    from .file_utils import batch_iterator, string_split_by_numbers
    
    new_names_sequences = {}
    counter = 1
    
    for i, batch in enumerate(batch_iterator(record_iter, 1)):
        seq_index = f"LOC_{i + 1}"
        
        # Write temporary file
        with open(f"{seq_index}.temp.fna", "w") as handle:
            SeqIO.write(batch, handle, "fasta")
        
        # Process and rename sequences
        with open(f"{seq_index}.temp.fna", "r") as original, \
             open(f"{seq_index}.fna", "w") as corrected:
            
            if IUPAC:
                sequences = SeqIO.parse(original, "fasta", IUPAC.ambiguous_dna)
            else:
                sequences = SeqIO.parse(original, "fasta")
                
            for record in sequences:
                original_name = record.id
                record.id = f"{args.locus}_{counter}"
                record.description = record.description
                counter += 1
                new_names_sequences[record.id] = original_name
                print(f"WARNING: {original_name} was renamed as {record.id}", file=os.sys.stderr)
                SeqIO.write(record, corrected, "fasta")
        
        # Write log file
        with open("logfile.txt", "w") as logfile:
            logfile.write("#Original\tNew\n")
            for oldname in sorted(new_names_sequences, key=string_split_by_numbers):
                logfile.write(f"{new_names_sequences[oldname]}\t{oldname}\n")
        
        # Clean up temporary file
        os.remove(f"{seq_index}.temp.fna")
    
    return new_names_sequences


def get_combined_seqs(record, args, IUPAC=None) -> SeqRecord:
    """
    Get combined sequences for circularity detection.
    
    Args:
        record: Sequence record
        args: Pipeline arguments
        IUPAC: IUPAC alphabet (optional)
        
    Returns:
        Combined sequence record
    """
    seq_beginning = str(record.seq[0:args.gc_skew_read_length])
    seq_ending = str(record.seq[len(record.seq)-args.gc_skew_read_length:len(record.seq)])
    
    if IUPAC:
        return SeqRecord(
            Seq(seq_beginning + seq_ending, IUPAC.ambiguous_dna), 
            id=record.description
        )
    else:
        return SeqRecord(
            Seq(seq_beginning + seq_ending), 
            id=record.description
        )


def write_temp_file(combined_seqs: SeqRecord, filename: str = "temporal_circular.fasta") -> None:
    """
    Write combined sequences to a temporary file.
    
    Args:
        combined_seqs: Combined sequence record
        filename: Output filename
    """
    SeqIO.write(combined_seqs, filename, "fasta")


def validate_sequence(sequence: str) -> bool:
    """
    Validate if a sequence contains valid DNA/RNA characters.
    
    Args:
        sequence: Sequence string to validate
        
    Returns:
        True if valid, False otherwise
    """
    valid_chars = set('ATCGNMRWSYKVHDBXatcgnmrwsykvhdbx')
    return all(char in valid_chars for char in sequence)


def reverse_complement(sequence: str) -> str:
    """
    Get the reverse complement of a DNA sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Reverse complement sequence
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content of a sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        GC content as a percentage
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence)
    return (gc_count / total_count) * 100 if total_count > 0 else 0.0


def extract_sequence_region(sequence: str, start: int, end: int, strand: int = 1) -> str:
    """
    Extract a region from a sequence.
    
    Args:
        sequence: Input sequence
        start: Start position (0-based)
        end: End position (0-based)
        strand: Strand direction (1 for forward, -1 for reverse)
        
    Returns:
        Extracted sequence region
    """
    region = sequence[start:end]
    if strand == -1:
        region = reverse_complement(region)
    return region


def merge_sequences(sequences: List[str], gap: str = "N" * 100) -> str:
    """
    Merge multiple sequences with gaps.
    
    Args:
        sequences: List of sequences to merge
        gap: Gap sequence to insert between sequences
        
    Returns:
        Merged sequence
    """
    return gap.join(sequences)


def split_sequence(sequence: str, chunk_size: int, overlap: int = 0) -> List[str]:
    """
    Split a sequence into overlapping chunks.
    
    Args:
        sequence: Input sequence
        chunk_size: Size of each chunk
        overlap: Overlap between chunks
        
    Returns:
        List of sequence chunks
    """
    chunks = []
    for i in range(0, len(sequence), chunk_size - overlap):
        chunk = sequence[i:i + chunk_size]
        if len(chunk) >= chunk_size // 2:  # Only include chunks that are at least half the size
            chunks.append(chunk)
    return chunks 