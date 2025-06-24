#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Assembly statistics and QUAST integration for viral metagenomics.
"""

import os
import subprocess
import json
from typing import Dict, Any, List, Optional
from pathlib import Path

from .file_utils import cmd_exists


class AssemblyStatsCalculator:
    """
    Calculate assembly statistics using QUAST and other tools.
    
    This module provides functionality to:
    - Run QUAST for assembly quality assessment
    - Generate comprehensive assembly reports
    - Calculate viral-specific statistics
    """
    
    def __init__(self, config, logger=None):
        """Initialize the assembly stats calculator."""
        self.config = config
        self.logger = logger
        self.quast_results = {}
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        required_tools = ["quast.py"]
        missing_tools = [tool for tool in required_tools if not cmd_exists(tool)]
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        return True
    
    def run_quast_analysis(self, 
                          contigs_file: str,
                          output_dir: str = "quast_results",
                          reference_genome: Optional[str] = None,
                          min_contig_length: int = 200) -> Dict[str, Any]:
        """
        Run QUAST analysis on assembled contigs.
        
        Args:
            contigs_file: Input FASTA file containing contigs
            output_dir: Output directory for QUAST results
            reference_genome: Optional reference genome for comparison
            min_contig_length: Minimum contig length to consider
            
        Returns:
            Dictionary containing QUAST results
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # Build QUAST command
            cmd = [
                "quast.py",
                "-o", output_dir,
                "--min-contig", str(min_contig_length),
                "--threads", str(self.config.ncpus or 1),
                "--no-icarus",  # Skip Icarus viewer for faster processing
                "--no-snps",    # Skip SNP detection for metagenomes
                contigs_file
            ]
            
            # Add reference genome if provided
            if reference_genome and os.path.exists(reference_genome):
                cmd.extend(["-r", reference_genome])
                self.logger.info(f"Using reference genome: {reference_genome}")
            
            # Run QUAST
            self.logger.info("Running QUAST analysis")
            subprocess.run(cmd, check=True)
            
            # Parse QUAST results
            quast_results = self._parse_quast_results(output_dir)
            
            # Generate viral-specific statistics
            viral_stats = self._calculate_viral_specific_stats(contigs_file)
            
            # Combine results
            results = {
                "success": True,
                "quast_dir": output_dir,
                "quast_results": quast_results,
                "viral_stats": viral_stats,
            }
            
            return results
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"QUAST analysis failed: {e}")
            return {"success": False, "error": str(e)}
        except Exception as e:
            self.logger.error(f"Assembly statistics calculation failed: {e}")
            return {"success": False, "error": str(e)}
    
    def _parse_quast_results(self, quast_dir: str) -> Dict[str, Any]:
        """Parse QUAST results from the output directory."""
        results = {}
        
        try:
            # Parse transposed_report.tsv
            report_file = Path(quast_dir) / "transposed_report.tsv"
            if report_file.exists():
                with open(report_file, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                metric = parts[0]
                                value = parts[1]
                                # Try to convert to appropriate type
                                try:
                                    if '.' in value:
                                        results[metric] = float(value)
                                    else:
                                        results[metric] = int(value)
                                except ValueError:
                                    results[metric] = value
            
            # Parse contigs_report.tsv for detailed contig information
            contigs_file = Path(quast_dir) / "contigs_report.tsv"
            if contigs_file.exists():
                contigs_data = []
                with open(contigs_file, 'r') as f:
                    header = f.readline().strip().split('\t')
                    for line in f:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= len(header):
                                contig_info = dict(zip(header, parts))
                                contigs_data.append(contig_info)
                results["contigs_detailed"] = contigs_data
            
            return results
            
        except Exception as e:
            self.logger.warning(f"Could not parse QUAST results: {e}")
            return {}
    
    def _calculate_viral_specific_stats(self, contigs_file: str) -> Dict[str, Any]:
        """Calculate viral-specific assembly statistics."""
        viral_stats = {
            "total_contigs": 0,
            "total_bases": 0,
            "contig_lengths": [],
            "gc_content": [],
            "n50": 0,
            "n90": 0,
            "largest_contig": 0,
            "smallest_contig": 0,
            "mean_contig_length": 0,
            "median_contig_length": 0
        }
        
        try:
            from Bio import SeqIO
            from Bio.SeqUtils import GC
            
            contig_lengths = []
            gc_contents = []
            
            for record in SeqIO.parse(contigs_file, "fasta"):
                length = len(record.seq)
                contig_lengths.append(length)
                
                # Safely calculate GC content with error handling
                try:
                    # Convert sequence to string and ensure it's valid DNA
                    seq_str = str(record.seq).upper()
                    # Remove any non-DNA characters
                    seq_str = ''.join(c for c in seq_str if c in 'ATCGN')
                    if seq_str:
                        gc_content = GC(seq_str)
                    else:
                        gc_content = 0.0
                except Exception as gc_error:
                    self.logger.warning(f"Could not calculate GC content for contig {record.id}: {gc_error}")
                    gc_content = 0.0
                
                gc_contents.append(gc_content)
                viral_stats["total_bases"] += length
            
            viral_stats["total_contigs"] = len(contig_lengths)
            viral_stats["contig_lengths"] = contig_lengths
            viral_stats["gc_content"] = gc_contents
            
            if contig_lengths:
                contig_lengths.sort(reverse=True)
                viral_stats["largest_contig"] = contig_lengths[0]
                viral_stats["smallest_contig"] = contig_lengths[-1]
                viral_stats["mean_contig_length"] = sum(contig_lengths) / len(contig_lengths)
                viral_stats["median_contig_length"] = contig_lengths[len(contig_lengths)//2]
                
                # Calculate N50 and N90
                total_bases = sum(contig_lengths)
                cumulative_bases = 0
                n50_found = False
                n90_found = False
                
                for length in contig_lengths:
                    cumulative_bases += length
                    if not n50_found and cumulative_bases >= total_bases * 0.5:
                        viral_stats["n50"] = length
                        n50_found = True
                    if not n90_found and cumulative_bases >= total_bases * 0.9:
                        viral_stats["n90"] = length
                        n90_found = True
                    if n50_found and n90_found:
                        break
            
            if gc_contents:
                viral_stats["mean_gc_content"] = sum(gc_contents) / len(gc_contents)
            
        except Exception as e:
            self.logger.warning(f"Could not calculate viral-specific stats: {e}")
        
        return viral_stats
    
    def get_assembly_summary(self, quast_results: Dict[str, Any], viral_stats: Dict[str, Any]) -> Dict[str, Any]:
        """Get summary statistics for assembly."""
        summary = {
            "total_contigs": viral_stats.get('total_contigs', 0),
            "total_bases": viral_stats.get('total_bases', 0),
            "n50": viral_stats.get('n50', 0),
            "largest_contig": viral_stats.get('largest_contig', 0),
            "mean_contig_length": viral_stats.get('mean_contig_length', 0),
            "gc_content": viral_stats.get('mean_gc_content', 0)
        }
        
        # Add QUAST metrics if available
        for key in ['GC (%)', 'N50', 'Largest contig', 'Total length']:
            if key in quast_results:
                summary[f"quast_{key.lower().replace(' ', '_').replace('(', '').replace(')', '')}"] = quast_results[key]
        
        return summary 