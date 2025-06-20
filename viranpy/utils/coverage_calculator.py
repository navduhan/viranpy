#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Coverage calculation utilities for viral metagenomic contigs.
"""

import os
import subprocess
import re
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
from collections import defaultdict

from .file_utils import cmd_exists


class CoverageCalculator:
    """
    Calculate coverage statistics for viral metagenomic contigs.
    
    This module provides functionality to:
    - Map reads to contigs using BWA
    - Calculate coverage statistics using samtools
    - Generate coverage reports
    """
    
    def __init__(self, config, logger=None):
        """Initialize the coverage calculator."""
        self.config = config
        self.logger = logger
        self.coverage_results = {}
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        required_tools = ["bwa", "samtools"]
        missing_tools = [tool for tool in required_tools if not cmd_exists(tool)]
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        return True
    
    def calculate_coverage(self, 
                          contigs_file: str,
                          read_files: List[str],
                          paired: bool = False,
                          output_dir: str = "coverage_results") -> Dict[str, Any]:
        """
        Calculate coverage for contigs using BWA and samtools.
        
        Args:
            contigs_file: Input FASTA file containing contigs
            read_files: List of input FASTQ files
            paired: Whether reads are paired-end
            output_dir: Output directory for coverage results
            
        Returns:
            Dictionary containing coverage results
        """
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # Step 1: Index contigs with BWA
            self.logger.info("Indexing contigs with BWA")
            index_result = self._index_contigs(contigs_file)
            if not index_result:
                return {"success": False, "error": "Failed to index contigs"}
            
            # Step 2: Map reads to contigs
            self.logger.info("Mapping reads to contigs")
            bam_file = self._map_reads_to_contigs(contigs_file, read_files, paired, output_dir)
            if not bam_file:
                return {"success": False, "error": "Failed to map reads"}
            
            # Step 3: Calculate coverage statistics
            self.logger.info("Calculating coverage statistics")
            coverage_stats = self._calculate_coverage_stats(bam_file, contigs_file)
            
            # Step 4: Generate coverage report
            coverage_report = self._generate_coverage_report(coverage_stats, output_dir)
            
            return {
                "success": True,
                "bam_file": bam_file,
                "coverage_stats": coverage_stats,
                "coverage_report": coverage_report
            }
            
        except Exception as e:
            self.logger.error(f"Coverage calculation failed: {e}")
            return {"success": False, "error": str(e)}
    
    def _index_contigs(self, contigs_file: str) -> bool:
        """Index contigs using BWA."""
        try:
            cmd = ["bwa", "index", contigs_file]
            subprocess.run(cmd, check=True)
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BWA indexing failed: {e}")
            return False
    
    def _map_reads_to_contigs(self, 
                             contigs_file: str,
                             read_files: List[str],
                             paired: bool,
                             output_dir: str) -> Optional[str]:
        """Map reads to contigs using BWA."""
        try:
            sam_file = Path(output_dir) / "mapped_reads.sam"
            bam_file = Path(output_dir) / "mapped_reads.bam"
            sorted_bam = Path(output_dir) / "mapped_reads.sorted.bam"
            
            # BWA mem command
            cmd = [
                "bwa", "mem",
                "-t", str(getattr(self.config, 'ncpus', 1)),
                contigs_file
            ]
            
            if paired and len(read_files) >= 2:
                cmd.extend(read_files)
            else:
                cmd.extend(read_files)
            
            # Run BWA mem
            with open(sam_file, 'w') as sam_out:
                subprocess.run(cmd, stdout=sam_out, check=True)
            
            # Convert SAM to BAM
            cmd = ["samtools", "view", "-bS", str(sam_file)]
            with open(bam_file, 'wb') as bam_out:
                subprocess.run(cmd, stdout=bam_out, check=True)
            
            # Sort BAM file
            cmd = ["samtools", "sort", str(bam_file), "-o", str(sorted_bam)]
            subprocess.run(cmd, check=True)
            
            # Index sorted BAM
            cmd = ["samtools", "index", str(sorted_bam)]
            subprocess.run(cmd, check=True)
            
            # Clean up intermediate files
            sam_file.unlink(missing_ok=True)
            bam_file.unlink(missing_ok=True)
            
            return str(sorted_bam)
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Read mapping failed: {e}")
            return None
    
    def _calculate_coverage_stats(self, bam_file: str, contigs_file: str) -> Dict[str, Any]:
        """Calculate coverage statistics using samtools."""
        coverage_stats = {}
        
        try:
            # Get coverage depth using samtools depth
            cmd = ["samtools", "depth", "-a", bam_file]
            depth_output = subprocess.check_output(cmd, text=True)
            
            # Parse depth output
            contig_coverage = defaultdict(list)
            for line in depth_output.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        contig = parts[0]
                        position = int(parts[1])
                        depth = int(parts[2])
                        contig_coverage[contig].append(depth)
            
            # Calculate statistics for each contig
            for contig, depths in contig_coverage.items():
                if depths:
                    coverage_stats[contig] = {
                        "mean_coverage": sum(depths) / len(depths),
                        "median_coverage": sorted(depths)[len(depths)//2],
                        "min_coverage": min(depths),
                        "max_coverage": max(depths),
                        "total_positions": len(depths),
                        "covered_positions": len([d for d in depths if d > 0]),
                        "coverage_breadth": len([d for d in depths if d > 0]) / len(depths) * 100
                    }
            
            # Get contig lengths
            contig_lengths = self._get_contig_lengths(contigs_file)
            for contig in coverage_stats:
                if contig in contig_lengths:
                    coverage_stats[contig]["contig_length"] = contig_lengths[contig]
            
            return coverage_stats
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Coverage calculation failed: {e}")
            return {}
    
    def _get_contig_lengths(self, contigs_file: str) -> Dict[str, int]:
        """Get contig lengths from FASTA file."""
        contig_lengths = {}
        
        try:
            from Bio import SeqIO
            for record in SeqIO.parse(contigs_file, "fasta"):
                contig_lengths[record.id] = len(record.seq)
        except Exception as e:
            self.logger.warning(f"Could not parse contig lengths: {e}")
        
        return contig_lengths
    
    def _generate_coverage_report(self, coverage_stats: Dict[str, Any], output_dir: str) -> str:
        """Generate coverage report."""
        report_file = Path(output_dir) / "coverage_report.tsv"
        
        with open(report_file, 'w') as f:
            # Write header
            f.write("Contig_ID\tContig_Length\tMean_Coverage\tMedian_Coverage\t"
                   "Min_Coverage\tMax_Coverage\tCovered_Positions\tTotal_Positions\t"
                   "Coverage_Breadth(%)\n")
            
            # Write data
            for contig, stats in coverage_stats.items():
                f.write(f"{contig}\t{stats.get('contig_length', 'N/A')}\t"
                       f"{stats.get('mean_coverage', 0):.2f}\t"
                       f"{stats.get('median_coverage', 0)}\t"
                       f"{stats.get('min_coverage', 0)}\t"
                       f"{stats.get('max_coverage', 0)}\t"
                       f"{stats.get('covered_positions', 0)}\t"
                       f"{stats.get('total_positions', 0)}\t"
                       f"{stats.get('coverage_breadth', 0):.2f}\n")
        
        return str(report_file)
    
    def get_coverage_summary(self, coverage_stats: Dict[str, Any]) -> Dict[str, Any]:
        """Get summary statistics for coverage."""
        if not coverage_stats:
            return {}
        
        all_mean_coverages = [stats.get('mean_coverage', 0) for stats in coverage_stats.values()]
        all_breadths = [stats.get('coverage_breadth', 0) for stats in coverage_stats.values()]
        all_lengths = [stats.get('contig_length', 0) for stats in coverage_stats.values()]
        
        summary = {
            "total_contigs": len(coverage_stats),
            "total_bases": sum(all_lengths),
            "mean_coverage_across_contigs": sum(all_mean_coverages) / len(all_mean_coverages) if all_mean_coverages else 0,
            "median_coverage_across_contigs": sorted(all_mean_coverages)[len(all_mean_coverages)//2] if all_mean_coverages else 0,
            "mean_breadth_across_contigs": sum(all_breadths) / len(all_breadths) if all_breadths else 0,
            "high_coverage_contigs": len([c for c in all_mean_coverages if c >= 10]),
            "low_coverage_contigs": len([c for c in all_mean_coverages if c < 5]),
            "well_covered_contigs": len([b for b in all_breadths if b >= 90])
        }
        
        return summary 