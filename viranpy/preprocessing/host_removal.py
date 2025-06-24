#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Host removal using Bowtie2 for viral genome annotation.
"""

import os
import subprocess
import shutil
from typing import Dict, Any, List, Optional
from pathlib import Path

from ..utils.file_utils import cmd_exists


class HostRemover:
    """
    Host removal using Bowtie2 alignment.
    
    This module handles:
    - Bowtie2 index building for host genomes
    - Bowtie2 alignment for host read filtering
    - Fast host removal based on alignment
    """
    
    def __init__(self, config, logger=None):
        """Initialize the host remover."""
        self.config = config
        self.logger = logger
        self.bowtie2_results = {}
        self.filtered_files = {}
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        required_tools = ["bowtie2", "bowtie2-build"]
        missing_tools = [tool for tool in required_tools if not cmd_exists(tool)]
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        return True
    
    def build_host_index(self, host_genome: str, index_prefix: str = "host_index") -> Dict[str, Any]:
        """
        Build Bowtie2 index for host genome.
        
        Args:
            host_genome: Path to host genome FASTA file
            index_prefix: Prefix for Bowtie2 index files
            
        Returns:
            Dictionary containing index building results
        """
        if not os.path.exists(host_genome):
            return {"success": False, "error": f"Host genome file not found: {host_genome}"}
        
        try:
            cmd = [
                "bowtie2-build",
                host_genome,
                index_prefix,
                "--threads", str(self.config.ncpus or 1)
            ]
            
            subprocess.run(cmd, check=True)
            self.logger.info(f"Bowtie2 index built successfully: {index_prefix}")
            
            return {"success": True, "index_prefix": index_prefix}
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Bowtie2 index building failed: {e}")
            return {"success": False, "error": str(e)}
    
    def run_bowtie2_alignment(self, input_files: List[str], index_prefix: str, 
                             paired: bool = False, output_dir: str = "bowtie2_results") -> Dict[str, Any]:
        """
        Run Bowtie2 alignment against host genome.
        
        Args:
            input_files: List of input FASTQ files
            index_prefix: Bowtie2 index prefix
            paired: Whether the reads are paired-end
            output_dir: Output directory for Bowtie2 results
            
        Returns:
            Dictionary containing Bowtie2 alignment results
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        results = {}
        
        for i, input_file in enumerate(input_files):
            base_name = Path(input_file).stem
            sam_file = Path(output_dir) / f"{base_name}_host_alignment.sam"
            bam_file = Path(output_dir) / f"{base_name}_host_alignment.bam"
            
            cmd = [
                "bowtie2",
                "-x", index_prefix,
                "--threads", str(self.config.ncpus or 1),
                "--sensitive",  # Use sensitive mode for better sensitivity
                "-S", str(sam_file)
            ]
            
            if paired and i % 2 == 0:
                # For paired reads, process both files together
                if i + 1 < len(input_files):
                    cmd.extend(["-1", input_file, "-2", input_files[i + 1]])
                else:
                    cmd.extend(["-U", input_file])
            else:
                cmd.extend(["-U", input_file])
            
            try:
                subprocess.run(cmd, check=True)
                self.logger.info(f"Bowtie2 alignment completed for {input_file}")
                
                # Convert SAM to BAM and sort
                self._convert_sam_to_bam(sam_file, bam_file)
                
                # Parse alignment statistics
                stats = self._parse_bowtie2_stats(sam_file)
                results[input_file] = {
                    "sam_file": str(sam_file),
                    "bam_file": str(bam_file),
                    "statistics": stats
                }
                
            except subprocess.CalledProcessError as e:
                self.logger.error(f"Bowtie2 alignment failed for {input_file}: {e}")
                results[input_file] = {"success": False, "error": str(e)}
        
        return {"success": True, "results": results}
    
    def filter_host_reads(self, input_files: List[str], bowtie2_results: Dict[str, Any],
                         paired: bool = False, output_dir: str = "filtered_reads") -> Dict[str, Any]:
        """
        Filter out host reads based on Bowtie2 alignment.
        
        Args:
            input_files: List of input FASTQ files
            bowtie2_results: Bowtie2 alignment results
            paired: Whether the reads are paired-end
            output_dir: Output directory for filtered reads
            
        Returns:
            Dictionary containing filtering results
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        filtered_files = []
        
        for input_file in input_files:
            base_name = Path(input_file).stem
            filtered_file = Path(output_dir) / f"{base_name}_filtered.fastq"
            
            # Get alignment file from results
            alignment_info = bowtie2_results.get("results", {}).get(input_file, {})
            sam_file = alignment_info.get("sam_file")
            
            if sam_file and os.path.exists(sam_file):
                # Filter reads that don't align to host
                self._filter_reads_by_alignment(input_file, sam_file, filtered_file, paired)
                filtered_files.append(str(filtered_file))
            else:
                # If no alignment file, copy original file
                shutil.copy2(input_file, filtered_file)
                filtered_files.append(str(filtered_file))
        
        return {"success": True, "filtered_files": filtered_files}
    
    def _convert_sam_to_bam(self, sam_file: Path, bam_file: Path) -> None:
        """Convert SAM to BAM and sort."""
        try:
            # Convert SAM to BAM
            cmd = ["samtools", "view", "-bS", str(sam_file)]
            with open(bam_file, 'wb') as bam_out:
                subprocess.run(cmd, stdout=bam_out, check=True)
            
            # Sort BAM file
            sorted_bam = bam_file.with_suffix('.sorted.bam')
            cmd = ["samtools", "sort", str(bam_file), "-o", str(sorted_bam)]
            subprocess.run(cmd, check=True)
            
            # Replace original BAM with sorted BAM
            bam_file.unlink()
            sorted_bam.rename(bam_file)
            
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"SAM to BAM conversion failed: {e}")
    
    def _parse_bowtie2_stats(self, sam_file: Path) -> Dict[str, Any]:
        """Parse Bowtie2 alignment statistics from SAM file."""
        stats = {
            "total_reads": 0,
            "aligned_reads": 0,
            "unaligned_reads": 0,
            "alignment_rate": 0.0
        }
        
        try:
            with open(sam_file, 'r') as f:
                for line in f:
                    if not line.startswith('@'):  # Skip header lines
                        stats["total_reads"] += 1
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            flag = int(parts[1])
                            # Check if read is aligned (not unmapped)
                            if not (flag & 0x4):  # 0x4 = unmapped
                                stats["aligned_reads"] += 1
            
            stats["unaligned_reads"] = stats["total_reads"] - stats["aligned_reads"]
            if stats["total_reads"] > 0:
                stats["alignment_rate"] = (stats["aligned_reads"] / stats["total_reads"]) * 100
                
        except Exception as e:
            self.logger.warning(f"Could not parse Bowtie2 statistics: {e}")
        
        return stats
    
    def _filter_reads_by_alignment(self, input_file: str, sam_file: str, 
                                  filtered_file: Path, paired: bool) -> None:
        """Filter reads based on Bowtie2 alignment."""
        # Read SAM file to get aligned read IDs
        aligned_reads = set()
        
        with open(sam_file, 'r') as f:
            for line in f:
                if not line.startswith('@'):  # Skip header lines
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        read_id = parts[0]
                        flag = int(parts[1])
                        # Keep reads that are unmapped (don't align to host)
                        if flag & 0x4:  # 0x4 = unmapped
                            aligned_reads.add(read_id)
        
        # Filter FASTQ file
        with open(input_file, 'r') as infile, open(filtered_file, 'w') as outfile:
            keep_read = False
            for line in infile:
                if line.startswith('@'):
                    read_id = line.strip().split()[0][1:]  # Remove @ and get read ID
                    keep_read = read_id in aligned_reads
                
                if keep_read:
                    outfile.write(line) 