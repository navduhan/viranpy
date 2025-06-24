#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Quality control and trimming for viral genome annotation.
"""

import os
import subprocess
import json
from typing import Dict, Any, List, Tuple, Optional
from pathlib import Path
import shutil

from ..core.base import BaseAnnotator
from ..core.results import AnnotationResult
from ..utils.file_utils import cmd_exists


class QualityController:
    """
    Quality control and trimming for viral genome sequencing data.
    
    This module handles:
    - FastQC quality assessment
    - Trim Galore quality trimming
    - Quality metrics reporting
    """
    
    def __init__(self, config, logger=None):
        """Initialize the quality controller."""
        self.config = config
        self.logger = logger
        self.quality_reports = {}
        self.trimmed_files = {}
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        required_tools = ["fastqc", "trim_galore"]
        missing_tools = [tool for tool in required_tools if not cmd_exists(tool)]
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        return True
    
    def run_fastqc(self, input_files: List[str], output_dir: str = "fastqc_reports") -> Dict[str, Any]:
        """
        Run FastQC quality assessment.
        
        Args:
            input_files: List of input FASTQ files
            output_dir: Output directory for FastQC reports
            
        Returns:
            Dictionary containing FastQC results
        """
        output_dir = os.path.join(self.config.root_output, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        cmd = ["fastqc", "-o", output_dir, "--noextract", "--threads", str(self.config.ncpus or 1)]
        cmd.extend(input_files)
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"FastQC completed for {len(input_files)} files")
            
            # Parse FastQC results
            results = self._parse_fastqc_results(output_dir, input_files)
            return results
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"FastQC failed: {e}")
            return {"success": False, "error": str(e)}
    
    def run_trim_galore(self, input_files: List[str], paired: bool = False, 
                       output_dir: str = "trimmed_reads") -> Dict[str, Any]:
        """
        Run Trim Galore quality trimming.
        
        Args:
            input_files: List of input FASTQ files
            paired: Whether the reads are paired-end
            output_dir: Output directory for trimmed reads
            
        Returns:
            Dictionary containing trimming results
        """
        output_dir = os.path.join(self.config.root_output, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        # Get trimming parameters from config
        quality = getattr(self.config, 'trim_quality', 20)
        length = getattr(self.config, 'min_length', 20)
        adapter = getattr(self.config, 'adapter', None)
        
        cmd = [
            "trim_galore",
            "--quality", str(quality),
            "--length", str(length),
            "--output_dir", output_dir,
            "--cores", str(self.config.ncpus or 1)
        ]
        
        if paired:
            cmd.append("--paired")
        
        if adapter:
            cmd.extend(["--adapter", adapter])
        
        cmd.extend(input_files)
        
        # Debug logging
        self.logger.info(f"Trim Galore input files: {input_files}")
        self.logger.info(f"Trim Galore paired mode: {paired}")
        self.logger.info(f"Trim Galore command: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"Trim Galore completed for {len(input_files)} files")
            
            # Get trimmed file paths
            trimmed_files = self._get_trimmed_files(output_dir, input_files, paired)
            self.logger.info(f"Trim Galore found trimmed files: {trimmed_files}")
            return {"success": True, "trimmed_files": trimmed_files}
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Trim Galore failed: {e}")
            return {"success": False, "error": str(e)}
    
    def _parse_fastqc_results(self, output_dir: str, input_files: List[str]) -> Dict[str, Any]:
        """Parse FastQC results and extract quality metrics."""
        results = {"success": True, "reports": {}}
        
        for input_file in input_files:
            base_name = Path(input_file).stem.replace('.fastq', '').replace('.fq', '')
            fastqc_dir = Path(output_dir) / f"{base_name}_fastqc"
            data_file = fastqc_dir / "fastqc_data.txt"
            
            if data_file.exists():
                metrics = self._extract_fastqc_metrics(data_file)
                results["reports"][input_file] = metrics
        
        return results
    
    def _extract_fastqc_metrics(self, data_file: Path) -> Dict[str, Any]:
        """Extract key metrics from FastQC data file."""
        metrics = {}
        
        with open(data_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('Total Sequences'):
                    metrics['total_sequences'] = int(line.split('\t')[1])
                elif line.startswith('Sequence length'):
                    metrics['sequence_length'] = line.split('\t')[1]
                elif line.startswith('%GC'):
                    metrics['gc_content'] = float(line.split('\t')[1])
                elif line.startswith('Total Deduplicated percentage'):
                    metrics['duplication_rate'] = 100 - float(line.split('\t')[1])
        
        return metrics
    
    def _get_trimmed_files(self, output_dir: str, input_files: List[str], paired: bool) -> List[str]:
        """Get paths to trimmed files."""
        trimmed_files = []
        
        if paired:
            # For paired-end, we need to find both R1 and R2 files
            # Trim Galore creates _val_1.fq.gz and _val_2.fq.gz for paired files
            # We'll look for all _val_*.fq.gz files and sort them
            output_path = Path(output_dir)
            val_files = list(output_path.glob("*_val_*.fq.gz"))
            val_files.sort()  # This will sort _val_1 before _val_2
            trimmed_files = [str(f) for f in val_files if f.exists()]
            self.logger.info(f"Found {len(trimmed_files)} paired trimmed files: {trimmed_files}")
        else:
            # For single-end, look for _trimmed.fq.gz files
            for input_file in input_files:
                base_name = Path(input_file).stem.replace('.fastq', '').replace('.fq', '')
                trimmed_file = Path(output_dir) / f"{base_name}_trimmed.fq.gz"
            
                if trimmed_file.exists():
                    trimmed_files.append(str(trimmed_file))
        
        return trimmed_files 