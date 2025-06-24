#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Taxonomic classification using Kraken2 and Krona for viral genome annotation.
"""

import os
import subprocess
from typing import Dict, Any, List, Optional
from pathlib import Path

from ..utils.file_utils import cmd_exists


class TaxonomicClassifier:
    """
    Taxonomic classification using Kraken2 and Krona visualization.
    
    This module handles:
    - Kraken2 taxonomic classification
    - Krona visualization
    - Taxonomic analysis reporting
    """
    
    def __init__(self, config, logger=None):
        """Initialize the taxonomic classifier."""
        self.config = config
        self.logger = logger
        self.kraken_results = {}
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        required_tools = ["kraken2", "ktImportTaxonomy"]
        missing_tools = [tool for tool in required_tools if not cmd_exists(tool)]
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        return True
    
    def run_kraken2(self, input_files: List[str], paired: bool = False, 
                   output_dir: str = "kraken2_results") -> Dict[str, Any]:
        """
        Run Kraken2 taxonomic classification.
        
        Args:
            input_files: List of input FASTQ/FASTA files
            paired: Whether the reads are paired-end
            output_dir: Output directory for Kraken2 results
            
        Returns:
            Dictionary containing Kraken2 results
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        # Debug logging
        self.logger.info(f"Kraken2 input files: {input_files}")
        self.logger.info(f"Kraken2 paired mode: {paired}")
        self.logger.info(f"Kraken2 number of input files: {len(input_files)}")
        
        kraken_db = getattr(self.config, 'kraken2_db', None)
        if not kraken_db:
            self.logger.error("Kraken2 database not specified")
            return {"success": False, "error": "Kraken2 database not specified"}
        
        results = {}
        
        if paired:
            # For paired-end, process files in pairs (R1, R2)
            for i in range(0, len(input_files), 2):
                if i + 1 < len(input_files):
                    # Process R1 and R2 together
                    r1_file = input_files[i]
                    r2_file = input_files[i + 1]
                    base_name = Path(r1_file).stem
                    report_file = Path(output_dir) / f"{base_name}_kraken2.report"
                    output_file = Path(output_dir) / f"{base_name}_kraken2.out"
                    
                    cmd = [
                        "kraken2",
                        "--db", kraken_db,
                        "--output", str(output_file),
                        "--report", str(report_file),
                        "--threads", str(self.config.ncpus or 1)
                    ]
                    
                    cmd.extend([r1_file, r2_file])
                    self.logger.info(f"Kraken2 paired-end command: R1={r1_file}, R2={r2_file}")
                    self.logger.info(f"Kraken2 full command: {' '.join(cmd)}")
                    
                    try:
                        subprocess.run(cmd, check=True)
                        self.logger.info(f"Kraken2 completed for paired files: {r1_file} and {r2_file}")
                        
                        # Parse Kraken2 report
                        report_data = self._parse_kraken2_report(report_file)
                        results[r1_file] = report_data
                        
                    except subprocess.CalledProcessError as e:
                        self.logger.error(f"Kraken2 failed for paired files {r1_file} and {r2_file}: {e}")
                        results[r1_file] = {"success": False, "error": str(e)}
                else:
                    # Handle odd number of files (shouldn't happen with proper paired input)
                    self.logger.warning(f"Kraken2 paired mode but odd number of files. Skipping: {input_files[i]}")
        else:
            # For single-end, process each file individually
            for input_file in input_files:
                base_name = Path(input_file).stem
                report_file = Path(output_dir) / f"{base_name}_kraken2.report"
                output_file = Path(output_dir) / f"{base_name}_kraken2.out"
                
                cmd = [
                    "kraken2",
                    "--db", kraken_db,
                    "--output", str(output_file),
                    "--report", str(report_file),
                    "--threads", str(self.config.ncpus or 1)
                ]
                
                cmd.append(input_file)
                self.logger.info(f"Kraken2 single-end command: {input_file}")
                self.logger.info(f"Kraken2 full command: {' '.join(cmd)}")
                
                try:
                    subprocess.run(cmd, check=True)
                    self.logger.info(f"Kraken2 completed for {input_file}")
                    
                    # Parse Kraken2 report
                    report_data = self._parse_kraken2_report(report_file)
                    results[input_file] = report_data
                    
                except subprocess.CalledProcessError as e:
                    self.logger.error(f"Kraken2 failed for {input_file}: {e}")
                    results[input_file] = {"success": False, "error": str(e)}
        
        return {"success": True, "results": results}
    
    def create_krona_visualization(self, kraken_results: Dict[str, Any], 
                                  output_dir: str = "kraken2_results") -> str:
        """
        Create Krona visualization of taxonomic classification.
        
        Args:
            kraken_results: Kraken2 classification results
            output_dir: Directory containing Kraken2 output files
            
        Returns:
            Path to Krona HTML file
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        kraken_files = []
        
        for file_path in Path(output_dir).glob("*_kraken2.out"):
            kraken_files.append(str(file_path))
        
        if not kraken_files:
            self.logger.warning("No Kraken2 output files found for Krona visualization")
            return None
        
        krona_html = Path(output_dir) / "taxonomy_krona.html"
        
        cmd = ["ktImportTaxonomy", "-t", "5", "-m", "3", "-o", str(krona_html), ]
        cmd.extend(kraken_files)
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"Krona visualization created: {krona_html}")
            return str(krona_html)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Krona visualization failed: {e}")
            # Check if it's a taxonomy database issue
            if "Taxonomy not found" in str(e):
                self.logger.warning("Krona taxonomy database not found. To fix this, run: updateTaxonomy.sh")
                self.logger.info("Krona visualization skipped due to missing taxonomy database")
            else:
                self.logger.error(f"Krona visualization failed with error: {e}")
            return None
    
    def _parse_kraken2_report(self, report_file: Path) -> Dict[str, Any]:
        """Parse Kraken2 report file and extract taxonomic information."""
        results = {"taxonomy": {}, "summary": {}}
        
        with open(report_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    percentage = float(parts[0])
                    reads = int(parts[1])
                    tax_reads = int(parts[2])
                    rank = parts[3]
                    taxid = parts[4]
                    name = parts[5]
                    
                    if percentage > 0.1:  # Only include taxa with >0.1% abundance
                        results["taxonomy"][taxid] = {
                            "name": name,
                            "rank": rank,
                            "percentage": percentage,
                            "reads": reads,
                            "tax_reads": tax_reads
                        }
        
        return results 