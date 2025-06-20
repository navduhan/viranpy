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
        
        kraken_db = getattr(self.config, 'kraken2_db', None)
        if not kraken_db:
            self.logger.error("Kraken2 database not specified")
            return {"success": False, "error": "Kraken2 database not specified"}
        
        results = {}
        
        for i, input_file in enumerate(input_files):
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
            
            if paired and i % 2 == 0:
                # For paired reads, process both files together
                if i + 1 < len(input_files):
                    cmd.extend([input_file, input_files[i + 1]])
                else:
                    cmd.append(input_file)
            elif not paired:
                cmd.append(input_file)
            
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
        
        cmd = ["ktImportTaxonomy", "-o", str(krona_html)]
        cmd.extend(kraken_files)
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"Krona visualization created: {krona_html}")
            return str(krona_html)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Krona visualization failed: {e}")
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
    
    def generate_taxonomic_report(self, kraken_results: Dict[str, Any], 
                                krona_file: str = None) -> str:
        """Generate a comprehensive taxonomic classification report."""
        report_file = "taxonomic_classification_report.html"
        
        html_content = self._create_taxonomic_html_report(kraken_results, krona_file)
        
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        return report_file
    
    def _create_taxonomic_html_report(self, kraken_results: Dict[str, Any], 
                                    krona_file: str = None) -> str:
        """Create HTML taxonomic classification report."""
        html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Taxonomic Classification Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                .section { margin: 20px 0; padding: 10px; border: 1px solid #ddd; }
                .taxonomy-item { margin: 5px 0; padding: 5px; background: #f9f9f9; }
                .viral { background: #e8f5e8; border-left: 4px solid #28a745; }
                .bacterial { background: #fff3cd; border-left: 4px solid #ffc107; }
                .other { background: #f8d7da; border-left: 4px solid #dc3545; }
            </style>
        </head>
        <body>
            <h1>Taxonomic Classification Report</h1>
        """
        
        # Add Kraken2 results
        if kraken_results.get("success"):
            html += "<div class='section'><h2>Kraken2 Classification Results</h2>"
            for file, results in kraken_results.get("results", {}).items():
                html += f"<h3>{file}</h3>"
                for taxid, info in results.get("taxonomy", {}).items():
                    # Determine classification type for styling
                    rank = info['rank'].lower()
                    if 'virus' in info['name'].lower() or rank in ['viruses', 'viral']:
                        css_class = "viral"
                    elif rank in ['bacteria', 'bacterial', 'genus', 'species']:
                        css_class = "bacterial"
                    else:
                        css_class = "other"
                    
                    html += f"""
                    <div class='taxonomy-item {css_class}'>
                        <strong>{info['name']}</strong> ({info['rank']}) - {info['percentage']:.2f}%
                        <br><small>Reads: {info['reads']:,} | Tax reads: {info['tax_reads']:,}</small>
                    </div>
                    """
            html += "</div>"
        
        # Add Krona visualization link
        if krona_file:
            html += f"""
            <div class='section'>
                <h2>Interactive Taxonomic Visualization</h2>
                <p><a href="{krona_file}" target="_blank">Open Krona Interactive Chart</a></p>
            </div>
            """
        
        html += "</body></html>"
        return html 