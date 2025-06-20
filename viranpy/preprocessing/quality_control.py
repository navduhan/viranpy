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
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"Trim Galore completed for {len(input_files)} files")
            
            # Get trimmed file paths
            trimmed_files = self._get_trimmed_files(output_dir, input_files, paired)
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
        
        for input_file in input_files:
            base_name = Path(input_file).stem.replace('.fastq', '').replace('.fq', '')
            if paired:
                trimmed_file = Path(output_dir) / f"{base_name}_val_1.fq.gz"
            else:
                trimmed_file = Path(output_dir) / f"{base_name}_trimmed.fq.gz"
            
            if trimmed_file.exists():
                trimmed_files.append(str(trimmed_file))
        
        return trimmed_files
    
    def generate_quality_report(self, fastqc_results: Dict[str, Any], 
                              trim_results: Dict[str, Any]) -> str:
        """Generate a comprehensive quality control report."""
        # Write report inside output directory
        report_file = os.path.join(self.config.root_output, "quality_control_report.html")
        html_content = self._create_comprehensive_quality_report(fastqc_results, trim_results)
        with open(report_file, 'w') as f:
            f.write(html_content)
        return report_file

    def _create_comprehensive_quality_report(self, fastqc_results: Dict[str, Any], trim_results: Dict[str, Any]) -> str:
        """Create a modern, comprehensive HTML quality control report."""
        import base64
        import io
        html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Quality Control Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }
                .container { max-width: 900px; margin: 0 auto; background: white; padding: 24px; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
                h1, h2 { color: #007acc; }
                .section { margin: 30px 0; padding: 18px; border: 1px solid #ddd; border-radius: 6px; background: #fafafa; }
                .metric { margin: 8px 0; padding: 6px; background: #fff; border-left: 4px solid #007acc; }
                .success { color: #28a745; }
                .warning { color: #ffc107; }
                .error { color: #dc3545; }
                .file-list { background: #f8f9fa; padding: 10px; border-radius: 4px; margin: 10px 0; }
                .summary { background: #e7f3ff; padding: 15px; border-radius: 5px; margin: 20px 0; }
                .plot-img { display: block; margin: 0 auto 20px auto; max-width: 400px; }
            </style>
        </head>
        <body>
            <div class='container'>
                <h1>ViRAnPy Quality Control Report</h1>
                <div class='summary'>
                    <h2>Summary</h2>
                    <p><strong>Status:</strong> <span class='success'>QC Completed</span></p>
                    <p><strong>Input files:</strong> """
        if hasattr(self.config, 'input_file') and self.config.input_file:
            html += f"{self.config.input_file}"
        elif hasattr(self.config, 'read_files') and self.config.read_files:
            html += ', '.join(self.config.read_files)
        else:
            html += "N/A"
        html += """</p>
                    <p><strong>Parameters:</strong> Quality threshold = {}</p>
                </div>
        """.format(getattr(self.config, 'trim_quality', 20))

        # --- Collect QC metrics ---
        qc_stats = []
        fastqc_reports = fastqc_results.get("reports", {}) if fastqc_results.get("success") else {}
        trimmed_files = trim_results.get('trimmed_files', []) if trim_results.get('success') else []
        trimming_stats = []
        paired = len(fastqc_reports) == 2
        # Try to match input/trimmed files by order
        for i, (infile, metrics) in enumerate(fastqc_reports.items()):
            stat = {
                'file': infile,
                'total_sequences': metrics.get('total_sequences', None),
                'gc_content': metrics.get('gc_content', None),
                'duplication_rate': metrics.get('duplication_rate', None),
                'sequence_length': metrics.get('sequence_length', None),
                'trimmed_file': trimmed_files[i] if i < len(trimmed_files) else None
            }
            # Try to get trimmed read count if possible
            trimmed_count = None
            if stat['trimmed_file'] and os.path.exists(stat['trimmed_file']):
                try:
                    with open(stat['trimmed_file'], 'rt') as f:
                        # Count reads in FASTQ (4 lines per read)
                        trimmed_count = sum(1 for _ in f) // 4
                except Exception:
                    trimmed_count = None
            stat['trimmed_sequences'] = trimmed_count
            qc_stats.append(stat)

        # --- Plot: Input vs Cleaned Reads ---
        plot_html = ""
        try:
            import matplotlib.pyplot as plt
            labels = []
            input_counts = []
            cleaned_counts = []
            for stat in qc_stats:
                labels.append(os.path.basename(stat['file']))
                input_counts.append(stat['total_sequences'] or 0)
                cleaned_counts.append(stat['trimmed_sequences'] or 0)
            fig, ax = plt.subplots(figsize=(6, 3.5))
            x = range(len(labels))
            ax.bar(x, input_counts, width=0.4, label='Input Reads', color='#007acc', alpha=0.7)
            ax.bar([i+0.4 for i in x], cleaned_counts, width=0.4, label='Cleaned Reads', color='#28a745', alpha=0.7)
            ax.set_xticks([i+0.2 for i in x])
            ax.set_xticklabels(labels, rotation=20)
            ax.set_ylabel('Read Count')
            ax.set_title('Input vs Cleaned Reads')
            ax.legend()
            plt.tight_layout()
            buf = io.BytesIO()
            plt.savefig(buf, format='png')
            plt.close(fig)
            buf.seek(0)
            img_base64 = base64.b64encode(buf.read()).decode('utf-8')
            plot_html = f"<img class='plot-img' src='data:image/png;base64,{img_base64}' alt='Input vs Cleaned Reads'/>"
        except Exception:
            plot_html = ""

        # --- FastQC Results Section ---
        html += "<div class='section'><h2>FastQC Results</h2>"
        if fastqc_results.get("success"):
            for stat in qc_stats:
                html += f"<div class='metric'><strong>{os.path.basename(stat['file'])}</strong>"
                if paired:
                    html += " (Paired-end)"
                html += "</div>"
                html += f"<div class='metric'>Input reads: {stat['total_sequences'] if stat['total_sequences'] is not None else 'N/A'}</div>"
                html += f"<div class='metric'>GC content: {stat['gc_content'] if stat['gc_content'] is not None else 'N/A'}%</div>"
                html += f"<div class='metric'>Duplication rate: {stat['duplication_rate'] if stat['duplication_rate'] is not None else 'N/A'}%</div>"
                html += f"<div class='metric'>Read length: {stat['sequence_length'] if stat['sequence_length'] is not None else 'N/A'}</div>"
        else:
            html += "<div class='error'>FastQC failed or not run.</div>"
        html += "</div>"

        # --- Trimming Results Section ---
        html += "<div class='section'><h2>Trimming Results (Trim Galore)</h2>"
        if trim_results.get("success"):
            html += plot_html
            for stat in qc_stats:
                html += f"<div class='metric'><strong>{os.path.basename(stat['file'])}</strong>"
                if paired:
                    html += " (Paired-end)"
                html += "</div>"
                html += f"<div class='metric'>Input reads: {stat['total_sequences'] if stat['total_sequences'] is not None else 'N/A'}</div>"
                html += f"<div class='metric'>Cleaned reads: {stat['trimmed_sequences'] if stat['trimmed_sequences'] is not None else 'N/A'}</div>"
                if stat['total_sequences'] and stat['trimmed_sequences']:
                    percent_retained = 100.0 * stat['trimmed_sequences'] / stat['total_sequences']
                    html += f"<div class='metric'>% Reads retained: {percent_retained:.1f}%</div>"
                    if percent_retained < 70.0:
                        html += f"<div class='warning'>Warning: Only {percent_retained:.1f}% reads retained after trimming. Check input quality and parameters.</div>"
                html += f"<div class='metric'>Trimmed file: {stat['trimmed_file'] if stat['trimmed_file'] else 'N/A'}</div>"
        else:
            html += f"<div class='error'>Trim Galore failed: {trim_results.get('error', 'Unknown error')}</div>"
        html += "</div>"

        # Warnings
        if not fastqc_results.get("success") or not trim_results.get("success"):
            html += "<div class='section warning'><strong>Warning:</strong> One or more QC steps failed. Please check the logs and input files.</div>"

        html += """
            </div>
        </body>
        </html>
        """
        return html 