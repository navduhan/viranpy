#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Comprehensive reporting for viral metagenomic analysis.
"""

import os
import json
from datetime import datetime
from typing import Dict, Any, List, Optional
from pathlib import Path
import pandas as pd

from .coverage_calculator import CoverageCalculator
from .assembly_stats import AssemblyStatsCalculator


class ComprehensiveReporter:
    """
    Generate comprehensive reports for viral metagenomic analysis.
    
    This module creates HTML reports with:
    - Pipeline summary and statistics
    - Quality control results
    - Host removal statistics
    - Assembly statistics (QUAST)
    - Coverage analysis
    - Taxonomic classification results
    - Annotation results
    """
    
    def __init__(self, config, logger=None):
        """Initialize the comprehensive reporter."""
        self.config = config
        self.logger = logger
        self.report_data = {}
    
    def generate_comprehensive_report(self,
                                    output_dir: str,
                                    pipeline_results: Dict[str, Any],
                                    sample_name: str) -> str:
        """
        Generate comprehensive HTML report for viral metagenomic analysis.
        
        Args:
            output_dir: Output directory for reports
            pipeline_results: Results from the pipeline
            sample_name: Name of the sample
            
        Returns:
            Path to the generated HTML report
        """
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            # Collect all data for the report
            self._collect_report_data(pipeline_results)
            
            # Generate HTML report
            html_report = self._generate_html_report(output_dir, sample_name)
            
            # Generate TSV files
            self._generate_tsv_files(output_dir, sample_name)
            
            # Generate annotated FASTA
            self._generate_annotated_fasta(output_dir, sample_name)
            
            return html_report
            
        except Exception as e:
            self.logger.error(f"Report generation failed: {e}")
            return ""
    
    def _collect_report_data(self, pipeline_results: Dict[str, Any]):
        """Collect all data needed for the comprehensive report."""
        self.report_data = {
            "pipeline_summary": self._extract_pipeline_summary(pipeline_results),
            "quality_control": self._extract_qc_data(pipeline_results),
            "host_removal": self._extract_host_removal_data(pipeline_results),
            "assembly": self._extract_assembly_data(pipeline_results),
            "coverage": self._extract_coverage_data(pipeline_results),
            "taxonomy": self._extract_taxonomy_data(pipeline_results),
            "annotation": self._extract_annotation_data(pipeline_results)
        }
    
    def _extract_pipeline_summary(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract pipeline summary information."""
        return {
            "sample_name": results.get("sample_name", "Unknown"),
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "pipeline_version": "ViRAnPy 0.1.0",
            "total_runtime": results.get("runtime", "Unknown"),
            "success": results.get("success", False)
        }
    
    def _extract_qc_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract quality control data."""
        qc_results = results.get("quality_control", {})
        return {
            "raw_reads": qc_results.get("raw_reads", 0),
            "cleaned_reads": qc_results.get("cleaned_reads", 0),
            "trimmed_reads": qc_results.get("trimmed_reads", 0),
            "quality_threshold": qc_results.get("quality_threshold", 20),
            "min_length": qc_results.get("min_length", 20),
            "fastqc_reports": qc_results.get("fastqc_reports", [])
        }
    
    def _extract_host_removal_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract host removal data."""
        host_results = results.get("host_removal", {})
        return {
            "input_reads": host_results.get("input_reads", 0),
            "host_reads": host_results.get("host_reads", 0),
            "viral_reads": host_results.get("viral_reads", 0),
            "host_percentage": host_results.get("host_percentage", 0),
            "viral_percentage": host_results.get("viral_percentage", 0),
            "bowtie2_stats": host_results.get("bowtie2_stats", {})
        }
    
    def _extract_assembly_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract assembly data."""
        assembly_results = results.get("assembly", {})
        return {
            "assembler": assembly_results.get("assembler", "Unknown"),
            "total_contigs": assembly_results.get("total_contigs", 0),
            "total_bases": assembly_results.get("total_bases", 0),
            "n50": assembly_results.get("n50", 0),
            "largest_contig": assembly_results.get("largest_contig", 0),
            "mean_contig_length": assembly_results.get("mean_contig_length", 0),
            "gc_content": assembly_results.get("gc_content", 0),
            "quast_results": assembly_results.get("quast_results", {})
        }
    
    def _extract_coverage_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract coverage data."""
        coverage_results = results.get("coverage", {})
        return {
            "coverage_stats": coverage_results.get("coverage_stats", {}),
            "coverage_summary": coverage_results.get("coverage_summary", {}),
            "coverage_report": coverage_results.get("coverage_report", "")
        }
    
    def _extract_taxonomy_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract taxonomic classification data."""
        taxonomy_results = results.get("taxonomy", {})
        return {
            "kraken2_results": taxonomy_results.get("kraken2_results", {}),
            "taxonomic_classifications": taxonomy_results.get("classifications", {}),
            "viral_taxa": taxonomy_results.get("viral_taxa", []),
            "confidence_scores": taxonomy_results.get("confidence_scores", {})
        }
    
    def _extract_annotation_data(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Extract annotation data."""
        annotation_results = results.get("annotation", {})
        return {
            "total_genes": annotation_results.get("total_genes", 0),
            "predicted_proteins": annotation_results.get("predicted_proteins", 0),
            "functional_annotations": annotation_results.get("functional_annotations", {}),
            "genbank_file": annotation_results.get("genbank_file", ""),
            "gff_file": annotation_results.get("gff_file", ""),
            "csv_file": annotation_results.get("csv_file", "")
        }
    
    def _generate_html_report(self, output_dir: str, sample_name: str) -> str:
        """Generate comprehensive HTML report."""
        html_file = Path(output_dir) / f"{sample_name}_comprehensive_report.html"
        
        html_content = self._create_html_template()
        
        with open(html_file, 'w') as f:
            f.write(html_content)
        
        return str(html_file)
    
    def _create_html_template(self) -> str:
        """Create HTML template for the comprehensive report."""
        return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ViRAnPy Viral Metagenomic Analysis Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }}
        .header {{
            text-align: center;
            border-bottom: 3px solid #2c3e50;
            padding-bottom: 20px;
            margin-bottom: 30px;
        }}
        .header h1 {{
            color: #2c3e50;
            margin: 0;
            font-size: 2.5em;
        }}
        .header p {{
            color: #7f8c8d;
            margin: 10px 0 0 0;
            font-size: 1.1em;
        }}
        .section {{
            margin-bottom: 40px;
            padding: 20px;
            border: 1px solid #ecf0f1;
            border-radius: 8px;
            background-color: #fafafa;
        }}
        .section h2 {{
            color: #34495e;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
            margin-top: 0;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #3498db;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        .stat-card h3 {{
            margin: 0 0 10px 0;
            color: #2c3e50;
            font-size: 1.1em;
        }}
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            color: #3498db;
        }}
        .stat-label {{
            color: #7f8c8d;
            font-size: 0.9em;
            margin-top: 5px;
        }}
        .table-container {{
            overflow-x: auto;
            margin: 20px 0;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ecf0f1;
        }}
        th {{
            background-color: #34495e;
            color: white;
            font-weight: 600;
        }}
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        .progress-bar {{
            width: 100%;
            height: 20px;
            background-color: #ecf0f1;
            border-radius: 10px;
            overflow: hidden;
            margin: 10px 0;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #3498db, #2ecc71);
            transition: width 0.3s ease;
        }}
        .warning {{
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            color: #856404;
            padding: 15px;
            border-radius: 5px;
            margin: 10px 0;
        }}
        .success {{
            background-color: #d4edda;
            border: 1px solid #c3e6cb;
            color: #155724;
            padding: 15px;
            border-radius: 5px;
            margin: 10px 0;
        }}
        .footer {{
            text-align: center;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
            color: #7f8c8d;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 ViRAnPy Viral Metagenomic Analysis</h1>
            <p>Comprehensive Analysis Report</p>
            <p><strong>Sample:</strong> {self.report_data['pipeline_summary'].get('sample_name', 'Unknown')} | 
               <strong>Date:</strong> {self.report_data['pipeline_summary'].get('analysis_date', 'Unknown')}</p>
        </div>

        {self._generate_pipeline_summary_section()}
        {self._generate_quality_control_section()}
        {self._generate_host_removal_section()}
        {self._generate_assembly_section()}
        {self._generate_coverage_section()}
        {self._generate_taxonomy_section()}
        {self._generate_annotation_section()}

        <div class="footer">
            <p>Report generated by ViRAnPy Viral Metagenomic Pipeline v1.0.0</p>
            <p>For more information, visit: <a href="https://github.com/viranpy">https://github.com/viranpy</a></p>
        </div>
    </div>

    <script>
        // Add interactive features
        document.addEventListener('DOMContentLoaded', function() {{
            // Animate progress bars
            const progressBars = document.querySelectorAll('.progress-fill');
            progressBars.forEach(bar => {{
                const width = bar.style.width;
                bar.style.width = '0%';
                setTimeout(() => {{
                    bar.style.width = width;
                }}, 500);
            }});
        }});
    </script>
</body>
</html>
        """
    
    def _generate_pipeline_summary_section(self) -> str:
        """Generate pipeline summary section."""
        summary = self.report_data['pipeline_summary']
        return f"""
        <div class="section">
            <h2>📊 Pipeline Summary</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <h3>Sample Name</h3>
                    <div class="stat-value">{summary.get('sample_name', 'Unknown')}</div>
                </div>
                <div class="stat-card">
                    <h3>Analysis Date</h3>
                    <div class="stat-value">{summary.get('analysis_date', 'Unknown')}</div>
                </div>
                <div class="stat-card">
                    <h3>Pipeline Version</h3>
                    <div class="stat-value">{summary.get('pipeline_version', 'Unknown')}</div>
                </div>
                <div class="stat-card">
                    <h3>Status</h3>
                    <div class="stat-value">{'✅ Success' if summary.get('success') else '❌ Failed'}</div>
                </div>
            </div>
        </div>
        """
    
    def _generate_quality_control_section(self) -> str:
        """Generate quality control section."""
        qc = self.report_data['quality_control']
        raw_reads = qc.get('raw_reads', 0)
        cleaned_reads = qc.get('cleaned_reads', 0)
        retention_rate = (cleaned_reads / raw_reads * 100) if raw_reads > 0 else 0
        
        return f"""
        <div class="section">
            <h2>🔍 Quality Control</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <h3>Raw Reads</h3>
                    <div class="stat-value">{raw_reads:,}</div>
                </div>
                <div class="stat-card">
                    <h3>Cleaned Reads</h3>
                    <div class="stat-value">{cleaned_reads:,}</div>
                </div>
                <div class="stat-card">
                    <h3>Retention Rate</h3>
                    <div class="stat-value">{retention_rate:.1f}%</div>
                    <div class="progress-bar">
                        <div class="progress-fill" style="width: {retention_rate}%"></div>
                    </div>
                </div>
                <div class="stat-card">
                    <h3>Quality Threshold</h3>
                    <div class="stat-value">{qc.get('quality_threshold', 20)}</div>
                </div>
            </div>
        </div>
        """
    
    def _generate_host_removal_section(self) -> str:
        """Generate host removal section."""
        host = self.report_data['host_removal']
        input_reads = host.get('input_reads', 0)
        viral_reads = host.get('viral_reads', 0)
        viral_percentage = host.get('viral_percentage', 0)
        
        return f"""
        <div class="section">
            <h2>🧬 Host Removal</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <h3>Input Reads</h3>
                    <div class="stat-value">{input_reads:,}</div>
                </div>
                <div class="stat-card">
                    <h3>Viral Reads</h3>
                    <div class="stat-value">{viral_reads:,}</div>
                </div>
                <div class="stat-card">
                    <h3>Viral Percentage</h3>
                    <div class="stat-value">{viral_percentage:.1f}%</div>
                    <div class="progress-bar">
                        <div class="progress-fill" style="width: {viral_percentage}%"></div>
                    </div>
                </div>
                <div class="stat-card">
                    <h3>Host Reads Removed</h3>
                    <div class="stat-value">{host.get('host_reads', 0):,}</div>
                </div>
            </div>
        </div>
        """
    
    def _generate_assembly_section(self) -> str:
        """Generate assembly section."""
        assembly = self.report_data['assembly']
        
        return f"""
        <div class="section">
            <h2>🧩 Assembly Statistics</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <h3>Total Contigs</h3>
                    <div class="stat-value">{assembly.get('total_contigs', 0):,}</div>
                </div>
                <div class="stat-card">
                    <h3>Total Bases</h3>
                    <div class="stat-value">{assembly.get('total_bases', 0):,}</div>
                </div>
                <div class="stat-card">
                    <h3>N50</h3>
                    <div class="stat-value">{assembly.get('n50', 0):,}</div>
                </div>
                <div class="stat-card">
                    <h3>Largest Contig</h3>
                    <div class="stat-value">{assembly.get('largest_contig', 0):,}</div>
                </div>
                <div class="stat-card">
                    <h3>Mean Contig Length</h3>
                    <div class="stat-value">{assembly.get('mean_contig_length', 0):.0f}</div>
                </div>
                <div class="stat-card">
                    <h3>GC Content</h3>
                    <div class="stat-value">{assembly.get('gc_content', 0):.1f}%</div>
                </div>
            </div>
        </div>
        """
    
    def _generate_coverage_section(self) -> str:
        """Generate coverage section."""
        coverage = self.report_data['coverage']
        summary = coverage.get('coverage_summary', {})
        
        return f"""
        <div class="section">
            <h2>📈 Coverage Analysis</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <h3>Total Contigs</h3>
                    <div class="stat-value">{summary.get('total_contigs', 0):,}</div>
                </div>
                <div class="stat-card">
                    <h3>Mean Coverage</h3>
                    <div class="stat-value">{summary.get('mean_coverage_across_contigs', 0):.1f}x</div>
                </div>
                <div class="stat-card">
                    <h3>High Coverage Contigs</h3>
                    <div class="stat-value">{summary.get('high_coverage_contigs', 0):,}</div>
                </div>
                <div class="stat-card">
                    <h3>Well Covered Contigs</h3>
                    <div class="stat-value">{summary.get('well_covered_contigs', 0):,}</div>
                </div>
            </div>
        </div>
        """
    
    def _generate_taxonomy_section(self) -> str:
        """Generate taxonomy section."""
        taxonomy = self.report_data['taxonomy']
        viral_taxa = taxonomy.get('viral_taxa', [])
        
        taxa_table = ""
        if viral_taxa:
            taxa_table = """
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Contig ID</th>
                            <th>Phylum</th>
                            <th>Class</th>
                            <th>Order</th>
                            <th>Family</th>
                            <th>Genus</th>
                            <th>Species</th>
                            <th>Confidence (%)</th>
                        </tr>
                    </thead>
                    <tbody>
            """
            for taxon in viral_taxa[:10]:  # Show first 10
                taxa_table += f"""
                        <tr>
                            <td>{taxon.get('contig_id', 'N/A')}</td>
                            <td>{taxon.get('phylum', 'N/A')}</td>
                            <td>{taxon.get('class', 'N/A')}</td>
                            <td>{taxon.get('order', 'N/A')}</td>
                            <td>{taxon.get('family', 'N/A')}</td>
                            <td>{taxon.get('genus', 'N/A')}</td>
                            <td>{taxon.get('species', 'N/A')}</td>
                            <td>{taxon.get('confidence', 0):.1f}</td>
                        </tr>
                """
            taxa_table += """
                    </tbody>
                </table>
            </div>
            """
        
        return f"""
        <div class="section">
            <h2>🔬 Taxonomic Classification</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <h3>Viral Contigs Identified</h3>
                    <div class="stat-value">{len(viral_taxa)}</div>
                </div>
                <div class="stat-card">
                    <h3>Unique Viral Families</h3>
                    <div class="stat-value">{len(set(t.get('family') for t in viral_taxa if t.get('family') != 'N/A'))}</div>
                </div>
                <div class="stat-card">
                    <h3>High Confidence Classifications</h3>
                    <div class="stat-value">{len([t for t in viral_taxa if t.get('confidence', 0) >= 80])}</div>
                </div>
            </div>
            {taxa_table}
        </div>
        """
    
    def _generate_annotation_section(self) -> str:
        """Generate annotation section."""
        annotation = self.report_data['annotation']
        
        return f"""
        <div class="section">
            <h2>🧬 Gene Annotation</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <h3>Total Genes</h3>
                    <div class="stat-value">{annotation.get('total_genes', 0):,}</div>
                </div>
                <div class="stat-card">
                    <h3>Predicted Proteins</h3>
                    <div class="stat-value">{annotation.get('predicted_proteins', 0):,}</div>
                </div>
                <div class="stat-card">
                    <h3>Functional Annotations</h3>
                    <div class="stat-value">{len(annotation.get('functional_annotations', {}))}</div>
                </div>
            </div>
        </div>
        """
    
    def _generate_tsv_files(self, output_dir: str, sample_name: str):
        """Generate TSV files with detailed results."""
        # Taxonomic classification TSV
        taxonomy = self.report_data['taxonomy']
        viral_taxa = taxonomy.get('viral_taxa', [])
        coverage = self.report_data['coverage']
        coverage_stats = coverage.get('coverage_stats', {})
        
        if viral_taxa:
            tsv_file = Path(output_dir) / f"{sample_name}_taxonomy_coverage.tsv"
            with open(tsv_file, 'w') as f:
                f.write("Contig_ID\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tConfidence(%)\tMean_Coverage\tCoverage_Breadth(%)\n")
                for taxon in viral_taxa:
                    contig_id = taxon.get('contig_id', 'N/A')
                    coverage_info = coverage_stats.get(contig_id, {})
                    f.write(f"{contig_id}\t{taxon.get('phylum', 'N/A')}\t{taxon.get('class', 'N/A')}\t"
                           f"{taxon.get('order', 'N/A')}\t{taxon.get('family', 'N/A')}\t"
                           f"{taxon.get('genus', 'N/A')}\t{taxon.get('species', 'N/A')}\t"
                           f"{taxon.get('confidence', 0):.1f}\t"
                           f"{coverage_info.get('mean_coverage', 0):.2f}\t"
                           f"{coverage_info.get('coverage_breadth', 0):.2f}\n")
    
    def _generate_annotated_fasta(self, output_dir: str, sample_name: str):
        """Generate annotated FASTA file with taxonomic information."""
        # This would integrate with the annotation results to create
        # a FASTA file with annotated sequences
        pass 