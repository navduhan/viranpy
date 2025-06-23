#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Assembly and contig processing for viral genome annotation.
"""

import os
import subprocess
import shutil
from typing import Dict, Any, List, Optional
from pathlib import Path

from Bio import SeqIO
from ..utils.file_utils import cmd_exists


class Assembler:
    """
    Assembly and contig processing for viral genomes.
    
    This module handles:
    - SPAdes assembly
    - MEGAHIT assembly
    - CD-HIT clustering
    - Hybrid assembly creation
    """
    
    def __init__(self, config, logger=None):
        """Initialize the assembler."""
        self.config = config
        self.logger = logger
        self.assembly_results = {}
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        required_tools = ["spades.py", "megahit", "cd-hit-est"]
        missing_tools = [tool for tool in required_tools if not cmd_exists(tool)]
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        return True
    
    def run_spades(self, input_files: List[str], paired: bool = False, 
                   output_dir: str = "spades_assembly") -> Dict[str, Any]:
        """
        Run SPAdes assembly.
        
        Args:
            input_files: List of input FASTQ files
            paired: Whether the reads are paired-end
            output_dir: Output directory for SPAdes results
            
        Returns:
            Dictionary containing SPAdes results
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        # Debug logging
        self.logger.info(f"SPAdes input files: {input_files}")
        self.logger.info(f"SPAdes paired mode: {paired}")
        self.logger.info(f"SPAdes number of input files: {len(input_files)}")
        
        cmd = [
            "spades.py",
            "--meta",
            "-o", output_dir,
            "--threads", str(self.config.ncpus or 1),
            "--memory", str(getattr(self.config, 'spades_memory', 16))
        ]
        
        if paired:
            if len(input_files) >= 2:
                cmd.extend(["-1", input_files[0], "-2", input_files[1]])
                self.logger.info(f"SPAdes paired-end command: R1={input_files[0]}, R2={input_files[1]}")
            else:
                cmd.extend(["-s", input_files[0]])
                self.logger.warning(f"SPAdes paired mode but only 1 file provided: {input_files[0]}")
        else:
            cmd.extend(["-s", input_files[0]])
            self.logger.info(f"SPAdes single-end command: {input_files[0]}")
        
        self.logger.info(f"SPAdes full command: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"SPAdes assembly completed: {output_dir}")
            
            # Get contigs file
            contigs_file = Path(output_dir) / "contigs.fasta"
            if contigs_file.exists():
                return {"success": True, "contigs_file": str(contigs_file)}
            else:
                return {"success": False, "error": "No contigs file found"}
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"SPAdes assembly failed: {e}")
            return {"success": False, "error": str(e)}
    
    def run_megahit(self, input_files: List[str], paired: bool = False, 
                    output_dir: str = "megahit_assembly") -> Dict[str, Any]:
        """
        Run MEGAHIT assembly.
        
        Args:
            input_files: List of input FASTQ files
            paired: Whether the reads are paired-end
            output_dir: Output directory for MEGAHIT results
            
        Returns:
            Dictionary containing MEGAHIT results
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        # MEGAHIT does not allow pre-existing output dir
        if os.path.exists(output_dir):
            self.logger.warning(f"MEGAHIT output directory {output_dir} already exists. Deleting it to avoid MEGAHIT error.")
            shutil.rmtree(output_dir)
        # Do NOT pre-create output_dir for MEGAHIT
        
        # Debug logging
        self.logger.info(f"MEGAHIT input files: {input_files}")
        self.logger.info(f"MEGAHIT paired mode: {paired}")
        self.logger.info(f"MEGAHIT number of input files: {len(input_files)}")
        
        cmd = [
            "megahit",
            "-o", output_dir,
            "-t", str(self.config.ncpus or 1)
        ]
        
        if paired:
            if len(input_files) >= 2:
                cmd.extend(["-1", input_files[0], "-2", input_files[1]])
                self.logger.info(f"MEGAHIT paired-end command: R1={input_files[0]}, R2={input_files[1]}")
            else:
                cmd.extend(["-r", input_files[0]])
                self.logger.warning(f"MEGAHIT paired mode but only 1 file provided: {input_files[0]}")
        else:
            cmd.extend(["-r", input_files[0]])
            self.logger.info(f"MEGAHIT single-end command: {input_files[0]}")
        
        self.logger.info(f"MEGAHIT full command: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"MEGAHIT assembly completed: {output_dir}")
            
            # Get contigs file
            contigs_file = Path(output_dir) / "final.contigs.fa"
            if contigs_file.exists():
                return {"success": True, "contigs_file": str(contigs_file)}
            else:
                return {"success": False, "error": "No contigs file found"}
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"MEGAHIT assembly failed: {e}")
            return {"success": False, "error": str(e)}
    
    def run_cdhit(self, input_files: List[str], output_file: str = "cdhit_contigs.fasta",
                  identity: float = 0.95) -> Dict[str, Any]:
        """
        Run CD-HIT clustering to remove redundant contigs.
        
        Args:
            input_files: List of input FASTA files
            output_file: Output file for clustered contigs
            identity: Identity threshold for clustering
            
        Returns:
            Dictionary containing CD-HIT results
        """
        # Always write to output_file inside self.config.root_output
        output_file = os.path.join(self.config.root_output, output_file)
        
        # Combine input files if multiple
        combined_input = "combined_contigs.fasta"
        self._combine_fasta_files(input_files, combined_input)
        
        cmd = [
            "cd-hit-est",
            "-i", combined_input,
            "-o", output_file.replace('.fasta', ''),
            "-c", str(identity),
            "-T", str(self.config.ncpus or 1)
        ]
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"CD-HIT clustering completed: {output_file}")
            
            # Get statistics
            stats = self._get_cdhit_statistics(output_file)
            return {"success": True, "output_file": output_file, "statistics": stats}
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"CD-HIT clustering failed: {e}")
            return {"success": False, "error": str(e)}
        finally:
            # Clean up combined file
            if os.path.exists(combined_input):
                os.remove(combined_input)
    
    def create_hybrid_assembly(self, spades_contigs: str, megahit_contigs: str,
                              output_file: str = "hybrid_assembly.fasta") -> Dict[str, Any]:
        """
        Create hybrid assembly from SPAdes and MEGAHIT contigs.
        
        Args:
            spades_contigs: SPAdes contigs file
            megahit_contigs: MEGAHIT contigs file
            output_file: Output file for hybrid assembly
            
        Returns:
            Dictionary containing hybrid assembly results
        """
        # Always write to output_file inside self.config.root_output
        output_file = os.path.join(self.config.root_output, output_file)
        
        try:
            # Combine contigs from both assemblers
            all_contigs = []
            
            # Read SPAdes contigs
            if os.path.exists(spades_contigs):
                for record in SeqIO.parse(spades_contigs, "fasta"):
                    record.id = f"spades_{record.id}"
                    all_contigs.append(record)
            
            # Read MEGAHIT contigs
            if os.path.exists(megahit_contigs):
                for record in SeqIO.parse(megahit_contigs, "fasta"):
                    record.id = f"megahit_{record.id}"
                    all_contigs.append(record)
            
            # Write combined contigs
            SeqIO.write(all_contigs, "temp_combined.fasta", "fasta")
            
            # Run CD-HIT to remove redundancy
            cdhit_result = self.run_cdhit(["temp_combined.fasta"], output_file)
            
            # Clean up temporary file
            if os.path.exists("temp_combined.fasta"):
                os.remove("temp_combined.fasta")
            
            if cdhit_result["success"]:
                self.logger.info(f"Hybrid assembly created: {output_file}")
                return {"success": True, "output_file": output_file}
            else:
                return cdhit_result
                
        except Exception as e:
            self.logger.error(f"Hybrid assembly failed: {e}")
            return {"success": False, "error": str(e)}
    
    def filter_contigs_by_length(self, input_file: str, min_length: int = 200,
                                output_file: str = None) -> Dict[str, Any]:
        """
        Filter contigs by minimum length.
        
        Args:
            input_file: Input FASTA file
            min_length: Minimum contig length
            output_file: Output file (if None, auto-generated)
            
        Returns:
            Dictionary containing filtering results
        """
        # Always write to output_file inside self.config.root_output
        if output_file is None:
            output_file = os.path.join(self.config.root_output, input_file.replace('.fasta', f'_filtered_{min_length}bp.fasta'))
        else:
            output_file = os.path.join(self.config.root_output, output_file)
        
        try:
            filtered_contigs = []
            total_contigs = 0
            
            for record in SeqIO.parse(input_file, "fasta"):
                total_contigs += 1
                if len(record.seq) >= min_length:
                    filtered_contigs.append(record)
            
            SeqIO.write(filtered_contigs, output_file, "fasta")
            
            stats = {
                "total_contigs": total_contigs,
                "filtered_contigs": len(filtered_contigs),
                "total_bases": sum(len(contig.seq) for contig in filtered_contigs)
            }
            
            self.logger.info(f"Contig filtering completed: {len(filtered_contigs)}/{total_contigs} contigs retained")
            return {"success": True, "output_file": output_file, "statistics": stats}
            
        except Exception as e:
            self.logger.error(f"Contig filtering failed: {e}")
            return {"success": False, "error": str(e)}
    
    def _combine_fasta_files(self, input_files: List[str], output_file: str) -> None:
        """Combine multiple FASTA files into one."""
        with open(output_file, 'w') as outfile:
            for input_file in input_files:
                if os.path.exists(input_file):
                    with open(input_file, 'r') as infile:
                        shutil.copyfileobj(infile, outfile)
    
    def _get_cdhit_statistics(self, output_file: str) -> Dict[str, Any]:
        """Get statistics from CD-HIT output."""
        stats = {"total_contigs": 0, "total_bases": 0}
        
        try:
            for record in SeqIO.parse(output_file, "fasta"):
                stats["total_contigs"] += 1
                stats["total_bases"] += len(record.seq)
        except Exception as e:
            self.logger.warning(f"Could not parse CD-HIT output for statistics: {e}")
        
        return stats
    
    def generate_assembly_report(self, assembly_results: Dict[str, Any]) -> str:
        """Generate a comprehensive assembly report."""
        report_file = "assembly_report.html"
        
        html_content = self._create_assembly_html_report(assembly_results)
        
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        return report_file
    
    def _create_assembly_html_report(self, assembly_results: Dict[str, Any]) -> str:
        """Create HTML assembly report."""
        html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>Assembly Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; }
                .section { margin: 20px 0; padding: 10px; border: 1px solid #ddd; }
                .metric { margin: 10px 0; }
                .assembly-result { margin: 15px 0; padding: 10px; background: #f9f9f9; }
            </style>
        </head>
        <body>
            <h1>Assembly Report</h1>
        """
        
        for assembler, result in assembly_results.items():
            html += f"<div class='section'><h2>{assembler.upper()} Assembly</h2>"
            if result.get("success"):
                html += f"<div class='assembly-result'>"
                html += f"<p><strong>Status:</strong> Success</p>"
                html += f"<p><strong>Contigs file:</strong> {result.get('contigs_file', 'N/A')}</p>"
                
                if "statistics" in result:
                    stats = result["statistics"]
                    html += f"<p><strong>Total contigs:</strong> {stats.get('total_contigs', 'N/A')}</p>"
                    html += f"<p><strong>Total bases:</strong> {stats.get('total_bases', 'N/A'):,}</p>"
                
                html += "</div>"
            else:
                html += f"<div class='assembly-result'>"
                html += f"<p><strong>Status:</strong> Failed</p>"
                html += f"<p><strong>Error:</strong> {result.get('error', 'Unknown error')}</p>"
                html += "</div>"
            html += "</div>"
        
        html += "</body></html>"
        return html 