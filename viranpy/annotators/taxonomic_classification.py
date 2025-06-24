#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Taxonomic classification annotator for viral genome annotation.
"""

import os
import subprocess
from typing import Dict, Any, List, Optional
from pathlib import Path

from ..core.base import BaseAnnotator
from ..core.results import AnnotationResult
from ..utils.file_utils import cmd_exists


class TaxonomicClassificationAnnotator(BaseAnnotator):
    """
    Taxonomic classification annotator using Kraken2 and Krona.
    
    This annotator provides:
    - Kraken2 taxonomic classification of contigs
    - Krona visualization of taxonomic composition
    - Taxonomic analysis and reporting
    """
    
    def __init__(self, config, logger=None):
        """Initialize the taxonomic classification annotator."""
        super().__init__(config, logger)
        self.name = "TaxonomicClassification"
        self.description = "Taxonomic classification using Kraken2 and Krona"
        self.kraken_results = {}
        self.krona_file = None
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        required_tools = ["kraken2", "ktImportTaxonomy"]
        missing_tools = [tool for tool in required_tools if not cmd_exists(tool)]
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        return True
    
    def annotate(self, input_file: str) -> AnnotationResult:
        """
        Run taxonomic classification on contigs.
        
        Args:
            input_file: Input FASTA file containing contigs
            
        Returns:
            AnnotationResult containing taxonomic classification results
        """
        result = AnnotationResult(
            annotator_name=self.name,
            input_file=input_file,
            success=False
        )
        
        try:
            # Check dependencies
            if not self.check_dependencies():
                result.error_message = "Missing required dependencies"
                return result
            
            # Check if Kraken2 database is specified
            kraken_db = getattr(self.config, 'kraken2_db', None)
            if not kraken_db:
                result.error_message = "Kraken2 database not specified in configuration"
                return result
            
            # Run Kraken2 classification
            self.logger.info("Running Kraken2 taxonomic classification on contigs")
            kraken_results = self._run_kraken2_classification(input_file)
            
            if not kraken_results.get("success"):
                result.error_message = kraken_results.get("error", "Kraken2 classification failed")
                return result
            
            self.kraken_results = kraken_results
            
            # Create Krona visualization
            self.logger.info("Creating Krona visualization")
            self.krona_file = self._create_krona_visualization()
            
            # Prepare results
            result.success = True
            result.output_files = {
                "kraken2_report": kraken_results.get("report_file"),
                "kraken2_output": kraken_results.get("output_file"),
                "krona_visualization": self.krona_file,
            }
            result.annotations = kraken_results.get("taxonomy", {})
            result.metadata = {
                "total_contigs": kraken_results.get("total_contigs", 0),
                "classified_contigs": kraken_results.get("classified_contigs", 0),
                "classification_rate": kraken_results.get("classification_rate", 0.0)
            }
            
            self.logger.info("Taxonomic classification completed successfully")
            
        except Exception as e:
            result.error_message = str(e)
            self.logger.error(f"Taxonomic classification failed: {e}")
        
        return result
    
    def _run_kraken2_classification(self, input_file: str) -> Dict[str, Any]:
        """Run Kraken2 classification on contigs."""
        output_dir = "kraken2_contigs_results"
        os.makedirs(output_dir, exist_ok=True)
        
        base_name = Path(input_file).stem
        report_file = Path(output_dir) / f"{base_name}_kraken2.report"
        output_file = Path(output_dir) / f"{base_name}_kraken2.out"
        
        kraken_db = getattr(self.config, 'kraken2_db', None)
        
        cmd = [
            "kraken2",
            "--db", kraken_db,
            "--output", str(output_file),
            "--report", str(report_file),
            "--threads", str(getattr(self.config, 'ncpus', 1)),
            input_file
        ]
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"Kraken2 classification completed: {input_file}")
            
            # Parse results
            taxonomy = self._parse_kraken2_report(report_file)
            stats = self._calculate_classification_stats(report_file)
            
            return {
                "success": True,
                "report_file": str(report_file),
                "output_file": str(output_file),
                "taxonomy": taxonomy,
                "total_contigs": stats["total_contigs"],
                "classified_contigs": stats["classified_contigs"],
                "classification_rate": stats["classification_rate"]
            }
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Kraken2 classification failed: {e}")
            return {"success": False, "error": str(e)}
    
    def _create_krona_visualization(self) -> Optional[str]:
        """Create Krona visualization of taxonomic classification."""
        output_dir = "kraken2_contigs_results"
        kraken_files = []
        
        for file_path in Path(output_dir).glob("*_kraken2.out"):
            kraken_files.append(str(file_path))
        
        if not kraken_files:
            self.logger.warning("No Kraken2 output files found for Krona visualization")
            return None
        
        krona_html = Path(output_dir) / "contigs_taxonomy_krona.html"
        
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
        taxonomy = {}
        
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
                        taxonomy[taxid] = {
                            "name": name,
                            "rank": rank,
                            "percentage": percentage,
                            "reads": reads,
                            "tax_reads": tax_reads
                        }
        
        return taxonomy
    
    def _calculate_classification_stats(self, report_file: Path) -> Dict[str, Any]:
        """Calculate classification statistics from Kraken2 report."""
        stats = {
            "total_contigs": 0,
            "classified_contigs": 0,
            "classification_rate": 0.0
        }
        
        try:
            with open(report_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        reads = int(parts[1])
                        rank = parts[3]
                        
                        if rank == "U":  # Unclassified
                            stats["total_contigs"] += reads
                        else:
                            stats["classified_contigs"] += reads
                            stats["total_contigs"] += reads
            
            if stats["total_contigs"] > 0:
                stats["classification_rate"] = (stats["classified_contigs"] / stats["total_contigs"]) * 100
                
        except Exception as e:
            self.logger.warning(f"Could not calculate classification statistics: {e}")
        
        return stats 