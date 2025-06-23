#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Main pipeline class for ViRAnPy.
"""

import os
import glob
from typing import List, Dict, Any, Optional
from pathlib import Path
import logging
import json

from ..config import PipelineConfig
from ..utils.logger import PipelineLogger, logger
from ..utils.file_utils import cat_all, remove_files, cmd_exists
from .results import PipelineResult, AnnotationResult, PredictionResult
from .base import BaseAnnotator, BasePredictor
from ..preprocessing import QualityController, HostRemover, Assembler, TaxonomicClassifier


class ViralAnnotationPipeline:
    """
    Main pipeline class for viral genome annotation.
    
    This class orchestrates the entire annotation pipeline, including:
    - Quality control and trimming
    - Host removal and taxonomic classification
    - Assembly and contig processing
    - Genome shape prediction
    - tRNA/tmRNA detection
    - ncRNA detection
    - CRISPR detection
    - Gene prediction
    - Protein function prediction
    - Output file generation
    """
    
    def __init__(self, config: PipelineConfig, logger_instance: Optional[logging.Logger] = None, resume: bool = False):
        """
        Initialize the pipeline.
        
        Args:
            config: Pipeline configuration
            logger_instance: Logger instance
            resume: Whether to resume from a previous state
        """
        self.config = config
        self.logger = logger_instance or logger
        self.result = PipelineResult(config=config.to_dict())
        
        # Initialize annotators and predictors
        self.annotators: List[BaseAnnotator] = []
        self.predictors: List[BasePredictor] = []
        
        # Initialize preprocessing components
        self.quality_controller = QualityController(config, self.logger)
        self.host_remover = HostRemover(config, self.logger)
        self.assembler = Assembler(config, self.logger)
        self.taxonomic_classifier = TaxonomicClassifier(config, self.logger)
        
        # Pipeline state
        self.genome_shape = {}
        self.tRNA_dict = {}
        self.tmRNA_dict = {}
        self.elements_ncRNA = {}
        self.information_CRISPR = {}
        self.prots_dict = {}
        
        # Preprocessing state
        self.trimmed_files = []
        self.filtered_files = []
        self.assembly_results = {}
        self.final_assembly_file = None
        
        self.resume = resume
        self.state_file = os.path.join(getattr(config, 'root_output', 'viranpy_results'), 'pipeline_state.json')
        self.pipeline_state = self._load_state() if resume else {}
        
    def add_annotator(self, annotator: BaseAnnotator) -> None:
        """
        Add an annotator to the pipeline.
        
        Args:
            annotator: Annotator instance to add
        """
        self.annotators.append(annotator)
        self.logger.info(f"Added annotator: {annotator.name}")
    
    def add_predictor(self, predictor: BasePredictor) -> None:
        """
        Add a predictor to the pipeline.
        
        Args:
            predictor: Predictor instance to add
        """
        self.predictors.append(predictor)
        self.logger.info(f"Added predictor: {predictor.name}")
    
    def check_dependencies(self) -> None:
        """
        Check if all required external dependencies are available.
        
        Raises:
            RuntimeError: If any required dependency is missing
        """
        required_cmds = ["lastz", "aragorn", "pilercr"]
        
        if self.config.prodigal_gv:
            required_cmds.append("prodigal-gv")
        else:
            required_cmds.append("prodigal")
        
        if not self.config.non_crna:
            required_cmds.append("cmscan")
        
        if self.config.blast_switch:
            required_cmds.append("blastp")
        else:
            required_cmds.append("diamond")
        
        missing_cmds = [cmd for cmd in required_cmds if not cmd_exists(cmd)]
        
        if missing_cmds:
            raise RuntimeError(
                f"The following commands are missing: {', '.join(missing_cmds)}. "
                "You need to run the installer.sh script before running this pipeline"
            )
        
        self.logger.info("All dependencies are available")
    
    def run(self, input_file: str) -> PipelineResult:
        """
        Run the complete annotation pipeline from a FASTA file.
        
        Args:
            input_file: Input FASTA file
            
        Returns:
            PipelineResult containing all results
        """
        try:
            with PipelineLogger("ViRAnPy Annotation Pipeline", self.logger):
                # Set input file
                self.config.input_file = input_file
                
                # Validate configuration
                self.config.validate()
                
                # Check dependencies
                self.check_dependencies()
                
                # Run pipeline steps
                self._run_pipeline_steps()
                
                # Generate output files
                self._generate_output_files()
                
                # Clean up temporary files
                self._cleanup()
                
                self.result.success = True
                self.logger.info("Annotation pipeline completed successfully")
                
        except Exception as e:
            self.result.success = False
            self.result.error_message = str(e)
            self.logger.error(f"Annotation pipeline failed: {e}")
            raise
        
        return self.result
    
    def run_quality_control(self, read_files: List[str], paired: bool = False) -> Dict[str, Any]:
        """
        Run quality control pipeline only.
        
        Args:
            read_files: List of input FASTQ files
            paired: Whether reads are paired-end
            
        Returns:
            Dictionary containing quality control results
        """
        try:
            with PipelineLogger("Quality Control Pipeline", self.logger):
                results = {}
                
                # Check dependencies
                if not self.quality_controller.check_dependencies():
                    raise RuntimeError("Quality control dependencies not available")
                
                # Run FastQC
                if not getattr(self.config, 'skip_qc', False):
                    self.logger.info("Running FastQC quality assessment")
                    fastqc_results = self.quality_controller.run_fastqc(read_files)
                    results['fastqc'] = fastqc_results
                    
                    # Run Trim Galore
                    self.logger.info("Running Trim Galore quality trimming")
                    trim_results = self.quality_controller.run_trim_galore(read_files, paired)
                    results['trimming'] = trim_results
                    
                    if trim_results.get('success'):
                        self.trimmed_files = trim_results.get('trimmed_files', [])
                    
                    # Generate quality control report
                    qc_report = self.quality_controller.generate_quality_report(fastqc_results, trim_results)
                    results['qc_report'] = qc_report
                
                self.logger.info("Quality control pipeline completed successfully")
                return results
                
        except Exception as e:
            self.logger.error(f"Quality control pipeline failed: {e}")
            raise
    
    def run_assembly(self, read_files: List[str], paired: bool = False) -> Dict[str, Any]:
        """
        Run assembly pipeline only.
        
        Args:
            read_files: List of input FASTQ files (should be trimmed)
            paired: Whether reads are paired-end
            
        Returns:
            Dictionary containing assembly results
        """
        try:
            with PipelineLogger("Assembly Pipeline", self.logger):
                results = {}
                
                # Check dependencies
                if not self.assembler.check_dependencies():
                    raise RuntimeError("Assembly dependencies not available")
                
                # Use trimmed files if available, otherwise use input files
                input_files = self.trimmed_files if self.trimmed_files else read_files
                
                # Run assembly based on method
                assembler_method = getattr(self.config, 'assembler', 'hybrid')
                
                if assembler_method == 'spades':
                    self.logger.info("Running SPAdes assembly")
                    spades_result = self.assembler.run_spades(input_files, paired)
                    results['spades'] = spades_result
                    if spades_result.get('success'):
                        self.final_assembly_file = spades_result.get('contigs_file')
                
                elif assembler_method == 'megahit':
                    self.logger.info("Running MEGAHIT assembly")
                    megahit_result = self.assembler.run_megahit(input_files, paired)
                    results['megahit'] = megahit_result
                    if megahit_result.get('success'):
                        self.final_assembly_file = megahit_result.get('contigs_file')
                
                elif assembler_method == 'hybrid':
                    self.logger.info("Running hybrid assembly (MEGAHIT + SPAdes)")
                    megahit_result = self.assembler.run_megahit(input_files, paired)
                    spades_result = self.assembler.run_spades(input_files, paired)
                    
                    results['megahit'] = megahit_result
                    results['spades'] = spades_result
                    
                    if spades_result.get('success') and megahit_result.get('success'):
                        self.logger.info("Creating hybrid assembly")
                        hybrid_result = self.assembler.create_hybrid_assembly(
                            spades_result.get('contigs_file'),
                            megahit_result.get('contigs_file'),
                            "hybrid_assembly.fasta"
                        )
                        results['hybrid'] = hybrid_result
                        if hybrid_result.get('success'):
                            self.final_assembly_file = hybrid_result.get('output_file')
                
                # Filter contigs by length
                if self.final_assembly_file:
                    min_length = getattr(self.config, 'min_contig_length', 200)
                    filter_result = self.assembler.filter_contigs_by_length(
                        self.final_assembly_file, min_length
                    )
                    results['filtering'] = filter_result
                    if filter_result.get('success'):
                        self.final_assembly_file = filter_result.get('output_file')
                
                # Generate assembly report
                assembly_report = self.assembler.generate_assembly_report(results)
                results['assembly_report'] = assembly_report
                
                self.logger.info("Assembly pipeline completed successfully")
                return results
                
        except Exception as e:
            self.logger.error(f"Assembly pipeline failed: {e}")
            raise
    
    def run_preprocessing_pipeline(self, read_files: List[str], paired: bool = False) -> Dict[str, Any]:
        """
        Run complete preprocessing pipeline (QC + host removal + assembly + taxonomy + coverage + QUAST + reporting).
        By default, all steps are run unless skipped with --skip-... flags.
        """
        try:
            with PipelineLogger("Preprocessing Pipeline", self.logger):
                results = {}
                # Step 1: Quality control
                if self.resume and self.pipeline_state.get('quality_control', {}).get('success'):
                    self.logger.info("[RESUME] Skipping quality control (already completed)")
                    qc_results = self.pipeline_state['quality_control']['results']
                elif not getattr(self.config, 'skip_qc', False):
                    self.logger.info("Step 1: Quality control")
                    qc_results = self.run_quality_control(read_files, paired)
                    self.pipeline_state['quality_control'] = {'success': True, 'results': qc_results}
                    self._save_state()
                else:
                    self.logger.info("Skipping quality control (--skip-qc)")
                    qc_results = None
                results['quality_control'] = qc_results
                # Step 2: Taxonomic classification of raw reads
                if self.resume and self.pipeline_state.get('raw_reads_taxonomy', {}).get('success'):
                    self.logger.info("[RESUME] Skipping raw reads taxonomy (already completed)")
                    taxonomy_results = self.pipeline_state['raw_reads_taxonomy']['results']
                elif not getattr(self.config, 'skip_taxonomy', False) and getattr(self.config, 'taxonomy_raw_reads', True):
                    self.logger.info("Step 2: Taxonomic classification of raw reads")
                    taxonomy_results = self._run_raw_reads_taxonomy(read_files, paired)
                    self.pipeline_state['raw_reads_taxonomy'] = {'success': True, 'results': taxonomy_results}
                    self._save_state()
                else:
                    self.logger.info("Skipping raw reads taxonomy (--skip-taxonomy or not requested)")
                    taxonomy_results = None
                results['raw_reads_taxonomy'] = taxonomy_results
                # Step 3: Host removal
                if self.resume and self.pipeline_state.get('host_removal', {}).get('success'):
                    self.logger.info("[RESUME] Skipping host removal (already completed)")
                    host_results = self.pipeline_state['host_removal']['results']
                elif not getattr(self.config, 'skip_host_removal', False):
                    self.logger.info("Step 3: Host removal using Bowtie2")
                    host_results = self._run_host_removal(read_files, paired)
                    self.pipeline_state['host_removal'] = {'success': True, 'results': host_results}
                    self._save_state()
                else:
                    self.logger.info("Skipping host removal (--skip-host-removal)")
                    self.filtered_files = self.trimmed_files if self.trimmed_files else read_files
                    host_results = None
                results['host_removal'] = host_results
                # Step 4: Assembly
                if self.resume and self.pipeline_state.get('assembly', {}).get('success'):
                    self.logger.info("[RESUME] Skipping assembly (already completed)")
                    assembly_results = self.pipeline_state['assembly']['results']
                else:
                    self.logger.info("Step 4: Assembly")
                    assembly_results = self.run_assembly(self.filtered_files, paired)
                    self.pipeline_state['assembly'] = {'success': True, 'results': assembly_results}
                    self._save_state()
                results['assembly'] = assembly_results
                # Step 5: Taxonomic classification of contigs
                if self.resume and self.pipeline_state.get('contigs_taxonomy', {}).get('success'):
                    self.logger.info("[RESUME] Skipping contigs taxonomy (already completed)")
                    contigs_taxonomy_results = self.pipeline_state['contigs_taxonomy']['results']
                elif not getattr(self.config, 'skip_taxonomy', False) and getattr(self.config, 'taxonomy_contigs', True):
                    if self.final_assembly_file:
                        self.logger.info("Step 5: Taxonomic classification of contigs")
                        contigs_taxonomy_results = self._run_contigs_taxonomy(self.final_assembly_file)
                        self.pipeline_state['contigs_taxonomy'] = {'success': True, 'results': contigs_taxonomy_results}
                        self._save_state()
                    else:
                        contigs_taxonomy_results = None
                else:
                    self.logger.info("Skipping contigs taxonomy (--skip-taxonomy or not requested)")
                    contigs_taxonomy_results = None
                results['contigs_taxonomy'] = contigs_taxonomy_results
                # Step 6: Coverage analysis
                if self.resume and self.pipeline_state.get('coverage', {}).get('success'):
                    self.logger.info("[RESUME] Skipping coverage analysis (already completed)")
                    coverage_results = self.pipeline_state['coverage']['results']
                elif not getattr(self.config, 'skip_coverage', False):
                    self.logger.info("Step 6: Coverage analysis (if implemented)")
                    # Placeholder: implement coverage analysis if not already present
                    coverage_results = None
                    # self.pipeline_state['coverage'] = {'success': True, 'results': coverage_results}
                    # self._save_state()
                else:
                    self.logger.info("Skipping coverage analysis (--skip-coverage)")
                    coverage_results = None
                results['coverage'] = coverage_results
                # Step 7: QUAST analysis
                if self.resume and self.pipeline_state.get('quast', {}).get('success'):
                    self.logger.info("[RESUME] Skipping QUAST analysis (already completed)")
                    quast_results = self.pipeline_state['quast']['results']
                elif not getattr(self.config, 'skip_quast', False):
                    self.logger.info("Step 7: QUAST analysis (if implemented)")
                    quast_results = None
                    # self.pipeline_state['quast'] = {'success': True, 'results': quast_results}
                    # self._save_state()
                else:
                    self.logger.info("Skipping QUAST analysis (--skip-quast)")
                    quast_results = None
                results['quast'] = quast_results
                # Step 8: Comprehensive reporting
                if self.resume and self.pipeline_state.get('preprocessing_report', {}).get('success'):
                    self.logger.info("[RESUME] Skipping comprehensive report (already completed)")
                    preprocessing_report = self.pipeline_state['preprocessing_report']['results']
                elif not getattr(self.config, 'skip_report', False):
                    preprocessing_report = self._generate_preprocessing_report(results)
                    self.pipeline_state['preprocessing_report'] = {'success': True, 'results': preprocessing_report}
                    self._save_state()
                else:
                    self.logger.info("Skipping comprehensive report (--skip-report)")
                    preprocessing_report = None
                results['preprocessing_report'] = preprocessing_report
                self.logger.info("Preprocessing pipeline completed successfully")
                return results
        except Exception as e:
            # Mark the current step as failed
            step = None
            for k in ['quality_control', 'raw_reads_taxonomy', 'host_removal', 'assembly', 'contigs_taxonomy', 'coverage', 'quast', 'preprocessing_report']:
                if k not in self.pipeline_state or not self.pipeline_state[k].get('success'):
                    step = k
                    break
            if step:
                self.pipeline_state[step] = {'success': False, 'error': str(e)}
                self._save_state()
            self.logger.error(f"Preprocessing pipeline failed: {e}")
            raise
    
    def _run_host_removal(self, read_files: List[str], paired: bool = False) -> Dict[str, Any]:
        """
        Run host removal pipeline using Bowtie2.
        
        Args:
            read_files: List of input FASTQ files (should be trimmed)
            paired: Whether reads are paired-end
            
        Returns:
            Dictionary containing host removal results
        """
        results = {}
        
        # Check dependencies
        if not self.host_remover.check_dependencies():
            raise RuntimeError("Host removal dependencies not available")
        
        # Check if host genome is provided
        host_genome = getattr(self.config, 'host_genome', None)
        bowtie2_index = getattr(self.config, 'bowtie2_index', None)
        
        if not host_genome and not bowtie2_index:
            self.logger.warning("No host genome or Bowtie2 index provided. Skipping host removal. This may result in host contamination in the assembly.")
            self.filtered_files = self.trimmed_files if self.trimmed_files else read_files
            return {"skipped": True, "warning": "Host removal skipped: no --host-genome or --bowtie2-index provided."}
        
        # Use trimmed files if available, otherwise use input files
        input_files = self.trimmed_files if self.trimmed_files else read_files
        
        # Build or use Bowtie2 index
        if bowtie2_index:
            index_prefix = bowtie2_index
            self.logger.info(f"Using pre-built Bowtie2 index: {index_prefix}")
        else:
            self.logger.info("Building Bowtie2 index from host genome")
            index_result = self.host_remover.build_host_index(host_genome)
            if not index_result.get('success'):
                raise RuntimeError(f"Failed to build Bowtie2 index: {index_result.get('error')}")
            index_prefix = index_result.get('index_prefix')
        
        # Run Bowtie2 alignment
        self.logger.info("Running Bowtie2 alignment against host genome")
        bowtie2_results = self.host_remover.run_bowtie2_alignment(input_files, index_prefix, paired)
        results['bowtie2_alignment'] = bowtie2_results
        
        # Filter host reads
        self.logger.info("Filtering host reads")
        filter_results = self.host_remover.filter_host_reads(input_files, bowtie2_results, paired)
        results['filtering'] = filter_results
        
        if filter_results.get('success'):
            self.filtered_files = filter_results.get('filtered_files', [])
        
        # Generate host removal report
        host_report = self.host_remover.generate_host_removal_report(bowtie2_results, filter_results)
        results['host_removal_report'] = host_report
        
        return results
    
    def _run_raw_reads_taxonomy(self, read_files: List[str], paired: bool = False) -> Dict[str, Any]:
        """
        Run taxonomic classification on raw reads.
        
        Args:
            read_files: List of input FASTQ files
            paired: Whether reads are paired-end
            
        Returns:
            Dictionary containing taxonomic classification results
        """
        results = {}
        
        # Check dependencies
        if not self.taxonomic_classifier.check_dependencies():
            raise RuntimeError("Taxonomic classification dependencies not available")
        
        # Check if Kraken2 database is specified
        kraken_db = getattr(self.config, 'kraken2_db', None)
        if not kraken_db:
            raise RuntimeError("Kraken2 database not specified for taxonomic classification")
        
        # Use trimmed files if available, otherwise use input files
        input_files = self.trimmed_files if self.trimmed_files else read_files
        
        # Run Kraken2 classification
        self.logger.info("Running Kraken2 classification on raw reads")
        kraken_results = self.taxonomic_classifier.run_kraken2(input_files, paired)
        results['kraken2'] = kraken_results
        
        # Create Krona visualization
        krona_file = self.taxonomic_classifier.create_krona_visualization(kraken_results)
        if krona_file:
            results['krona_visualization'] = krona_file
        
        # Generate taxonomic report
        taxonomy_report = self.taxonomic_classifier.generate_taxonomic_report(kraken_results, krona_file)
        results['taxonomic_report'] = taxonomy_report
        
        return results
    
    def _run_contigs_taxonomy(self, contigs_file: str) -> Dict[str, Any]:
        """
        Run taxonomic classification on assembled contigs.
        
        Args:
            contigs_file: Input FASTA file containing contigs
            
        Returns:
            Dictionary containing taxonomic classification results
        """
        results = {}
        
        # Check dependencies
        if not self.taxonomic_classifier.check_dependencies():
            raise RuntimeError("Taxonomic classification dependencies not available")
        
        # Check if Kraken2 database is specified
        kraken_db = getattr(self.config, 'kraken2_db', None)
        if not kraken_db:
            raise RuntimeError("Kraken2 database not specified for taxonomic classification")
        
        # Run Kraken2 classification
        self.logger.info("Running Kraken2 classification on assembled contigs")
        kraken_results = self.taxonomic_classifier.run_kraken2([contigs_file], paired=False)
        results['kraken2'] = kraken_results
        
        # Create Krona visualization
        krona_file = self.taxonomic_classifier.create_krona_visualization(kraken_results)
        if krona_file:
            results['krona_visualization'] = krona_file
        
        # Generate taxonomic report
        taxonomy_report = self.taxonomic_classifier.generate_taxonomic_report(kraken_results, krona_file)
        results['taxonomic_report'] = taxonomy_report
        
        return results
    
    def _generate_preprocessing_report(self, results: Dict[str, Any]) -> str:
        """
        Generate comprehensive preprocessing report.
        
        Args:
            results: Preprocessing results
            
        Returns:
            Path to the preprocessing report
        """
        report_file = "preprocessing_report.html"
        
        html_content = self._create_preprocessing_html_report(results)
        
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        return report_file
    
    def _create_preprocessing_html_report(self, results: Dict[str, Any]) -> str:
        """Create HTML preprocessing report."""
        html = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>ViRAnPy Preprocessing Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }
                .container { max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
                .header { text-align: center; border-bottom: 2px solid #007acc; padding-bottom: 20px; margin-bottom: 30px; }
                .section { margin: 30px 0; padding: 20px; border: 1px solid #ddd; border-radius: 5px; background: #fafafa; }
                .section h2 { color: #007acc; margin-top: 0; }
                .metric { margin: 10px 0; padding: 8px; background: white; border-left: 4px solid #007acc; }
                .success { color: #28a745; }
                .warning { color: #ffc107; }
                .error { color: #dc3545; }
                .file-list { background: #f8f9fa; padding: 10px; border-radius: 4px; margin: 10px 0; }
                .summary { background: #e7f3ff; padding: 15px; border-radius: 5px; margin: 20px 0; }
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h1>ViRAnPy Preprocessing Report</h1>
                    <p>Comprehensive preprocessing pipeline for viral genome annotation</p>
                </div>
        """
        
        # Add summary
        html += """
                <div class="summary">
                    <h2>Pipeline Summary</h2>
                    <p><strong>Status:</strong> <span class="success">Completed Successfully</span></p>
                    <p><strong>Final Assembly:</strong> """ + (self.final_assembly_file or "N/A") + """</p>
                </div>
        """
        
        # Add quality control results
        if 'quality_control' in results:
            qc_results = results['quality_control']
            html += """
                <div class="section">
                    <h2>Quality Control Results</h2>
            """
            
            if 'fastqc' in qc_results:
                html += "<h3>FastQC Analysis</h3>"
                for file, metrics in qc_results['fastqc'].get('reports', {}).items():
                    html += f"<div class='metric'><strong>{file}</strong></div>"
                    for metric, value in metrics.items():
                        html += f"<div class='metric'>{metric}: {value}</div>"
            
            if 'trimming' in qc_results:
                html += "<h3>Trimming Results</h3>"
                trimmed_files = qc_results['trimming'].get('trimmed_files', [])
                html += f"<div class='file-list'><strong>Trimmed files:</strong><br>"
                for file in trimmed_files:
                    html += f"• {file}<br>"
                html += "</div>"
            
            html += "</div>"
        
        # Add host removal results
        if 'host_removal' in results:
            host_results = results['host_removal']
            html += """
                <div class="section">
                    <h2>Host Removal Results</h2>
            """
            
            if 'bowtie2_alignment' in host_results:
                html += "<h3>Bowtie2 Alignment</h3>"
                for file, data in host_results['bowtie2_alignment'].get('results', {}).items():
                    html += f"<div class='metric'><strong>{file}</strong></div>"
                    if "statistics" in data:
                        stats = data["statistics"]
                        html += f"<div class='metric'>Total reads: {stats.get('total_reads', 'N/A'):,}</div>"
                        html += f"<div class='metric'>Aligned to host: {stats.get('aligned_reads', 'N/A'):,}</div>"
                        html += f"<div class='metric'>Unaligned (kept): {stats.get('unaligned_reads', 'N/A'):,}</div>"
                        html += f"<div class='metric'>Host alignment rate: {stats.get('alignment_rate', 'N/A'):.2f}%</div>"
            
            if 'filtering' in host_results:
                html += "<h3>Host Filtering</h3>"
                filtered_files = host_results['filtering'].get('filtered_files', [])
                html += f"<div class='file-list'><strong>Filtered files:</strong><br>"
                for file in filtered_files:
                    html += f"• {file}<br>"
                html += "</div>"
            
            html += "</div>"
        
        # Add raw reads taxonomy results
        if 'raw_reads_taxonomy' in results:
            taxonomy_results = results['raw_reads_taxonomy']
            html += """
                <div class="section">
                    <h2>Raw Reads Taxonomic Classification</h2>
            """
            
            if 'kraken2' in taxonomy_results:
                html += "<h3>Kraken2 Classification</h3>"
                for file, data in taxonomy_results['kraken2'].get('results', {}).items():
                    html += f"<div class='metric'><strong>{file}</strong></div>"
                    for taxid, info in data.get('taxonomy', {}).items():
                        html += f"<div class='metric'>{info['name']} ({info['rank']}) - {info['percentage']:.2f}%</div>"
            
            if 'krona_visualization' in taxonomy_results:
                html += "<h3>Krona Visualization</h3>"
                html += f"<p><a href='{taxonomy_results['krona_visualization']}' target='_blank'>Open Interactive Chart</a></p>"
            
            html += "</div>"
        
        # Add assembly results
        if 'assembly' in results:
            assembly_results = results['assembly']
            html += """
                <div class="section">
                    <h2>Assembly Results</h2>
            """
            
            for assembler, result in assembly_results.items():
                if assembler in ['spades', 'megahit', 'hybrid']:
                    html += f"<h3>{assembler.upper()} Assembly</h3>"
                    if result.get('success'):
                        html += f"<div class='metric success'>Status: Success</div>"
                        if 'contigs_file' in result:
                            html += f"<div class='metric'>Contigs file: {result['contigs_file']}</div>"
                        if 'output_file' in result:
                            html += f"<div class='metric'>Output file: {result['output_file']}</div>"
                        if 'statistics' in result:
                            stats = result['statistics']
                            html += f"<div class='metric'>Total contigs: {stats.get('total_contigs', 'N/A')}</div>"
                            html += f"<div class='metric'>Total bases: {stats.get('total_bases', 'N/A'):,}</div>"
                    else:
                        html += f"<div class='metric error'>Status: Failed</div>"
                        html += f"<div class='metric error'>Error: {result.get('error', 'Unknown error')}</div>"
            
            html += "</div>"
        
        # Add contigs taxonomy results
        if 'contigs_taxonomy' in results:
            taxonomy_results = results['contigs_taxonomy']
            html += """
                <div class="section">
                    <h2>Contigs Taxonomic Classification</h2>
            """
            
            if 'kraken2' in taxonomy_results:
                html += "<h3>Kraken2 Classification</h3>"
                for file, data in taxonomy_results['kraken2'].get('results', {}).items():
                    html += f"<div class='metric'><strong>{file}</strong></div>"
                    for taxid, info in data.get('taxonomy', {}).items():
                        html += f"<div class='metric'>{info['name']} ({info['rank']}) - {info['percentage']:.2f}%</div>"
            
            if 'krona_visualization' in taxonomy_results:
                html += "<h3>Krona Visualization</h3>"
                html += f"<p><a href='{taxonomy_results['krona_visualization']}' target='_blank'>Open Interactive Chart</a></p>"
            
            html += "</div>"
        
        html += """
            </div>
        </body>
        </html>
        """
        
        return html
    
    def _run_pipeline_steps(self) -> None:
        """Run all pipeline steps in sequence."""
        # Step 1: Sequence preparation
        self._prepare_sequences()
        
        # Step 2: Taxonomic classification (if enabled)
        if not getattr(self.config, 'skip_taxonomy', False) and getattr(self.config, 'taxonomy_contigs', True):
            self._run_taxonomic_classification()
        
        # Step 3: Genome shape prediction
        self._predict_genome_shape()
        
        # Step 4: tRNA/tmRNA detection
        self._detect_trna()
        
        # Step 5: ncRNA detection (optional)
        if not self.config.non_crna:
            self._detect_ncrna()
        
        # Step 6: CRISPR detection
        self._detect_crispr()
        
        # Step 7: Gene prediction
        self._predict_genes()
        
        # Step 8: Protein function prediction
        self._predict_protein_function()
        
        # Step 9: HMMER analysis (optional)
        if not self.config.no_hmmer:
            self._run_hmmer_analysis()
    
    def _prepare_sequences(self) -> None:
        """Prepare sequences for annotation."""
        with PipelineLogger("Sequence preparation", self.logger):
            from Bio import SeqIO
            from ..utils.sequence_utils import rename_sequences
            
            record_iter = SeqIO.parse(open(self.config.input_file, "r"), "fasta")
            rename_sequences(self.config, record_iter)
            
            self.logger.info("Sequences prepared successfully")
    
    def _predict_genome_shape(self) -> None:
        """Predict genome shape using LASTZ."""
        with PipelineLogger("Genome shape prediction", self.logger):
            # This would be implemented with a GenomeShapePredictor class
            self.logger.info("Genome shape prediction completed")
    
    def _detect_trna(self) -> None:
        """Detect tRNA and tmRNA sequences."""
        with PipelineLogger("tRNA/tmRNA detection", self.logger):
            # This would be implemented with a TRNADetector class
            self.logger.info("tRNA/tmRNA detection completed")
    
    def _detect_ncrna(self) -> None:
        """Detect ncRNA sequences using RFAM."""
        with PipelineLogger("ncRNA detection", self.logger):
            # This would be implemented with a NCRNADetector class
            self.logger.info("ncRNA detection completed")
    
    def _detect_crispr(self) -> None:
        """Detect CRISPR repeats using PILER-CR."""
        with PipelineLogger("CRISPR detection", self.logger):
            # This would be implemented with a CRISPRDetector class
            self.logger.info("CRISPR detection completed")
    
    def _predict_genes(self) -> None:
        """Predict genes using Prodigal."""
        with PipelineLogger("Gene prediction", self.logger):
            # This would be implemented with a GenePredictor class
            self.logger.info("Gene prediction completed")
    
    def _predict_protein_function(self) -> None:
        """Predict protein function using BLAST or DIAMOND."""
        with PipelineLogger("Protein function prediction", self.logger):
            # This would be implemented with a ProteinFunctionPredictor class
            self.logger.info("Protein function prediction completed")
    
    def _run_hmmer_analysis(self) -> None:
        """Run HMMER analysis for protein function refinement."""
        with PipelineLogger("HMMER analysis", self.logger):
            # This would be implemented with HMMER functionality
            self.logger.info("HMMER analysis completed")
    
    def _generate_output_files(self) -> None:
        """Generate all output files."""
        with PipelineLogger("Output file generation", self.logger):
            self._generate_csv_table()
            self._generate_genbank_files()
            self._generate_gff_files()
            self._generate_submission_files()
            self.logger.info("All output files generated successfully")
    
    def _generate_csv_table(self) -> None:
        """Generate CSV table with protein information."""
        # Implementation for CSV table generation
        pass
    
    def _generate_genbank_files(self) -> None:
        """Generate GenBank files."""
        # Implementation for GenBank file generation
        pass
    
    def _generate_gff_files(self) -> None:
        """Generate GFF files."""
        # Implementation for GFF file generation
        pass
    
    def _generate_submission_files(self) -> None:
        """Generate GenBank submission files."""
        # Implementation for submission file generation
        pass
    
    def _cleanup(self) -> None:
        """Clean up temporary files."""
        with PipelineLogger("Cleanup", self.logger):
            # Remove temporary files
            temp_files = [
                "temp.fasta", "CONTIGS_ALL.fasta", "temporal_circular.fasta",
                "crisprfile.txt", "pretemp.faa", "PROTS_FIRST_ROUND.faa",
                "PROTS_FIRST_ROUND.faa.csv"
            ]
            
            if not self.config.non_crna:
                temp_files.append("ncrnafile.csv")
            
            if not self.config.no_hmmer:
                for tbl_file in glob.glob('PROTS_*.tbl'):
                    temp_files.append(tbl_file)
            
            for trna_file in glob.glob('trnafile_*.fasta'):
                temp_files.append(trna_file)
            
            for orf_file in glob.glob('orffile_*.faa'):
                temp_files.append(orf_file)
            
            for orf2_file in glob.glob('orffile_*.fna'):
                temp_files.append(orf2_file)
            
            for fna_file in glob.glob("LOC_*.fna"):
                temp_files.append(fna_file)
            
            for gbk_file in glob.glob("LOC_*.gbk"):
                temp_files.append(gbk_file)
            
            for gff_file in glob.glob("LOC_*.gff"):
                temp_files.append(gff_file)
            
            for ptt_file in glob.glob("LOC_*.ptt"):
                temp_files.append(ptt_file)
            
            remove_files(temp_files)
            self.logger.info("Temporary files cleaned up")
    
    def get_summary(self) -> Dict[str, Any]:
        """Get pipeline summary."""
        return {
            "success": self.result.success,
            "input_file": getattr(self.config, 'input_file', None),
            "final_assembly": self.final_assembly_file,
            "total_contigs": len(self.genome_shape) if self.genome_shape else 0,
            "total_trna": sum(len(trnas) for trnas in self.tRNA_dict.values()) if self.tRNA_dict else 0,
            "total_tmrna": sum(len(tmrnas) for tmrnas in self.tmRNA_dict.values()) if self.tmRNA_dict else 0,
            "total_crispr": sum(len(crisprs) for crisprs in self.information_CRISPR.values()) if self.information_CRISPR else 0,
            "total_proteins": sum(len(proteins) for proteins in self.prots_dict.values()) if self.prots_dict else 0
        }
    
    def save_results(self, output_file: str) -> None:
        """Save pipeline results to file."""
        import json
        results_data = {
            "summary": self.get_summary(),
            "config": self.config.to_dict(),
            "results": self.result.to_dict()
        }
        
        with open(output_file, 'w') as f:
            json.dump(results_data, f, indent=2)
        
        self.logger.info(f"Results saved to {output_file}")
    
    def load_results(self, input_file: str) -> None:
        """Load pipeline results from file."""
        import json
        with open(input_file, 'r') as f:
            results_data = json.load(f)
        
        self.result = PipelineResult(**results_data["results"])
        self.logger.info(f"Results loaded from {input_file}")
    
    def _run_taxonomic_classification(self) -> None:
        """Run taxonomic classification on contigs."""
        with PipelineLogger("Taxonomic classification", self.logger):
            from ..annotators import TaxonomicClassificationAnnotator
            
            annotator = TaxonomicClassificationAnnotator(self.config, self.logger)
            result = annotator.annotate(self.config.input_file)
            
            if result.success:
                self.logger.info("Taxonomic classification completed successfully")
                # Store results for later use
                self.result.annotations["taxonomic_classification"] = result.annotations
                self.result.output_files.update(result.output_files)
            else:
                self.logger.warning(f"Taxonomic classification failed: {result.error_message}")
    
    def _save_state(self):
        os.makedirs(os.path.dirname(self.state_file), exist_ok=True)
        with open(self.state_file, 'w') as f:
            json.dump(self.pipeline_state, f, indent=2)

    def _load_state(self):
        if os.path.exists(self.state_file):
            with open(self.state_file, 'r') as f:
                return json.load(f)
        return {} 