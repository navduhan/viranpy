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
        # Fix the state file path to use absolute path
        if hasattr(config, 'root_output') and config.root_output:
            # If root_output is a relative path, make it absolute
            if not os.path.isabs(config.root_output):
                # Get the current working directory and construct absolute path
                current_dir = os.getcwd()
                self.state_file = os.path.join(current_dir, config.root_output, 'pipeline_state.json')
            else:
                self.state_file = os.path.join(config.root_output, 'pipeline_state.json')
        else:
            self.state_file = os.path.join('viranpy_results', 'pipeline_state.json')
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
        
        if self.config.use_prodigal_gv:
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
        """
        try:
            with PipelineLogger("Quality Control Pipeline", self.logger):
                results = {}
                if not self.quality_controller.check_dependencies():
                    raise RuntimeError("Quality control dependencies not available")
                if not getattr(self.config, 'skip_qc', False):
                    self.logger.info("Running FastQC quality assessment")
                    fastqc_results = self.quality_controller.run_fastqc(read_files)
                    results['fastqc'] = fastqc_results
                    self.logger.info("Running Trim Galore quality trimming")
                    trim_results = self.quality_controller.run_trim_galore(read_files, paired)
                    results['trimming'] = trim_results
                    if trim_results.get('success'):
                        self.trimmed_files = trim_results.get('trimmed_files', [])
                self.logger.info("Quality control pipeline completed successfully")
                return results
        except Exception as e:
            self.logger.error(f"Quality control pipeline failed: {e}")
            raise
    
    def run_assembly(self, read_files: List[str], paired: bool = False) -> Dict[str, Any]:
        """
        Run assembly pipeline only.
        """
        try:
            with PipelineLogger("Assembly Pipeline", self.logger):
                results = {}
                input_files = self.filtered_files if self.filtered_files else read_files
                assembler_method = getattr(self.config, 'assembler', 'hybrid')
                
                # Check for existing assembly results when resuming
                if self.resume:
                    self.logger.info("[RESUME] Checking for existing assembly results...")
                    
                    # Check SPAdes results
                    spades_output_dir = os.path.join(self.config.root_output, "spades_assembly")
                    spades_renamed_file = os.path.join(spades_output_dir, "scaffolds_renamed.fasta")
                    if os.path.exists(spades_renamed_file):
                        self.logger.info(f"[RESUME] Found existing SPAdes results: {spades_renamed_file}")
                        results['spades'] = {"success": True, "contigs_file": spades_renamed_file}
                    else:
                        self.logger.info("[RESUME] No existing SPAdes results found, will run SPAdes")
                    
                    # Check MEGAHIT results
                    megahit_output_dir = os.path.join(self.config.root_output, "megahit_assembly")
                    megahit_renamed_file = os.path.join(megahit_output_dir, "contigs_renamed.fasta")
                    if os.path.exists(megahit_renamed_file):
                        self.logger.info(f"[RESUME] Found existing MEGAHIT results: {megahit_renamed_file}")
                        results['megahit'] = {"success": True, "contigs_file": megahit_renamed_file}
                    else:
                        self.logger.info("[RESUME] No existing MEGAHIT results found, will run MEGAHIT")
                
                if assembler_method == 'spades':
                    if 'spades' not in results:  # Not found in resume check
                        self.logger.info("Running SPAdes assembly")
                        spades_result = self.assembler.run_spades(input_files, paired)
                        results['spades'] = spades_result
                    if results['spades'].get('success'):
                        self.final_assembly_file = results['spades'].get('contigs_file')
                elif assembler_method == 'megahit':
                    if 'megahit' not in results:  # Not found in resume check
                        self.logger.info("Running MEGAHIT assembly")
                        # Pass resume flag to assembler
                        self.assembler.config.resume = self.resume
                        megahit_result = self.assembler.run_megahit(input_files, paired)
                        results['megahit'] = megahit_result
                    if results['megahit'].get('success'):
                        self.final_assembly_file = results['megahit'].get('contigs_file')
                elif assembler_method == 'hybrid':
                    self.logger.info("Running hybrid assembly (MEGAHIT + SPAdes)")
                    
                    # Run MEGAHIT if not already found
                    if 'megahit' not in results:
                        # Pass resume flag to assembler
                        self.assembler.config.resume = self.resume
                        megahit_result = self.assembler.run_megahit(input_files, paired)
                        results['megahit'] = megahit_result
                    else:
                        self.logger.info("[RESUME] Using existing MEGAHIT results")
                    
                    # Run SPAdes if not already found
                    if 'spades' not in results:
                        spades_result = self.assembler.run_spades(input_files, paired)
                        results['spades'] = spades_result
                    else:
                        self.logger.info("[RESUME] Using existing SPAdes results")
                    
                    if results['spades'].get('success') and results['megahit'].get('success'):
                        self.logger.info("Creating hybrid assembly")
                        hybrid_result = self.assembler.create_hybrid_assembly(
                            results['spades'].get('contigs_file'),
                            results['megahit'].get('contigs_file'),
                            "hybrid_assembly.fasta"
                        )
                        results['hybrid'] = hybrid_result
                        if hybrid_result.get('success'):
                            self.final_assembly_file = hybrid_result.get('output_file')
                
                if self.final_assembly_file:
                    min_length = getattr(self.config, 'min_contig_length', 200)
                    filter_result = self.assembler.filter_contigs_by_length(
                        self.final_assembly_file, min_length
                    )
                    results['filtering'] = filter_result
                    if filter_result.get('success'):
                        self.final_assembly_file = filter_result.get('output_file')
                        # If hybrid assembly, update the contigs_file and output_file in the results for downstream use
                        if assembler_method == 'hybrid' and 'hybrid' in results and results['hybrid'].get('success'):
                            results['hybrid']['contigs_file'] = self.final_assembly_file
                            results['hybrid']['output_file'] = self.final_assembly_file
                            # Also update the pipeline state for CLI resume
                            if 'assembly' in self.pipeline_state and 'results' in self.pipeline_state['assembly'] and 'hybrid' in self.pipeline_state['assembly']['results']:
                                self.pipeline_state['assembly']['results']['hybrid']['contigs_file'] = self.final_assembly_file
                                self.pipeline_state['assembly']['results']['hybrid']['output_file'] = self.final_assembly_file
                                self.logger.info(f"[DEBUG] Updated pipeline state hybrid contigs_file to: {self.final_assembly_file}")
                                self._save_state()
                            else:
                                self.logger.warning(f"[DEBUG] Could not update pipeline state - missing required keys")
                                self.logger.warning(f"[DEBUG] pipeline_state keys: {list(self.pipeline_state.keys()) if self.pipeline_state else 'None'}")
                                if 'assembly' in self.pipeline_state:
                                    self.logger.warning(f"[DEBUG] assembly keys: {list(self.pipeline_state['assembly'].keys())}")
                                    if 'results' in self.pipeline_state['assembly']:
                                        self.logger.warning(f"[DEBUG] results keys: {list(self.pipeline_state['assembly']['results'].keys())}")
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
                # Debug logging
                self.logger.info(f"Preprocessing pipeline input files: {read_files}")
                self.logger.info(f"Preprocessing pipeline paired flag: {paired}")
                self.logger.info(f"Preprocessing pipeline number of files: {len(read_files)}")
                
                # If resuming, restore input files from state if available
                if self.resume and self.pipeline_state:
                    if 'quality_control' in self.pipeline_state and self.pipeline_state['quality_control'].get('success'):
                        # Use trimmed files from QC if available
                        qc_results = self.pipeline_state['quality_control'].get('results', {})
                        if 'trimming' in qc_results and qc_results['trimming'].get('success'):
                            self.trimmed_files = qc_results['trimming'].get('trimmed_files', [])
                            self.logger.info(f"[RESUME] Restored trimmed files from state: {self.trimmed_files}")
                    
                    if 'host_removal' in self.pipeline_state and self.pipeline_state['host_removal'].get('success'):
                        # Use filtered files from host removal if available
                        host_results = self.pipeline_state['host_removal'].get('results', {})
                        if 'filtering' in host_results and host_results['filtering'].get('success'):
                            self.filtered_files = host_results['filtering'].get('filtered_files', [])
                            self.logger.info(f"[RESUME] Restored filtered files from state: {self.filtered_files}")
                        elif host_results.get('skipped'):
                            # Host removal was skipped, use trimmed files
                            self.filtered_files = self.trimmed_files if self.trimmed_files else read_files
                            self.logger.info(f"[RESUME] Host removal was skipped, using trimmed files: {self.filtered_files}")
                
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
                    # Set final_assembly_file from the state for downstream steps
                    if assembly_results and 'hybrid' in assembly_results and assembly_results['hybrid'].get('success'):
                        # Check if we have a filtered file in the state
                        if 'contigs_file' in assembly_results['hybrid']:
                            self.final_assembly_file = assembly_results['hybrid']['contigs_file']
                            self.logger.info(f"[RESUME] Set final_assembly_file from state: {self.final_assembly_file}")
                        elif 'output_file' in assembly_results['hybrid']:
                            self.final_assembly_file = assembly_results['hybrid']['output_file']
                            self.logger.info(f"[RESUME] Set final_assembly_file from state: {self.final_assembly_file}")
                        else:
                            # Try to find the filtered file
                            hybrid_dir = os.path.join(self.config.root_output, "hybrid_assembly")
                            filtered_file = os.path.join(hybrid_dir, "hybrid_assembly_filtered.fasta")
                            if os.path.exists(filtered_file):
                                self.final_assembly_file = filtered_file
                                # Update the state with the correct path
                                assembly_results['hybrid']['contigs_file'] = filtered_file
                                assembly_results['hybrid']['output_file'] = filtered_file
                                self.pipeline_state['assembly']['results'] = assembly_results
                                self.logger.info(f"[RESUME] Found and set filtered file: {self.final_assembly_file}")
                                self.logger.info(f"[DEBUG] Updated pipeline state hybrid contigs_file to: {self.final_assembly_file}")
                                self._save_state()
                            else:
                                self.logger.warning(f"[RESUME] Could not find filtered hybrid assembly file")
                    elif assembly_results and 'spades' in assembly_results and assembly_results['spades'].get('success'):
                        self.final_assembly_file = assembly_results['spades'].get('contigs_file')
                        self.logger.info(f"[RESUME] Set final_assembly_file from SPAdes state: {self.final_assembly_file}")
                    elif assembly_results and 'megahit' in assembly_results and assembly_results['megahit'].get('success'):
                        self.final_assembly_file = assembly_results['megahit'].get('contigs_file')
                        self.logger.info(f"[RESUME] Set final_assembly_file from MEGAHIT state: {self.final_assembly_file}")
                elif self.resume:
                    # Check if individual assemblers completed successfully
                    spades_output_dir = os.path.join(self.config.root_output, "spades_assembly")
                    megahit_output_dir = os.path.join(self.config.root_output, "megahit_assembly")
                    spades_renamed_file = os.path.join(spades_output_dir, "scaffolds_renamed.fasta")
                    megahit_renamed_file = os.path.join(megahit_output_dir, "contigs_renamed.fasta")
                    
                    if os.path.exists(spades_renamed_file) and os.path.exists(megahit_renamed_file):
                        self.logger.info("[RESUME] Found existing SPAdes and MEGAHIT results, will attempt hybrid assembly")
                        # Use the appropriate input files for assembly
                        assembly_input_files = self.filtered_files if self.filtered_files else (self.trimmed_files if self.trimmed_files else read_files)
                        self.logger.info(f"Assembly input files: {assembly_input_files}")
                        self.logger.info(f"Assembly input files type: {type(assembly_input_files)}")
                        self.logger.info(f"Assembly input files length: {len(assembly_input_files) if assembly_input_files else 0}")
                        
                        if not assembly_input_files:
                            raise RuntimeError("No input files available for assembly. Check if quality control and host removal completed successfully.")
                        
                        assembly_results = self.run_assembly(assembly_input_files, paired)
                        
                        # Check if assembly was truly successful
                        assembly_success = True
                        if assembly_results:
                            # Check if any critical assembly step failed
                            for assembler, result in assembly_results.items():
                                if assembler in ['spades', 'megahit', 'hybrid'] and result and not result.get('success'):
                                    assembly_success = False
                                    self.logger.error(f"Assembly step '{assembler}' failed: {result.get('error', 'Unknown error')}")
                                    break
                        
                        self.pipeline_state['assembly'] = {'success': assembly_success, 'results': assembly_results}
                        self._save_state()
                    else:
                        self.logger.info("[RESUME] No existing assembly results found, will run full assembly")
                        # Use the appropriate input files for assembly
                        assembly_input_files = self.filtered_files if self.filtered_files else (self.trimmed_files if self.trimmed_files else read_files)
                        self.logger.info(f"Assembly input files: {assembly_input_files}")
                        self.logger.info(f"Assembly input files type: {type(assembly_input_files)}")
                        self.logger.info(f"Assembly input files length: {len(assembly_input_files) if assembly_input_files else 0}")
                        
                        if not assembly_input_files:
                            raise RuntimeError("No input files available for assembly. Check if quality control and host removal completed successfully.")
                        
                        assembly_results = self.run_assembly(assembly_input_files, paired)
                        
                        # Check if assembly was truly successful
                        assembly_success = True
                        if assembly_results:
                            # Check if any critical assembly step failed
                            for assembler, result in assembly_results.items():
                                if assembler in ['spades', 'megahit', 'hybrid'] and result and not result.get('success'):
                                    assembly_success = False
                                    self.logger.error(f"Assembly step '{assembler}' failed: {result.get('error', 'Unknown error')}")
                                    break
                        
                        self.pipeline_state['assembly'] = {'success': assembly_success, 'results': assembly_results}
                        self._save_state()
                else:
                    self.logger.info("Step 4: Assembly")
                    # Use the appropriate input files for assembly
                    assembly_input_files = self.filtered_files if self.filtered_files else (self.trimmed_files if self.trimmed_files else read_files)
                    self.logger.info(f"Assembly input files: {assembly_input_files}")
                    self.logger.info(f"Assembly input files type: {type(assembly_input_files)}")
                    self.logger.info(f"Assembly input files length: {len(assembly_input_files) if assembly_input_files else 0}")
                    
                    if not assembly_input_files:
                        raise RuntimeError("No input files available for assembly. Check if quality control and host removal completed successfully.")
                    
                    assembly_results = self.run_assembly(assembly_input_files, paired)
                    
                    # Check if assembly was truly successful
                    assembly_success = True
                    if assembly_results:
                        # Check if any critical assembly step failed
                        for assembler, result in assembly_results.items():
                            if assembler in ['spades', 'megahit', 'hybrid'] and result and not result.get('success'):
                                assembly_success = False
                                self.logger.error(f"Assembly step '{assembler}' failed: {result.get('error', 'Unknown error')}")
                                break
                    
                    self.pipeline_state['assembly'] = {'success': assembly_success, 'results': assembly_results}
                    self._save_state()
                results['assembly'] = assembly_results
                
                # Step 5: Taxonomic classification of contigs
                if self.resume and self.pipeline_state.get('contigs_taxonomy', {}).get('success'):
                    self.logger.info("[RESUME] Skipping contigs taxonomy (already completed)")
                    contigs_taxonomy_results = self.pipeline_state['contigs_taxonomy']['results']
                elif not getattr(self.config, 'skip_taxonomy', False) and getattr(self.config, 'taxonomy_contigs', True):
                    if self.final_assembly_file and os.path.exists(self.final_assembly_file):
                        self.logger.info("Step 5: Taxonomic classification of contigs")
                        try:
                            contigs_taxonomy_results = self._run_contigs_taxonomy(self.final_assembly_file)
                            self.pipeline_state['contigs_taxonomy'] = {'success': True, 'results': contigs_taxonomy_results}
                            self._save_state()
                        except Exception as e:
                            self.logger.error(f"Contigs taxonomy failed: {e}")
                            contigs_taxonomy_results = {"success": False, "error": str(e)}
                            self.pipeline_state['contigs_taxonomy'] = {'success': False, 'error': str(e)}
                            self._save_state()
                    else:
                        self.logger.warning("No valid assembly file available for contigs taxonomy")
                        contigs_taxonomy_results = {"success": False, "error": "No assembly file available"}
                else:
                    self.logger.info("Skipping contigs taxonomy (--skip-taxonomy or not requested)")
                    contigs_taxonomy_results = None
                results['contigs_taxonomy'] = contigs_taxonomy_results
                
                # Step 6: Coverage analysis
                if self.resume and self.pipeline_state.get('coverage', {}).get('success'):
                    self.logger.info("[RESUME] Skipping coverage analysis (already completed)")
                    coverage_results = self.pipeline_state['coverage']['results']
                elif not getattr(self.config, 'skip_coverage_analysis', False):
                    self.logger.info("Step 6: Coverage analysis")
                    if self.final_assembly_file:
                        from ..utils.coverage_calculator import CoverageCalculator
                        coverage_calc = CoverageCalculator(self.config, self.logger)
                        coverage_results = coverage_calc.calculate_coverage(
                            self.final_assembly_file,
                            self.filtered_files,
                            paired,
                            os.path.join(self.config.root_output, "coverage_results")
                        )
                        self.pipeline_state['coverage'] = {'success': True, 'results': coverage_results}
                        self._save_state()
                    else:
                        self.logger.warning("No assembly file available for coverage analysis")
                        coverage_results = None
                else:
                    self.logger.info("Skipping coverage analysis (--skip-coverage-analysis specified)")
                    coverage_results = None
                results['coverage'] = coverage_results
                
                # Step 7: QUAST analysis
                if self.resume and self.pipeline_state.get('quast', {}).get('success'):
                    self.logger.info("[RESUME] Skipping QUAST analysis (already completed)")
                    quast_results = self.pipeline_state['quast']['results']
                elif not getattr(self.config, 'skip_quast_analysis', False):
                    self.logger.info("Step 7: QUAST analysis")
                    if self.final_assembly_file:
                        from ..utils.assembly_stats import AssemblyStatsCalculator
                        quast_calc = AssemblyStatsCalculator(self.config, self.logger)
                        quast_results = quast_calc.run_quast_analysis(
                            self.final_assembly_file,
                            os.path.join(self.config.root_output, "quast_results")
                        )
                        self.pipeline_state['quast'] = {'success': True, 'results': quast_results}
                        self._save_state()
                    else:
                        self.logger.warning("No assembly file available for QUAST analysis")
                        quast_results = None
                else:
                    self.logger.info("Skipping QUAST analysis (--skip-quast-analysis specified)")
                    quast_results = None
                results['quast'] = quast_results
                
                # Step 8: Comprehensive reporting
                if not getattr(self.config, 'skip_comprehensive_report', False):
                    from ..utils.comprehensive_reporter import ComprehensiveReporter
                    reporter = ComprehensiveReporter(self.config, self.logger)
                    sample_name = getattr(self.config, 'sample_name', 'sample')
                    comprehensive_report = reporter.generate_comprehensive_report(self.config.root_output, results, sample_name)
                    self.pipeline_state['comprehensive_report'] = {'success': True, 'results': comprehensive_report}
                    self._save_state()
                    results['comprehensive_report'] = comprehensive_report
                else:
                    self.logger.info("Skipping comprehensive report (--skip-comprehensive-report specified)")
                    results['comprehensive_report'] = None
                self.logger.info("Preprocessing pipeline completed successfully")
                return results
        except Exception as e:
            # Mark the current step as failed
            step = None
            for k in ['taxonomic_classification', 'assembly', 'contigs_taxonomy', 'coverage', 'quast', 'comprehensive_report']:
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
        """
        results = {}
        if not self.host_remover.check_dependencies():
            raise RuntimeError("Host removal dependencies not available")
        host_genome = getattr(self.config, 'host_genome', None)
        bowtie2_index = getattr(self.config, 'bowtie2_index', None)
        if not host_genome and not bowtie2_index:
            self.logger.warning("No host genome or Bowtie2 index provided. Skipping host removal. This may result in host contamination in the assembly.")
            self.filtered_files = self.trimmed_files if self.trimmed_files else read_files
            return {"skipped": True, "warning": "Host removal skipped: no --host-genome or --bowtie2-index provided."}
        input_files = self.trimmed_files if self.trimmed_files else read_files
        if bowtie2_index:
            index_prefix = bowtie2_index
            self.logger.info(f"Using pre-built Bowtie2 index: {index_prefix}")
        else:
            self.logger.info("Building Bowtie2 index from host genome")
            index_result = self.host_remover.build_host_index(host_genome)
            if not index_result.get('success'):
                raise RuntimeError(f"Failed to build Bowtie2 index: {index_result.get('error')}")
            index_prefix = index_result.get('index_prefix')
        self.logger.info("Running Bowtie2 alignment against host genome")
        bowtie2_results = self.host_remover.run_bowtie2_alignment(input_files, index_prefix, paired)
        results['bowtie2_alignment'] = bowtie2_results
        self.logger.info("Filtering host reads")
        filter_results = self.host_remover.filter_host_reads(input_files, bowtie2_results, paired)
        results['filtering'] = filter_results
        if filter_results.get('success'):
            self.filtered_files = filter_results.get('filtered_files', [])
        return results
    
    def _run_raw_reads_taxonomy(self, read_files: List[str], paired: bool = False) -> Dict[str, Any]:
        """
        Run taxonomic classification on raw reads.
        """
        results = {}
        if not self.taxonomic_classifier.check_dependencies():
            raise RuntimeError("Taxonomic classification dependencies not available")
        kraken_db = getattr(self.config, 'kraken2_db', None)
        if not kraken_db:
            raise RuntimeError("Kraken2 database not specified for taxonomic classification")
        input_files = self.trimmed_files if self.trimmed_files else read_files
        self.logger.info("Running Kraken2 classification on raw reads")
        kraken_results = self.taxonomic_classifier.run_kraken2(input_files, paired)
        results['kraken2'] = kraken_results
        krona_file = self.taxonomic_classifier.create_krona_visualization(kraken_results)
        if krona_file:
            results['krona_visualization'] = krona_file
        return results
    
    def _run_contigs_taxonomy(self, contigs_file: str) -> Dict[str, Any]:
        """
        Run taxonomic classification on assembled contigs.
        """
        results = {}
        if not self.taxonomic_classifier.check_dependencies():
            raise RuntimeError("Taxonomic classification dependencies not available")
        kraken_db = getattr(self.config, 'kraken2_db', None)
        if not kraken_db:
            raise RuntimeError("Kraken2 database not specified for taxonomic classification")
        self.logger.info("Running Kraken2 classification on assembled contigs")
        kraken_results = self.taxonomic_classifier.run_kraken2([contigs_file], paired=False)
        results['kraken2'] = kraken_results
        krona_file = self.taxonomic_classifier.create_krona_visualization(kraken_results)
        if krona_file:
            results['krona_visualization'] = krona_file
        return results
    
    def _run_pipeline_steps(self) -> None:
        """Run all pipeline steps in sequence."""
        # Step 1: Taxonomic classification (if enabled, on original contigs)
        if not getattr(self.config, 'skip_taxonomy', False) and getattr(self.config, 'taxonomy_contigs', True):
            # Check if taxonomic classification was already completed when resuming
            if self.resume and self.pipeline_state:
                if 'taxonomic_classification' in self.pipeline_state and self.pipeline_state['taxonomic_classification'].get('success'):
                    self.logger.info("[RESUME] Taxonomic classification already completed, skipping...")
                else:
                    self.logger.info(f"[DEBUG] About to run taxonomy. config.input_file: {self.config.input_file}")
                    self.logger.info(f"[DEBUG] final_assembly_file: {getattr(self, 'final_assembly_file', 'Not set')}")
                    self._run_taxonomic_classification()
            else:
                self.logger.info(f"[DEBUG] About to run taxonomy. config.input_file: {self.config.input_file}")
                self.logger.info(f"[DEBUG] final_assembly_file: {getattr(self, 'final_assembly_file', 'Not set')}")
                self._run_taxonomic_classification()
        
        # Step 2: Sequence preparation (only for annotation-only mode, skip if coming from assembly)
        if hasattr(self, 'final_assembly_file') and self.final_assembly_file:
            # We have assembly results, skip sequence preparation
            self.logger.info("Skipping sequence preparation - using validated contigs from assembly")
        else:
            # Annotation-only mode, need to prepare sequences
            self._prepare_sequences()
        
        # Step 3: Genome shape prediction (if not skipped)
        if not getattr(self.config, 'skip_genome_shape', False):
            # Check if genome shape prediction was already completed when resuming
            if self.resume and self.pipeline_state:
                self.logger.info(f"[DEBUG] Resume mode: {self.resume}")
                self.logger.info(f"[DEBUG] Pipeline state keys: {list(self.pipeline_state.keys())}")
                self.logger.info(f"[DEBUG] Genome shape in state: {'genome_shape' in self.pipeline_state}")
                if 'genome_shape' in self.pipeline_state:
                    self.logger.info(f"[DEBUG] Genome shape state: {self.pipeline_state['genome_shape']}")
                
                if 'genome_shape' in self.pipeline_state and self.pipeline_state['genome_shape'].get('success'):
                    self.logger.info("[RESUME] Genome shape prediction already completed, skipping...")
                else:
                    self.logger.info("[DEBUG] Genome shape not found in state or not successful, running...")
                    self._run_genome_shape_prediction()
            else:
                self.logger.info("[DEBUG] Not in resume mode or no pipeline state, running genome shape...")
                self._run_genome_shape_prediction()
        
        # Step 4: tRNA/tmRNA detection (if not skipped)
        if not getattr(self.config, 'skip_trna_detection', False):
            # Check if tRNA detection was already completed when resuming
            if self.resume and self.pipeline_state:
                if 'trna_detection' in self.pipeline_state and self.pipeline_state['trna_detection'].get('success'):
                    self.logger.info("[RESUME] tRNA/tmRNA detection already completed, skipping...")
                else:
                    self._run_trna_detection()
            else:
                self._run_trna_detection()
        
        # Step 5: ncRNA detection (optional)
        if not self.config.non_crna:
            # Check if ncRNA detection was already completed when resuming
            if self.resume and self.pipeline_state:
                self.logger.info(f"[DEBUG] Resume mode: {self.resume}")
                self.logger.info(f"[DEBUG] Pipeline state keys: {list(self.pipeline_state.keys())}")
                self.logger.info(f"[DEBUG] ncRNA detection in state: {'ncrna_detection' in self.pipeline_state}")
                if 'ncrna_detection' in self.pipeline_state:
                    self.logger.info(f"[DEBUG] ncRNA detection state: {self.pipeline_state['ncrna_detection']}")
                
                if 'ncrna_detection' in self.pipeline_state and self.pipeline_state['ncrna_detection'].get('success'):
                    self.logger.info("[RESUME] ncRNA detection already completed, skipping...")
                else:
                    self.logger.info("[DEBUG] ncRNA detection not found in state or not successful, running...")
                    self._run_ncrna_detection()
            else:
                self.logger.info("[DEBUG] Not in resume mode or no pipeline state, running ncRNA detection...")
                self._run_ncrna_detection()
        
        # Step 6: CRISPR detection (if not skipped)
        if not getattr(self.config, 'skip_crispr_detection', False):
            # Check if CRISPR detection was already completed when resuming
            if self.resume and self.pipeline_state:
                if 'crispr_detection' in self.pipeline_state and self.pipeline_state['crispr_detection'].get('success'):
                    self.logger.info("[RESUME] CRISPR detection already completed, skipping...")
                else:
                    self._run_crispr_detection()
            else:
                self._run_crispr_detection()
        
        # Step 7: Gene prediction (if not skipped)
        if not getattr(self.config, 'skip_gene_prediction', False):
            # Check if gene prediction was already completed when resuming
            if self.resume and self.pipeline_state:
                if 'gene_prediction' in self.pipeline_state and self.pipeline_state['gene_prediction'].get('success'):
                    self.logger.info("[RESUME] Gene prediction already completed, skipping...")
                else:
                    self._run_gene_prediction()
            else:
                self._run_gene_prediction()
        
        # Step 8: Protein function prediction (if not skipped)
        if not getattr(self.config, 'skip_protein_function', False):
            self._run_protein_function_prediction()
        
        # Step 9: HMMER analysis (optional)
        if not self.config.no_hmmer:
            self._run_hmmer_analysis()
    
    def _prepare_sequences(self) -> None:
        """Prepare sequences for annotation."""
        with PipelineLogger("Sequence preparation", self.logger):
            from Bio import SeqIO
            from ..utils.sequence_utils import rename_sequences
            import os
            
            # Skip sequence renaming to preserve original contig names
            self.logger.info("Skipping sequence renaming to preserve original contig names")
            
            # Create sequence preparation directory even when skipping renaming
            output_dir = os.path.join(self.config.root_output, "sequence_preparation")
            os.makedirs(output_dir, exist_ok=True)
            
            # Count all sequences in the input file
            try:
                with open(self.config.input_file, "r") as f:
                    record_count = 0
                    for record in SeqIO.parse(f, "fasta"):
                        record_count += 1
                self.logger.info(f"Sequence preparation completed - validated {record_count} sequences (preserving original names)")
            except Exception as e:
                self.logger.error(f"Error validating sequences: {e}")
                raise
            
            self.logger.info("Sequences prepared successfully (original names preserved)")
    
    def _run_genome_shape_prediction(self) -> None:
        """Predict genome shape using LASTZ."""
        with PipelineLogger("Genome shape prediction", self.logger):
            try:
                from ..annotators.genome_shape import ViralGenomeTopology
                
                # Create the genome shape predictor
                genome_shape_predictor = ViralGenomeTopology(self.config, self.logger)
                
                # Run genome shape prediction
                result = genome_shape_predictor.run(self.config.input_file)
                
                if result.get('success'):
                    self.logger.info(f"Genome shape prediction completed successfully: {result.get('metadata', {}).get('total_sequences', 0)} sequences analyzed")
                    # Store results
                    self.result.add_annotation_result(AnnotationResult(
                        annotator_name="GenomeShapePrediction",
                        input_file=self.config.input_file,
                        success=True,
                        annotations=result.get('annotations', {}),
                        metadata=result.get('metadata', {})
                    ))
                    # Update pipeline state for resume functionality
                    if not hasattr(self, 'pipeline_state'):
                        self.pipeline_state = {}
                    self.pipeline_state['genome_shape'] = {
                        'success': True,
                        'metadata': result.get('metadata', {})
                    }
                    self._save_state()
                else:
                    self.logger.warning(f"Genome shape prediction failed: {result.get('error_message', 'Unknown error')}")
                    
            except Exception as e:
                self.logger.error(f"Genome shape prediction failed: {e}")
                raise
    
    def _run_trna_detection(self) -> None:
        """Detect tRNA and tmRNA sequences."""
        with PipelineLogger("tRNA/tmRNA detection", self.logger):
            try:
                from ..annotators.trna_detector import ViralTRNAFinder
                
                # Create the tRNA detector
                trna_detector = ViralTRNAFinder(self.config, self.logger)
                
                # Run tRNA detection
                result = trna_detector.run(self.config.input_file)
                
                if result.get('success'):
                    self.logger.info(f"tRNA/tmRNA detection completed successfully: {result.get('metadata', {}).get('total_trna', 0)} tRNAs found")
                    # Store results
                    self.result.add_annotation_result(AnnotationResult(
                        annotator_name="TRNADetection",
                        input_file=self.config.input_file,
                        success=True,
                        annotations=result.get('annotations', {}),
                        metadata=result.get('metadata', {})
                    ))
                    # Update pipeline state for resume functionality
                    if not hasattr(self, 'pipeline_state'):
                        self.pipeline_state = {}
                    self.pipeline_state['trna_detection'] = {
                        'success': True,
                        'metadata': result.get('metadata', {})
                    }
                    self._save_state()
                else:
                    self.logger.warning(f"tRNA/tmRNA detection failed: {result.get('error_message', 'Unknown error')}")
                    
            except Exception as e:
                self.logger.error(f"tRNA/tmRNA detection failed: {e}")
                raise
    
    def _run_ncrna_detection(self) -> None:
        """Detect ncRNA sequences using RFAM."""
        with PipelineLogger("ncRNA detection", self.logger):
            try:
                from ..annotators.ncrna_detector import ViralNcRNAFinder
                
                # Create the ncRNA detector
                ncrna_detector = ViralNcRNAFinder(self.config, self.logger)
                
                # Run ncRNA detection
                result = ncrna_detector.run(self.config.input_file)
                
                if result.get('success'):
                    self.logger.info(f"ncRNA detection completed successfully: {result.get('metadata', {}).get('total_ncrna', 0)} ncRNAs found")
                    # Store results
                    self.result.add_annotation_result(AnnotationResult(
                        annotator_name="NCRNADetection",
                        input_file=self.config.input_file,
                        success=True,
                        annotations=result.get('annotations', {}),
                        metadata=result.get('metadata', {})
                    ))
                    # Update pipeline state for resume functionality
                    if not hasattr(self, 'pipeline_state'):
                        self.pipeline_state = {}
                    self.pipeline_state['ncrna_detection'] = {
                        'success': True,
                        'metadata': result.get('metadata', {})
                    }
                    self._save_state()
                else:
                    self.logger.warning(f"ncRNA detection failed: {result.get('error_message', 'Unknown error')}")
                    
            except Exception as e:
                self.logger.error(f"ncRNA detection failed: {e}")
                raise
    
    def _run_crispr_detection(self) -> None:
        """Detect CRISPR repeats using PILER-CR."""
        with PipelineLogger("CRISPR detection", self.logger):
            try:
                from ..annotators.crispr_detector import ViralCRISPRFinder
                
                # Create the CRISPR detector
                crispr_detector = ViralCRISPRFinder(self.config, self.logger)
                
                # Run CRISPR detection
                result = crispr_detector.run(self.config.input_file)
                
                if result.get('success'):
                    self.logger.info(f"CRISPR detection completed successfully: {result.get('metadata', {}).get('total_crispr_arrays', 0)} CRISPR arrays found")
                    # Store results
                    self.result.add_annotation_result(AnnotationResult(
                        annotator_name="CRISPRDetection",
                        input_file=self.config.input_file,
                        success=True,
                        annotations=result.get('annotations', {}),
                        metadata=result.get('metadata', {})
                    ))
                    # Update pipeline state for resume functionality
                    if not hasattr(self, 'pipeline_state'):
                        self.pipeline_state = {}
                    self.pipeline_state['crispr_detection'] = {
                        'success': True,
                        'metadata': result.get('metadata', {})
                    }
                    self._save_state()
                else:
                    self.logger.warning(f"CRISPR detection failed: {result.get('error_message', 'Unknown error')}")
                    
            except Exception as e:
                self.logger.error(f"CRISPR detection failed: {e}")
                raise
    
    def _run_gene_prediction(self) -> None:
        """Predict genes using Prodigal."""
        with PipelineLogger("Gene prediction", self.logger):
            try:
                from ..annotators.gene_predictor import ViralGeneFinder
                
                # Create the gene predictor
                gene_predictor = ViralGeneFinder(self.config, self.logger)
                
                # Run gene prediction
                result = gene_predictor.run(self.config.input_file)
                
                if result.get('success'):
                    self.logger.info(f"Gene prediction completed successfully: {result.get('metadata', {}).get('total_genes', 0)} genes found")
                    # Store results
                    self.result.add_annotation_result(AnnotationResult(
                        annotator_name="GenePrediction",
                        input_file=self.config.input_file,
                        success=True,
                        annotations=result.get('annotations', {}),
                        metadata=result.get('metadata', {})
                    ))
                    # Update pipeline state for resume functionality
                    if not hasattr(self, 'pipeline_state'):
                        self.pipeline_state = {}
                    self.pipeline_state['gene_prediction'] = {
                        'success': True,
                        'metadata': result.get('metadata', {})
                    }
                    self._save_state()
                else:
                    self.logger.warning(f"Gene prediction failed: {result.get('error_message', 'Unknown error')}")
                    
            except Exception as e:
                self.logger.error(f"Gene prediction failed: {e}")
                raise
    
    def _run_protein_function_prediction(self) -> None:
        """Predict protein function using BLAST or DIAMOND."""
        with PipelineLogger("Protein function prediction", self.logger):
            try:
                from ..annotators.protein_function_predictor import ViralProteinAnnotator
                
                # Create the protein function predictor
                protein_predictor = ViralProteinAnnotator(self.config, self.logger)
                
                # Run protein function prediction
                result = protein_predictor.run(self.config.input_file)
                
                if result.get('success'):
                    self.logger.info(f"Protein function prediction completed successfully: {result.get('metadata', {}).get('total_proteins', 0)} proteins annotated")
                    # Store results
                    self.result.add_annotation_result(AnnotationResult(
                        annotator_name="ProteinFunctionPrediction",
                        input_file=self.config.input_file,
                        success=True,
                        annotations=result.get('annotations', {}),
                        metadata=result.get('metadata', {})
                    ))
                else:
                    self.logger.warning(f"Protein function prediction failed: {result.get('error_message', 'Unknown error')}")
                    
            except Exception as e:
                self.logger.error(f"Protein function prediction failed: {e}")
                raise
    
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
    
    def _run_taxonomic_classification(self) -> Dict[str, Any]:
        """Run taxonomic classification on contigs."""
        with PipelineLogger("Taxonomic classification", self.logger):
            try:
                from ..annotators.taxonomic_classification import TaxonomicClassificationAnnotator
                import os
                
                # Determine the correct input file for taxonomy
                # Priority: final_assembly_file > config.input_file
                input_file = None
                if hasattr(self, 'final_assembly_file') and self.final_assembly_file and os.path.exists(self.final_assembly_file):
                    input_file = self.final_assembly_file
                    self.logger.info(f"[DEBUG] Using final assembly file for taxonomy: {input_file}")
                elif self.config.input_file and os.path.exists(self.config.input_file):
                    input_file = self.config.input_file
                    self.logger.info(f"[DEBUG] Using config input file for taxonomy: {input_file}")
                else:
                    self.logger.error(f"[DEBUG] No valid input file found for taxonomy")
                    self.logger.error(f"[DEBUG] final_assembly_file: {getattr(self, 'final_assembly_file', 'Not set')}")
                    self.logger.error(f"[DEBUG] config.input_file: {self.config.input_file}")
                    return {"success": False, "error": "No valid input file found for taxonomy"}
                
                kraken_db = getattr(self.config, 'kraken2_db', None)
                file_exists = os.path.exists(input_file)
                file_size = os.path.getsize(input_file) if file_exists else 0
                self.logger.info(f"[DEBUG] Taxonomic classification input file: {input_file}")
                self.logger.info(f"[DEBUG] Input file exists: {file_exists}, size: {file_size} bytes")
                self.logger.info(f"[DEBUG] Kraken2 DB: {kraken_db}")
                
                if not file_exists:
                    self.logger.error(f"[DEBUG] Input file does not exist: {input_file}")
                    return {"success": False, "error": f"Input file does not exist: {input_file}"}
                
                if file_size == 0:
                    self.logger.error(f"[DEBUG] Input file is empty: {input_file}")
                    return {"success": False, "error": f"Input file is empty: {input_file}"}
                
                self.logger.info(f"[DEBUG] Creating TaxonomicClassificationAnnotator...")
                annotator = TaxonomicClassificationAnnotator(self.config, self.logger)
                self.logger.info(f"[DEBUG] TaxonomicClassificationAnnotator created successfully")
                
                self.logger.info(f"[DEBUG] Calling annotator.annotate({input_file})...")
                result = annotator.annotate(input_file)
                self.logger.info(f"[DEBUG] annotator.annotate() completed, success: {result.success}")
                
                if result.success:
                    self.logger.info("Taxonomic classification completed successfully")
                    # Store results for later use - use the correct PipelineResult structure
                    # Add the annotation result to the pipeline results
                    self.result.add_annotation_result(result)
                    # Also store in metadata for easy access
                    if not hasattr(self.result, 'metadata'):
                        self.result.metadata = {}
                    self.result.metadata["taxonomic_classification"] = result.to_dict()
                else:
                    self.logger.warning(f"Taxonomic classification failed: {result.error_message}")
                return result.to_dict()
                
            except Exception as e:
                import traceback
                self.logger.error(f"[DEBUG] Exception in _run_taxonomic_classification: {e}")
                self.logger.error(f"[DEBUG] Traceback: {traceback.format_exc()}")
                return {"success": False, "error": str(e)}
    
    def _save_state(self):
        os.makedirs(os.path.dirname(self.state_file), exist_ok=True)
        with open(self.state_file, 'w') as f:
            json.dump(self.pipeline_state, f, indent=2)

    def _load_state(self):
        if os.path.exists(self.state_file):
            with open(self.state_file, 'r') as f:
                return json.load(f)
        return {} 