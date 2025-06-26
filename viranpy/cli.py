#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Command-line interface for ViRAnPy pipeline.
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Optional
import os
import shutil

from .config import PipelineConfig
from .core.pipeline import ViralAnnotationPipeline
from .utils.logger import setup_logger, logger


def create_parser() -> argparse.ArgumentParser:
    """
    Create the command-line argument parser.
    
    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description='ViRAnPy - Viral Metagenomic Analysis Pipeline\n\n'
                   'Primary Workflow: ViRAnPy is designed for viral metagenomic analysis starting from raw reads. '
                   'The default workflow includes quality control, host removal, assembly, and viral annotation.\n\n'
                   'Workflow Options:\n'
                   '  • Read files + no flags: Full pipeline (preprocessing + assembly + annotation)\n'
                   '  • Read files + --assemble-only: Preprocessing + assembly only\n'
                   '  • Read files + --qc-only: Quality control only\n'
                   '  • --input FASTA: Direct annotation of pre-assembled contigs\n\n'
                   'Database Management: ViRAnPy automatically manages all required databases. '
                   'Use --build-databases to install databases and --check-databases to verify installation.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Install all required databases (run first time)
  viranpy --build-databases
  
  # Check database status
  viranpy --check-databases
  
  # Generate viral metadata for metagenomic contigs
  viranpy --generate-metadata --sample-name VIROME001 --viral-family Myoviridae --viral-genus T4virus --source gut --location USA --collection-date 2024-01-15
  
  # Generate metadata for segmented viruses
  viranpy --generate-metadata --sample-name SEGMENTED001 --viral-family Orthomyxoviridae --viral-genus Influenzavirus --segments segment1 segment2 segment3 --source human --location Canada
  
  # Full viral metagenomic pipeline with paired-end reads (explicit)
  viranpy --pe1 R1.fastq --pe2 R2.fastq --generate-metadata --sample-name VIROME001 --host-genome host.fasta --taxonomy-contigs --coverage-analysis --quast-analysis
  
  # Full viral metagenomic pipeline with paired-end reads (automatic detection)
  viranpy --reads R1.fastq R2.fastq --generate-metadata --sample-name VIROME001 --host-genome host.fasta --taxonomy-contigs --coverage-analysis --quast-analysis
  
  # Full pipeline with single-end reads (explicit)
  viranpy --single reads.fastq --viral-metadata viral_metadata.txt --host-genome host.fasta --taxonomy-contigs --coverage-analysis --quast-analysis
  
  # Full pipeline with single-end reads (automatic detection)
  viranpy --reads reads.fastq --viral-metadata viral_metadata.txt --host-genome host.fasta --taxonomy-contigs --coverage-analysis --quast-analysis
  
  # Full pipeline with comprehensive reporting (annotation runs by default)
  viranpy --pe1 R1.fastq --pe2 R2.fastq --viral-metadata viral_metadata.txt --host-genome host.fasta --taxonomy-raw-reads --taxonomy-contigs --coverage-analysis --quast-analysis --comprehensive-report
  
  # Assembly only (skip annotation)
  viranpy --pe1 R1.fastq --pe2 R2.fastq --viral-metadata viral_metadata.txt --assemble-only --host-genome host.fasta --coverage-analysis --quast-analysis
  
  # Assembly with selective annotation (skip specific steps)
  viranpy --pe1 R1.fastq --pe2 R2.fastq --viral-metadata viral_metadata.txt --host-genome host.fasta --skip-crispr-detection --skip-trna-detection
  
  # Quality control and host removal only
  viranpy --pe1 R1.fastq --pe2 R2.fastq --viral-metadata viral_metadata.txt --qc-only --host-genome host.fasta
  
  # Direct annotation of pre-assembled contigs (ADVANCED)
  viranpy --input genome.fasta --viral-metadata viral_metadata.txt
  
  # Using pre-built Bowtie2 index
  viranpy --pe1 R1.fastq --pe2 R2.fastq --viral-metadata viral_metadata.txt --bowtie2-index /path/to/host_index --coverage-analysis

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
🎉 Thank you for using ViRAnPy! 🦠

📖 Please cite:
  Duhan N. et al., ViRAnPy
  https://github.com/naveenduhan/viranpy

💡 "In bioinformatics, the answers you get depend on the questions you ask."
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        """
    )
    
    # Workflow options
    workflow_group = parser.add_argument_group('Workflow options [REQUIRED]')
    workflow_group.add_argument(
        "--reads", dest="read_files", type=str, nargs='+',
        help='Input FASTQ files for preprocessing and assembly. For paired-end: provide 2 files in order (R1.fastq R2.fastq). For single-end: provide 1 file. Note: Use --pe1/--pe2 for explicit paired-end control.', metavar="FASTQFILES"
    )
    workflow_group.add_argument(
        "--pe1", dest="pe1_file", type=str,
        help='Forward/Read 1 FASTQ file for paired-end reads', metavar="FORWARD_FASTQ"
    )
    workflow_group.add_argument(
        "--pe2", dest="pe2_file", type=str,
        help='Reverse/Read 2 FASTQ file for paired-end reads', metavar="REVERSE_FASTQ"
    )
    workflow_group.add_argument(
        "--single", dest="single_file", type=str,
        help='Single-end FASTQ file', metavar="SINGLE_FASTQ"
    )
    workflow_group.add_argument(
        "--qc-only", dest="qc_only", action='store_true', default=False,
        help='Run quality control only (Default: False)'
    )
    workflow_group.add_argument(
        "--assemble-only", dest="assemble_only", action='store_true', default=False,
        help='Run preprocessing and assembly only (skip viral annotation) (Default: False)'
    )
    
    # Direct annotation options (for pre-assembled contigs)
    direct_annotation_group = parser.add_argument_group('Direct annotation options [OPTIONAL]')
    direct_annotation_group.add_argument(
        "--input", dest="input_file", type=str,
        help='Input FASTA file for direct annotation (use for pre-assembled contigs)', metavar="FASTAFILE"
    )
    
    # Basic options [REQUIRED]
    basic_group = parser.add_argument_group('Basic options [REQUIRED]')
    basic_group.add_argument(
        "--viral-metadata", dest="viral_metadata_file", type=str,
        help='Input file as a plain text file with the viral metadata per every FASTA header',
        metavar="TEXTFILE"
    )
    basic_group.add_argument(
        "--generate-metadata", dest="generate_metadata", action='store_true', default=False,
        help='Generate viral metadata file for metagenomic contigs (use with --sample-name)'
    )
    basic_group.add_argument(
        "--build-databases", dest="build_databases", action='store_true', default=False,
        help='Download and install all required databases (RFAM, DIAMOND, BLAST, VOGS, RVDB, PHROGS, VFAM)'
    )
    basic_group.add_argument(
        "--check-databases", dest="check_databases", action='store_true', default=False,
        help='Check if all required databases are installed and ready for use'
    )
    basic_group.add_argument(
        "--fix-missing", dest="fix_missing", action='store_true', default=False,
        help='When used with --build-databases, attempt to repair or generate missing database index files if possible.'
    )
    basic_group.add_argument(
        "--sample-name", dest="sample_name", type=str,
        help='Sample name/ID for viral metagenomic metadata generation', metavar="NAME"
    )
    basic_group.add_argument(
        "--viral-family", dest="viral_family", type=str,
        help='Viral family if known', metavar="FAMILY"
    )
    basic_group.add_argument(
        "--viral-genus", dest="viral_genus", type=str,
        help='Viral genus if known', metavar="GENUS"
    )
    basic_group.add_argument(
        "--viral-species", dest="viral_species", type=str,
        help='Viral species if known', metavar="SPECIES"
    )
    basic_group.add_argument(
        "--segments", dest="segments", type=str, nargs='+',
        help='Segment identifiers for segmented viruses', metavar="SEGMENTS"
    )
    basic_group.add_argument(
        "--host", dest="host", type=str,
        help='Host organism if known', metavar="HOST"
    )
    basic_group.add_argument(
        "--collection-date", dest="collection_date", type=str,
        help='Collection date in YYYY-MM-DD format', metavar="DATE"
    )
    basic_group.add_argument(
        "--location", dest="location", type=str,
        help='Geographic location/country', metavar="LOCATION"
    )
    basic_group.add_argument(
        "--source", dest="source", type=str,
        help='Sample source (soil, water, gut, etc.)', metavar="SOURCE"
    )
    
    # Analysis options
    analysis_group = parser.add_argument_group('Analysis options [OPTIONAL]')
    analysis_group.add_argument(
        "--skip-coverage-analysis", dest="skip_coverage_analysis", action='store_true', default=False,
        help='Skip coverage analysis for contigs using BWA and samtools (Default: False - coverage analysis runs by default)'
    )
    analysis_group.add_argument(
        "--skip-quast-analysis", dest="skip_quast_analysis", action='store_true', default=False,
        help='Skip QUAST for assembly quality assessment (Default: False - QUAST analysis runs by default)'
    )
    analysis_group.add_argument(
        "--skip-comprehensive-report", dest="skip_comprehensive_report", action='store_true', default=False,
        help='Skip comprehensive HTML report generation (Default: False - comprehensive report generated by default)'
    )
    
    # Quality control options
    qc_group = parser.add_argument_group('Quality control options [OPTIONAL]')
    qc_group.add_argument(
        "--skip-qc", dest="skip_qc", action='store_true', default=False,
        help='Skip quality control (Default: False)'
    )
    qc_group.add_argument(
        "--trim-quality", dest="trim_quality", type=int, default=20,
        help='Quality threshold for trimming (Default: 20)', metavar="INT"
    )
    qc_group.add_argument(
        "--min-length", dest="min_length", type=int, default=20,
        help='Minimum read length after trimming (Default: 20)', metavar="INT"
    )
    qc_group.add_argument(
        "--adapter", dest="adapter", type=str,
        help='Adapter sequence for trimming', metavar="SEQUENCE"
    )
    
    # Host removal options
    host_group = parser.add_argument_group('Host removal options [OPTIONAL]')
    host_group.add_argument(
        "--skip-host-removal", dest="skip_host_removal", action='store_true', default=False,
        help='Skip host removal (Default: False)'
    )
    host_group.add_argument(
        "--host-genome", dest="host_genome", type=str,
        help='Host genome FASTA file for Bowtie2 alignment', metavar="FASTAFILE"
    )
    host_group.add_argument(
        "--bowtie2-index", dest="bowtie2_index", type=str,
        help='Pre-built Bowtie2 index prefix (if not provided, will build from host genome)', metavar="PREFIX"
    )
    
    # Taxonomic classification options
    taxonomy_group = parser.add_argument_group('Taxonomic classification options [OPTIONAL]')
    taxonomy_group.add_argument(
        "--skip-taxonomy", dest="skip_taxonomy", action='store_true', default=False,
        help='Skip taxonomic classification (Default: False)'
    )
    taxonomy_group.add_argument(
        "--taxonomy-raw-reads", dest="taxonomy_raw_reads", action='store_true', default=False,
        help='Run taxonomic classification on raw reads (Default: False)'
    )
    taxonomy_group.add_argument(
        "--taxonomy-contigs", dest="taxonomy_contigs", action='store_true', default=True,
        help='Run taxonomic classification on assembled contigs (Default: True)'
    )
    
    # Assembly options
    assembly_group = parser.add_argument_group('Assembly options [OPTIONAL]')
    assembly_group.add_argument(
        "--assembler", dest="assembler", type=str, choices=['spades', 'megahit', 'hybrid'],
        default='hybrid', help='Assembly method (Default: hybrid)', metavar="METHOD"
    )
    assembly_group.add_argument(
        "--spades-memory", dest="spades_memory", type=int, default=16,
        help='SPAdes memory limit in GB (Default: 16). Use --memory for other tools if needed.', metavar="INT"
    )
    assembly_group.add_argument(
        "--cdhit-identity", dest="cdhit_identity", type=float, default=0.95,
        help='CD-HIT identity threshold (Default: 0.95)', metavar="FLOAT"
    )
    assembly_group.add_argument(
        "--min-contig-length", dest="min_contig_length", type=int, default=200,
        help='Minimum contig length to keep (Default: 200)', metavar="INT"
    )
    
    # Annotation options
    annotation_group = parser.add_argument_group('Annotation options [OPTIONAL]')
    annotation_group.add_argument(
        "--skip-annotation", dest="skip_annotation", action='store_true', default=False,
        help='Skip viral annotation pipeline after assembly (Default: False - annotation runs by default)'
    )
    annotation_group.add_argument(
        "--skip-genome-shape", dest="skip_genome_shape", action='store_true', default=False,
        help='Skip genome shape prediction (Default: False)'
    )
    annotation_group.add_argument(
        "--skip-trna-detection", dest="skip_trna_detection", action='store_true', default=False,
        help='Skip tRNA/tmRNA detection (Default: False)'
    )
    annotation_group.add_argument(
        "--skip-crispr-detection", dest="skip_crispr_detection", action='store_true', default=False,
        help='Skip CRISPR detection (Default: False)'
    )
    annotation_group.add_argument(
        "--skip-gene-prediction", dest="skip_gene_prediction", action='store_true', default=False,
        help='Skip gene prediction (Default: False)'
    )
    annotation_group.add_argument(
        "--skip-protein-function", dest="skip_protein_function", action='store_true', default=False,
        help='Skip protein function prediction (Default: False)'
    )
    
    # Advanced general options [OPTIONAL]
    advanced_general_group = parser.add_argument_group('Advanced general options [OPTIONAL]')
    advanced_general_group.add_argument(
        "--out", dest="root_output", type=str,
        help='Name of the outputs files (without extension)', metavar="OUTPUTNAME"
    )
    advanced_general_group.add_argument(
        "--locus", dest="locus", type=str, default='LOC',
        help='Name of the sequences (Default: LOC)', metavar="STRING"
    )
    ncpus_default = os.cpu_count() or 1
    advanced_general_group.add_argument(
        "--threads", dest="ncpus", type=int, default=ncpus_default,
        help='Number of threads/CPUs to use for all steps (Default: all available)', metavar="INT"
    )
    advanced_general_group.add_argument(
        "--mincontigsize", dest="min_contig_size", type=int, default=200,
        help='Minimum contig length to be considered (Default: 200 bp)', metavar="INT"
    )
    advanced_general_group.add_argument(
        "--blast", dest="blast_switch", action='store_true', default=False,
        help='Using BLAST to predict protein function based on homology (Default: False)'
    )
    advanced_general_group.add_argument(
        "--memory", dest="memory", type=str, default=None,
        help='Memory limit for tools (e.g., "16G", "8GB", "16384M"). For MEGAHIT on Mac, this is required and converted to bytes. Default: system-based (16GB if >=16GB RAM, 75%% of system RAM otherwise)', metavar="MEMORY"
    )
    
    # Advanced circularity options
    circularity_group = parser.add_argument_group('Advanced options for origin/terminus prediction [OPTIONAL]')
    circularity_group.add_argument(
        "--gc-skew-read-length", dest="gc_skew_read_length", type=int, default=101,
        help='Read length for GC skew-based origin/terminus prediction (default: 101 bp)', metavar="INT"
    )
    circularity_group.add_argument(
        "--gc-skew-window", dest="gc_skew_window", type=int, default=100,
        help='Window size for GC skew calculation (default: 100 bp)', metavar="INT"
    )
    circularity_group.add_argument(
        "--gc-skew-slide", dest="gc_skew_slide", type=int, default=10,
        help='Sliding window size for GC skew analysis (default: 10 bp)', metavar="INT"
    )
    
    # Advanced CRISPR/repeats options
    crispr_group = parser.add_argument_group('Advanced options for CRISPR/repeat detection [OPTIONAL]')
    crispr_group.add_argument(
        "--min-crispr-repeat", dest="min_crispr_repeat", type=int, default=16,
        help="Minimum CRISPR repeat length (Default: 16)", metavar="INT"
    )
    crispr_group.add_argument(
        "--max-crispr-repeat", dest="max_crispr_repeat", type=int, default=64,
        help="Maximum CRISPR repeat length (Default: 64)", metavar="INT"
    )
    crispr_group.add_argument(
        "--min-crispr-spacer", dest="min_crispr_spacer", type=int, default=8,
        help="Minimum CRISPR spacer length (Default: 8)", metavar="INT"
    )
    crispr_group.add_argument(
        "--max-crispr-spacer", dest="max_crispr_spacer", type=int, default=64,
        help="Maximum CRISPR spacer length (Default: 64)", metavar="INT"
    )
    
    # Advanced genetic code options
    gcode_group = parser.add_argument_group('Advanced options for gene prediction and translation table [OPTIONAL]')
    gcode_group.add_argument(
        "--use-prodigal-gv", dest="use_prodigal_gv", action='store_true', default=False,
        help='Use Prodigal-GV for gene prediction (Default: False)'
    )
    gcode_group.add_argument(
        "--genetic-code-table", dest="genetic_code_table", type=str, default='11',
        help='NCBI genetic code table number (Default: 11)', metavar="NUMBER"
    )
    
    # Advanced GenBank division options
    division_group = parser.add_argument_group('Advanced options for GenBank division [OPTIONAL]')
    division_group.add_argument(
        "--genbank-division", dest="genbank_division", type=str, default='CON',
        help='GenBank Division: BCT|CON|VRL|PHG (Default: CON)', metavar="BCT|CON|VRL|PHG"
    )
    
    # Advanced ncRNA options
    rfam_group = parser.add_argument_group('Advanced options for ncRNA [OPTIONAL]\n'
                                         'Note: Database paths are automatically managed by ViRAnPy')
    rfam_group.add_argument(
        "--noncrna", dest="non_crna", action='store_true', default=False,
        help="Don't use --ncrna option in Prodigal (Default: False)"
    )
    rfam_group.add_argument(
        "--norfam", dest="nor_fam", action='store_true', default=False,
        help="Don't use --rfam option in Infernal (Default: False)"
    )
    rfam_group.add_argument(
        "--hmmonly", dest="hmmonly", action='store_true', default=False,
        help="Use --hmmonly instead of --nohmmonly (Default: False)"
    )
    
    # Advanced DIAMOND options
    diamond_group = parser.add_argument_group('Advanced options for DIAMOND [OPTIONAL]\n'
                                            'Note: Database paths are automatically managed by ViRAnPy')
    diamond_group.add_argument(
        "--diamondevalue", dest="diamond_evalue", type=float, default=0.00001,
        help='DIAMOND e-value threshold (Default: 0.00001)', metavar="FLOAT"
    )
    diamond_group.add_argument(
        "--diamondwidthr", dest="diamond_width_threshold", type=float, default=50.0,
        help='DIAMOND ID threshold (Default: 50.0)', metavar="FLOAT"
    )
    diamond_group.add_argument(
        "--diamondcoverthr", dest="diamond_cov_threshold", type=float, default=50.0,
        help='DIAMOND Coverage threshold (Default: 50.0)', metavar="FLOAT"
    )
    
    # Advanced BLAST options
    blast_group = parser.add_argument_group('Advanced options for BLAST [OPTIONAL]\n'
                                          'Note: Database paths are automatically managed by ViRAnPy')
    blast_group.add_argument(
        "--blastexh", dest="blast_exh", action='store_true', default=False,
        help='Use exhaustive BLAST mode (Default: False)'
    )
    blast_group.add_argument(
        "--blastevalue", dest="blast_evalue", type=float, default=0.00001,
        help='BLAST e-value threshold (Default: 0.00001)', metavar="FLOAT"
    )
    blast_group.add_argument(
        "--blastwidthr", dest="blast_width_threshold", type=float, default=50.0,
        help='BLAST ID threshold (Default: 50.0)', metavar="FLOAT"
    )
    blast_group.add_argument(
        "--blastcoverthr", dest="blast_cov_threshold", type=float, default=50.0,
        help='BLAST Coverage threshold (Default: 50.0)', metavar="FLOAT"
    )
    
    # Advanced HMMER options
    hmmer_group = parser.add_argument_group('Advanced options for HMMER [OPTIONAL]\n'
                                          'Note: Database paths are automatically managed by ViRAnPy')
    hmmer_group.add_argument(
        "--nohmmer", dest="no_hmmer", action='store_true', default=False,
        help='Skip HMMER analysis (Default: False)'
    )
    hmmer_group.add_argument(
        "--novogs", dest="no_vogs", action='store_true', default=False,
        help='Skip VOGS database (Default: False)'
    )
    hmmer_group.add_argument(
        "--norvdb", dest="no_rvdb", action='store_true', default=False,
        help='Skip RVDB database (Default: False)'
    )
    hmmer_group.add_argument(
        "--nophrogs", dest="no_phrogs", action='store_true', default=False,
        help='Skip PHROGS database (Default: False)'
    )
    hmmer_group.add_argument(
        "--novfam", dest="no_vfam", action='store_true', default=False,
        help='Skip VFAM database (Default: False)'
    )
    hmmer_group.add_argument(
        "--hmmerevalue", dest="hmmer_evalue", type=float, default=0.001,
        help='HMMER e-value threshold (Default: 0.001)', metavar="FLOAT"
    )
    hmmer_group.add_argument(
        "--hmmercoverthr", dest="hmmer_cov_threshold", type=float, default=50.0,
        help='HMMER Coverage threshold (Default: 50.0)', metavar="FLOAT"
    )
    
    # Resume option
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume the pipeline from the last successful step using the pipeline_state.json file."
    )
    
    return parser


def parse_args(args: Optional[list] = None) -> argparse.Namespace:
    """
    Parse command-line arguments.
    
    Args:
        args: Command-line arguments (if None, uses sys.argv)
        
    Returns:
        Parsed arguments
    """
    parser = create_parser()
    parsed_args = parser.parse_args(args)
    
    # Handle database management options (these can run independently)
    if parsed_args.build_databases or parsed_args.check_databases:
        return parsed_args
    
    # Validate workflow options
    if not parsed_args.input_file and not parsed_args.read_files and not (parsed_args.pe1_file or parsed_args.pe2_file or parsed_args.single_file):
        parser.error("Must specify input files using --reads, --pe1/--pe2, or --single for viral metagenomic analysis (or use --input for direct annotation of pre-assembled contigs)")
    
    if parsed_args.input_file and (parsed_args.read_files or parsed_args.pe1_file or parsed_args.pe2_file or parsed_args.single_file):
        parser.error("Cannot specify both --input (direct annotation) and read files (metagenomic workflow)")
    
    # Handle explicit paired-end or single-end options
    if parsed_args.pe1_file or parsed_args.pe2_file or parsed_args.single_file:
        if parsed_args.read_files:
            parser.error("Cannot use --reads with --pe1/--pe2 or --single options")
        if parsed_args.pe1_file and parsed_args.pe2_file:
            parsed_args.read_files = [parsed_args.pe1_file, parsed_args.pe2_file]
            parsed_args.paired = True
            parsed_args.single = False
        elif parsed_args.single_file:
            parsed_args.read_files = [parsed_args.single_file]
            parsed_args.paired = False
            parsed_args.single = True
        elif parsed_args.pe1_file or parsed_args.pe2_file:
            parser.error("Both --pe1 and --pe2 must be specified for paired-end reads")
    
    # Auto-detect paired-end vs single-end based on number of files (for --reads option)
    elif parsed_args.read_files:
        if len(parsed_args.read_files) == 1:
            parsed_args.paired = False
            parsed_args.single = True
        elif len(parsed_args.read_files) == 2:
            parsed_args.paired = True
            parsed_args.single = False
            print("[INFO] Using --reads with 2 files. Ensure files are in correct order: R1.fastq R2.fastq")
            print("[INFO] For explicit control, consider using --pe1 R1.fastq --pe2 R2.fastq instead")
        else:
            parser.error("--reads must specify either 1 file (single-end) or 2 files (paired-end)")
    
    # No error here: if user provides read files and no flag, run full pipeline by default
    
    # Only require sample name for metadata generation, otherwise infer from file
    if parsed_args.generate_metadata and not parsed_args.sample_name:
        # Try to infer sample name from reads or input file
        sample_source = None
        if parsed_args.read_files:
            sample_source = Path(parsed_args.read_files[0]).stem
        elif parsed_args.input_file:
            sample_source = Path(parsed_args.input_file).stem
        if sample_source:
            parsed_args.sample_name = sample_source
        else:
            parser.error("--sample-name is required when using --generate-metadata and could not be inferred from input files.")

    # Only require viral_metadata_file for actual analysis/annotation, not for db management or metadata generation
    if not parsed_args.generate_metadata and not (parsed_args.build_databases or parsed_args.check_databases):
        # Only require for annotation/analysis, not for assembly-only or QC-only
        if not parsed_args.viral_metadata_file and not getattr(parsed_args, 'assemble_only', False) and not getattr(parsed_args, 'qc_only', False):
            parser.error("--viral-metadata is required for analysis/annotation runs.")
    
    # Warn if host genome is not provided but moving to assembly or QC
    if (
        (getattr(parsed_args, 'assemble_only', False) or getattr(parsed_args, 'qc_only', False))
        and not getattr(parsed_args, 'skip_host_removal', False)
        and not getattr(parsed_args, 'host_genome', None)
    ):
        print("[WARNING] No host genome provided (--host-genome not set). Proceeding without host removal. This may result in host contamination in the assembly.")
        if hasattr(logger, 'warning'):
            logger.warning("No host genome provided. Proceeding without host removal.")
    
    return parsed_args


def create_config_from_args(args: argparse.Namespace) -> PipelineConfig:
    """
    Create pipeline configuration from command-line arguments.
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Pipeline configuration object
    """
    # Check if this is a workflow that doesn't need databases
    skip_db_check = getattr(args, 'assemble_only', False) or getattr(args, 'qc_only', False)
    
    config = PipelineConfig()
    
    # Map arguments to config attributes
    for arg_name, arg_value in vars(args).items():
        if hasattr(config, arg_name):
            setattr(config, arg_name, arg_value)
    
    # Attach resume flag to config for pipeline usage
    config.resume = getattr(args, 'resume', False)
    
    # Skip database checks for assemble-only and qc-only workflows
    if skip_db_check:
        # Temporarily disable database checks by setting non_crna to True
        # This prevents RFAM database check during config initialization
        config.non_crna = True
        # Also disable HMMER to avoid other database checks
        config.no_hmmer = True
    
    return config


def main(args: Optional[list] = None) -> int:
    """
    Main entry point for the command-line interface.
    """
    try:
        # Parse arguments
        parsed_args = parse_args(args)

        # Handle database management (these don't need output directories)
        if parsed_args.build_databases:
            # Set up simple console logging for database operations
            logger = setup_logger(
                name="viranpy",
                level=logging.INFO,
                console_output=True
            )
            
            from .utils.database_manager import DatabaseManager
            db_manager = DatabaseManager(logger)
            success = db_manager.install_all_databases(fix_missing=parsed_args.fix_missing)
            if success:
                logger.info("All databases installed successfully!")
                return 0
            else:
                logger.error("Database installation failed!")
                return 1
        
        if parsed_args.check_databases:
            # Set up simple console logging for database operations
            logger = setup_logger(
                name="viranpy",
                level=logging.INFO,
                console_output=True
            )
            
            from .utils.database_manager import DatabaseManager
            db_manager = DatabaseManager(logger)
            status = db_manager.check_all_databases()
            if status['all_ready']:
                logger.info("All databases are installed and ready!")
                return 0
            else:
                logger.warning("Some databases are missing or incomplete:")
                for db_name, status_info in status['databases'].items():
                    if not status_info['ready']:
                        logger.warning(f"  - {db_name}: {status_info['message']}")
                logger.info("Run 'viranpy --build-databases' to install missing databases")
                return 1
        
        # Handle viral metadata generation (this doesn't need output directories)
        if parsed_args.generate_metadata:
            # Set up simple console logging for metadata operations
            logger = setup_logger(
                name="viranpy",
                level=logging.INFO,
                console_output=True
            )
            
            if not parsed_args.sample_name:
                logger.error("--sample-name is required when using --generate-metadata")
                return 1
            
            from .utils.metagenomic_metadata import create_viral_metagenomic_metadata
            
            # Handle segmented viruses
            if parsed_args.segments:
                from .utils.metagenomic_metadata import ViralMetagenomicMetadataGenerator
                generator = ViralMetagenomicMetadataGenerator()
                segment_metadata = generator.create_segmented_virus_metadata(
                    sample_name=parsed_args.sample_name,
                    segments=parsed_args.segments,
                    viral_family=parsed_args.viral_family,
                    viral_genus=parsed_args.viral_genus,
                    viral_species=parsed_args.viral_species,
                    collection_date=parsed_args.collection_date,
                    location=parsed_args.location,
                    source=parsed_args.source,
                    host=parsed_args.host
                )
                
                # Create separate metadata files for each segment
                for segment, metadata in segment_metadata.items():
                    output_file = f"viral_metadata_{segment}.txt"
                    with open(output_file, 'w') as f:
                        f.write(metadata + '\n')
                    logger.info(f"Generated metadata file for {segment}: {output_file}")
            else:
                # Generate single viral metadata file
                metadata_file = create_viral_metagenomic_metadata(
                    sample_name=parsed_args.sample_name,
                    output_file="viral_metadata.txt",
                    viral_family=parsed_args.viral_family,
                    viral_genus=parsed_args.viral_genus,
                    viral_species=parsed_args.viral_species,
                    collection_date=parsed_args.collection_date,
                    location=parsed_args.location,
                    source=parsed_args.source,
                    host=parsed_args.host
                )
                
                logger.info(f"Generated viral metadata file: {metadata_file}")
            
            logger.info("You can now use these files with the --viral-metadata option")
            return 0

        # --- Output directory logic (only for analysis/annotation workflows) ---
        # Determine output directory
        output_dir = getattr(parsed_args, 'root_output', None)
        if not output_dir:
            output_dir = 'viranpy_results'
            parsed_args.root_output = output_dir
        output_dir = Path(output_dir)
        orig_output_dir = output_dir
        
        # Handle resume mode differently
        if getattr(parsed_args, 'resume', False):
            # In resume mode, use the existing directory or create it if it doesn't exist
            if output_dir.exists():
                print(f"Resuming pipeline in existing directory: {output_dir}")
                # Check if pipeline state file exists
                state_file = output_dir / "pipeline_state.json"
                if state_file.exists():
                    print(f"Found existing pipeline state file: {state_file}")
                else:
                    print(f"Warning: No pipeline state file found in {output_dir}. This may be a fresh run.")
            else:
                print(f"Creating new directory for resume: {output_dir}")
                output_dir.mkdir(parents=True, exist_ok=True)
        else:
            # Normal mode: handle existing directories
            suffix = 1
            while output_dir.exists():
                # If running interactively, ask user; else, auto-increment
                if sys.stdin.isatty():
                    response = input(f"Output directory '{output_dir}' already exists. Overwrite? [y/N] (or type 'new' to create a new directory): ").strip().lower()
                    if response == 'y':
                        shutil.rmtree(output_dir)
                        break
                    elif response == 'new' or response == 'n' or response == '':
                        output_dir = Path(f"{orig_output_dir}.{suffix}")
                        suffix += 1
                    else:
                        print("Invalid response. Please answer 'y' to overwrite or 'new' to create a new directory.")
                else:
                    # Non-interactive: always create a new suffixed directory
                    output_dir = Path(f"{orig_output_dir}.{suffix}")
                    suffix += 1
            output_dir.mkdir(parents=True, exist_ok=True)
        
        parsed_args.root_output = str(output_dir)

        # --- Logging to output directory ---
        log_file = output_dir / "viranpy.log"
        logger = setup_logger(
            name="viranpy",
            level=logging.INFO,
            log_file=str(log_file),
            console_output=True
        )
        
        # Create configuration
        config = create_config_from_args(parsed_args)
        
        # Create and run pipeline
        pipeline = ViralAnnotationPipeline(config, logger_instance=logger, resume=getattr(config, 'resume', False))
        
        if parsed_args.input_file:
            # Direct annotation of pre-assembled contigs
            logger.info("Running viral annotation pipeline on pre-assembled contigs...")
            annotation_result = pipeline.run(parsed_args.input_file)
            logger.info("Viral annotation pipeline completed successfully")
        elif parsed_args.qc_only:
            # Run quality control only
            pipeline.run_quality_control(parsed_args.read_files, parsed_args.paired)
        elif parsed_args.assemble_only:
            # Run assembly only (without preprocessing steps)
            logger.info("Running assembly pipeline only...")
            
            # Determine input files for assembly
            if parsed_args.read_files:
                assembly_input_files = parsed_args.read_files
                paired = parsed_args.paired
            else:
                logger.error("No input files specified for assembly")
                return 1
            
            # Run assembly directly
            assembly_results = pipeline.run_assembly(assembly_input_files, paired)
            logger.info("Assembly pipeline completed successfully")
        else:
            # Run full metagenomic pipeline: preprocessing + assembly + annotation
            logger.info("Running full metagenomic pipeline: preprocessing, assembly, and viral annotation...")
            preprocessing_results = pipeline.run_preprocessing_pipeline(parsed_args.read_files, parsed_args.paired)
            
            # Run annotation by default after assembly (unless skipped)
            if not getattr(parsed_args, 'skip_annotation', False):
                logger.info("Running viral annotation pipeline on assembled contigs...")
                
                # Find the best assembly result to use for annotation
                best_assembly = None
                if preprocessing_results and 'assembly' in preprocessing_results:
                    assembly_results = preprocessing_results['assembly']
                    # Prefer hybrid, then spades, then megahit
                    for assembler in ['hybrid', 'spades', 'megahit']:
                        if assembler in assembly_results and assembly_results[assembler] and assembly_results[assembler].get('success'):
                            best_assembly = assembly_results[assembler].get('contigs_file')
                            logger.info(f"Using {assembler} assembly for annotation: {best_assembly}")
                            break
                
                if best_assembly and os.path.exists(best_assembly):
                    # Set the input file for annotation
                    config.input_file = best_assembly
                    # Run the annotation pipeline
                    annotation_result = pipeline.run(best_assembly)
                    logger.info("Viral annotation pipeline completed successfully")
                else:
                    logger.warning("No successful assembly found for annotation. Skipping annotation pipeline.")
            else:
                logger.info("Skipping viral annotation pipeline as requested (--skip-annotation)")
        
        logger.info("Pipeline completed successfully")
        print("\nThank you for using ViRAnPy!")
        print("Please cite: Duhan N. et al., ViRAnPy, https://github.com/naveenduhan/viranpy")
        print('"In bioinformatics, the answers you get depend on the questions you ask."')
        return 0
        
    except KeyboardInterrupt:
        if 'logger' in locals():
            logger.info("Pipeline interrupted by user")
        else:
            print("Pipeline interrupted by user")
        return 130
    except Exception as e:
        if 'logger' in locals():
            logger.error(f"Pipeline failed: {e}")
        else:
            print(f"Pipeline failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main()) 