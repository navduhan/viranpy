#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Configuration management for ViRAnPy pipeline.
"""

import os
import multiprocessing
from pathlib import Path
from typing import Optional, Dict, Any
from dataclasses import dataclass, field


@dataclass
class PipelineConfig:
    """
    Configuration class for ViRAnPy pipeline parameters.
    
    This class holds all configuration parameters for the viral metagenomic analysis pipeline,
    including input/output settings, tool parameters, database paths, and preprocessing options.
    """
    
    # Basic options [REQUIRED]
    input_file: Optional[str] = None
    viral_metadata_file: Optional[str] = None
    
    # Advanced general options [OPTIONAL]
    root_output: Optional[str] = None
    locus: str = "LOC"
    ncpus: int = field(default_factory=lambda: multiprocessing.cpu_count())
    min_contig_size: int = 200
    blast_switch: bool = False
    
    # Advanced circularity options
    gc_skew_read_length: int = 101
    gc_skew_window: int = 100
    gc_skew_slide: int = 10
    
    # Advanced CRISPR/repeats options
    min_crispr_repeat: int = 16
    max_crispr_repeat: int = 64
    min_crispr_spacer: int = 8
    max_crispr_spacer: int = 64
    
    # Advanced genetic code options
    use_prodigal_gv: bool = False
    genetic_code_table: str = "11"
    
    # Advanced GenBank division options
    genbank_division: str = "CON"
    
    # Advanced ncRNA options
    non_crna: bool = False
    nor_fam: bool = False
    hmmonly: bool = False
    rfam_database: Optional[str] = None
    
    # Advanced DIAMOND options
    diamond_database: Optional[str] = None
    diamond_evalue: float = 0.00001
    diamond_width_threshold: float = 50.0
    diamond_cov_threshold: float = 50.0
    
    # Advanced BLAST options
    blast_database: Optional[str] = None
    blast_exh: bool = False
    blast_evalue: float = 0.00001
    blast_width_threshold: float = 50.0
    blast_cov_threshold: float = 50.0
    
    # Advanced HMMER options
    no_hmmer: bool = False
    no_vogs: bool = False
    no_rvdb: bool = False
    no_phrogs: bool = False
    no_vfam: bool = False
    vogs_database: Optional[str] = None
    rvdb_database: Optional[str] = None
    phrogs_database: Optional[str] = None
    vfam_database: Optional[str] = None
    hmmer_evalue: float = 0.001
    hmmer_cov_threshold: float = 50.0
    
    # Kraken2 database options
    kraken2_db: Optional[str] = None
    
    def __post_init__(self):
        """Post-initialization processing."""
        # Set default output directory if not provided
        if self.root_output is None:
            self.root_output = 'viranpy_results'
        # Set default database paths
        self._set_default_databases()
    
    def _set_default_databases(self):
        """Set default database paths based on package location."""
        package_dir = Path(__file__).parent
        
        # RFAM database
        if not self.non_crna and self.rfam_database is None:
            rfam_path = package_dir / "databases" / "rfam" / "Rfam.cm"
            if rfam_path.exists():
                self.rfam_database = str(rfam_path)
            else:
                raise FileNotFoundError(
                    "RFAM database not found. Please run 'viranpy --build-databases' "
                    "to install all required databases."
                )
        
        # DIAMOND database
        if not self.blast_switch and self.diamond_database is None:
            diamond_path = package_dir / "databases" / "RefSeq_Viral_DIAMOND" / "refseq_viral_proteins.dmnd"
            if diamond_path.exists():
                self.diamond_database = str(diamond_path)
            else:
                raise FileNotFoundError(
                    "DIAMOND database not found. Please run 'viranpy --build-databases' "
                    "to install all required databases."
                )
        
        # BLAST database
        if self.blast_switch and self.blast_database is None:
            blast_path = package_dir / "databases" / "RefSeq_Viral_BLAST" / "refseq_viral_proteins"
            if blast_path.with_suffix('.pin').exists():
                self.blast_database = str(blast_path)
            else:
                raise FileNotFoundError(
                    "BLAST database not found. Please run 'viranpy --build-databases' "
                    "to install all required databases."
                )
        
        # Kraken2 database
        if self.kraken2_db is None:
            kraken2_path = package_dir / "databases" / "kraken2_viral" / "k2_viral_20250402"
            if kraken2_path.exists():
                self.kraken2_db = str(kraken2_path)
            else:
                raise FileNotFoundError(
                    "Kraken2 database not found. Please run 'viranpy --build-databases' "
                    "to install all required databases."
                )
        
        # HMMER databases
        if not self.no_hmmer:
            self._set_hmmer_databases(package_dir)
    
    def _set_hmmer_databases(self, package_dir: Path):
        """Set HMMER database paths."""
        # VOGS database
        if not self.no_vogs and self.vogs_database is None:
            vogs_path = package_dir / "databases" / "vogs" / "vog_latest.hmm"
            if vogs_path.exists():
                self.vogs_database = str(vogs_path)
            else:
                raise FileNotFoundError(
                    "VOGS database not found. Please run 'viranpy --build-databases' "
                    "to install all required databases."
                )
        
        # RVDB database
        if not self.no_rvdb and self.rvdb_database is None:
            rvdb_path = package_dir / "databases" / "rvdb" / "RVDB_28.0.hmm"
            if rvdb_path.exists():
                self.rvdb_database = str(rvdb_path)
            else:
                raise FileNotFoundError(
                    "RVDB database not found. Please run 'viranpy --build-databases' "
                    "to install all required databases."
                )
        
        # PHROGS database
        if not self.no_phrogs and self.phrogs_database is None:
            phrogs_path = package_dir / "databases" / "phrogs" / "phrogs_v4.hmm"
            if phrogs_path.exists():
                self.phrogs_database = str(phrogs_path)
            else:
                raise FileNotFoundError(
                    "PHROGS database not found. Please run 'viranpy --build-databases' "
                    "to install all required databases."
                )
        
        # VFAM database
        if not self.no_vfam and self.vfam_database is None:
            vfam_path = package_dir / "databases" / "vfam" / "vfam_latest.hmm"
            if vfam_path.exists():
                self.vfam_database = str(vfam_path)
            else:
                raise FileNotFoundError(
                    "VFAM database not found. Please run 'viranpy --build-databases' "
                    "to install all required databases."
                )
    
    def validate(self) -> None:
        """Validate configuration parameters."""
        # Check required files exist only if set
        if self.input_file and not os.path.exists(self.input_file):
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        if self.viral_metadata_file and not os.path.exists(self.viral_metadata_file):
            raise FileNotFoundError(f"Viral metadata file not found: {self.viral_metadata_file}")
        
        # Validate numeric parameters
        if self.ncpus < 1:
            raise ValueError("Number of CPUs must be at least 1")
        
        if self.min_contig_size < 1:
            raise ValueError("Minimum contig size must be at least 1")
        
        # Validate HMMER configuration
        if not self.no_hmmer:
            if all([self.no_vogs, self.no_rvdb, self.no_phrogs, self.no_vfam]):
                raise ValueError("At least one HMMER database must be enabled")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'input_file': self.input_file,
            'viral_metadata_file': self.viral_metadata_file,
            'root_output': self.root_output,
            'locus': self.locus,
            'ncpus': self.ncpus,
            'min_contig_size': self.min_contig_size,
            'blast_switch': self.blast_switch,
            'gc_skew_read_length': self.gc_skew_read_length,
            'gc_skew_window': self.gc_skew_window,
            'gc_skew_slide': self.gc_skew_slide,
            'min_crispr_repeat': self.min_crispr_repeat,
            'max_crispr_repeat': self.max_crispr_repeat,
            'min_crispr_spacer': self.min_crispr_spacer,
            'max_crispr_spacer': self.max_crispr_spacer,
            'use_prodigal_gv': self.use_prodigal_gv,
            'genetic_code_table': self.genetic_code_table,
            'genbank_division': self.genbank_division,
            'non_crna': self.non_crna,
            'nor_fam': self.nor_fam,
            'hmmonly': self.hmmonly,
            'rfam_database': self.rfam_database,
            'diamond_database': self.diamond_database,
            'diamond_evalue': self.diamond_evalue,
            'diamond_width_threshold': self.diamond_width_threshold,
            'diamond_cov_threshold': self.diamond_cov_threshold,
            'blast_database': self.blast_database,
            'blast_exh': self.blast_exh,
            'blast_evalue': self.blast_evalue,
            'blast_width_threshold': self.blast_width_threshold,
            'blast_cov_threshold': self.blast_cov_threshold,
            'no_hmmer': self.no_hmmer,
            'no_vogs': self.no_vogs,
            'no_rvdb': self.no_rvdb,
            'no_phrogs': self.no_phrogs,
            'no_vfam': self.no_vfam,
            'vogs_database': self.vogs_database,
            'rvdb_database': self.rvdb_database,
            'phrogs_database': self.phrogs_database,
            'vfam_database': self.vfam_database,
            'hmmer_evalue': self.hmmer_evalue,
            'hmmer_cov_threshold': self.hmmer_cov_threshold,
            'kraken2_db': self.kraken2_db,
        } 