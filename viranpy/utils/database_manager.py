#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Database management for ViRAnPy pipeline.
"""

import os
import sys
import subprocess
import logging
import shutil
import urllib.request
import tarfile
import gzip
from pathlib import Path
from typing import Dict, Any, List, Optional
import tempfile
from tqdm import tqdm


class DatabaseManager:
    """
    Manages installation and checking of all required databases for ViRAnPy.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        Initialize the database manager.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger or logging.getLogger(__name__)
        self.package_dir = Path(__file__).parent.parent
        self.databases_dir = self.package_dir / "databases"
        
        # Database URLs and configurations
        self.databases = {
            'rfam': {
                'url': 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz',
                'filename': 'Rfam.cm.gz',
                'extracted': 'Rfam.cm',
                'final': 'Rfam.cm',
                'type': 'infernal',
                'required_tools': ['cmpress'],
                'description': 'RFAM database for ncRNA detection'
            },
            'refseq_viral': {
                'url': 'https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz',
                'filename': 'viral.1.protein.faa.gz',
                'extracted': 'viral.1.protein.faa',
                'final': 'refseq_viral_proteins.faa',
                'type': 'protein',
                'required_tools': ['diamond', 'makeblastdb'],
                'description': 'RefSeq Viral Proteins database'
            },
            'kraken2_viral': {
                'url': 'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20250402.tar.gz',
                'filename': 'k2_viral_20250402.tar.gz',
                'extracted': 'k2_viral_20250402',
                'final': 'k2_viral_20250402',
                'type': 'kraken2',
                'required_tools': ['kraken2'],
                'description': 'Kraken2 viral-only database for taxonomic classification',
                'archive_size_gb': 0.5,
                'index_size_gb': 0.6
            },
            'krona_taxonomy': {
                'url': None,  # No download needed since Krona is installed via conda
                'filename': None,
                'extracted': None,
                'final': 'taxonomy',
                'type': 'krona',
                'required_tools': ['ktUpdateTaxonomy.sh'],
                'description': 'Krona taxonomy database for visualization',
                'conda_installed': True
            },
            'rvdb': {
                'url': 'https://rvdb-prot.pasteur.fr/files/U-RVDBv28.0-prot.hmm.xz',
                'filename': 'U-RVDBv28.0-prot.hmm.xz',
                'extracted': 'U-RVDBv28.0-prot.hmm',
                'final': 'RVDB_28.0.hmm',
                'type': 'hmmer',
                'required_tools': ['hmmpress'],
                'description': 'RVDB database for viral protein families'
            },
            'vogs': {
                'url': 'https://fileshare.lisc.univie.ac.at/vog/vog225/vog.hmm.tar.gz',
                'filename': 'vog.hmm.tar.gz',
                'extracted': 'vog_latest.hmm',
                'final': 'vog_latest.hmm',
                'type': 'hmmer',
                'required_tools': ['hmmpress'],
                'description': 'VOGS database for viral orthologous groups'
            },
            'vfam': {
                'url': 'https://fileshare.lisc.univie.ac.at/vog/vog225/vfam.hmm.tar.gz',
                'filename': 'vfam.hmm.tar.gz',
                'extracted': 'vfam_latest.hmm',
                'final': 'vfam_latest.hmm',
                'type': 'hmmer',
                'required_tools': ['hmmpress'],
                'description': 'VFAM database for viral protein families'
            },
            'phrogs': {
                'url': 'https://phrogs.lmge.uca.fr/downloads_from_website/MSA_phrogs.tar.gz',
                'filename': 'MSA_phrogs.tar.gz',
                'extracted': 'phrogs_v4.hmm',
                'final': 'phrogs_v4.hmm',
                'type': 'hmmer',
                'required_tools': ['esl-reformat', 'hmmbuild', 'hmmpress'],
                'description': 'PHROGS database for phage protein families'
            }
        }
    
    def check_required_tools(self) -> Dict[str, bool]:
        """
        Check if all required tools are available.
        
        Returns:
            Dictionary mapping tool names to availability status
        """
        tools_status = {}
        required_tools = set()
        
        # Collect all required tools
        for db_info in self.databases.values():
            required_tools.update(db_info['required_tools'])
        
        # Check each tool
        for tool in required_tools:
            tools_status[tool] = shutil.which(tool) is not None
        
        return tools_status
    
    def install_all_databases(self, fix_missing: bool = False) -> bool:
        """
        Install all required databases.
        
        Returns:
            True if all databases were installed successfully
        """
        self.logger.info("ViRAnPy Database Installation")
        self.logger.info("=" * 50)
        self.logger.info("Database installation may take 1-1.5 hours on an 8-core computer with a fast internet connection. This process is resource- and network-intensive. Sit back, relax, and let ViRAnPy handle everything automatically!")
        
        # Check required tools
        tools_status = self.check_required_tools()
        missing_tools = [tool for tool, available in tools_status.items() if not available]
        
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            self.logger.error("Please install these tools first using conda or your system package manager:")
            self.logger.error("conda install -c bioconda infernal diamond blast hmmer kraken2 krona")
            return False
        
        # Create databases directory
        self.databases_dir.mkdir(exist_ok=True)
        
        success_count = 0
        total_count = len(self.databases)
        optional_count = 0
        
        for db_name, db_info in self.databases.items():
            # Check if database is already installed and ready
            db_status = self._check_database(db_name, db_info, fix_missing=fix_missing)
            if db_status.get('ready', False):
                self.logger.info(f"{db_name.upper()} database already installed. Skipping.")
                success_count += 1
                continue
            self.logger.info(f"\nInstalling {db_name.upper()} database...")
            self.logger.info(f"Description: {db_info['description']}")
            
            # Show size information for large databases
            if 'archive_size_gb' in db_info:
                self.logger.info(f"Archive size: {db_info['archive_size_gb']} GB")
                self.logger.info(f"Index size: {db_info['index_size_gb']} GB")
            
            # Skip optional databases if tools are not available
            if db_info.get('optional', False):
                optional_count += 1
                missing_optional_tools = [tool for tool in db_info['required_tools'] 
                                        if not tools_status.get(tool, False)]
                if missing_optional_tools:
                    self.logger.warning(f"Skipping {db_name} (optional): Missing tools: {', '.join(missing_optional_tools)}")
                    continue
            
            try:
                if self._install_database(db_name, db_info):
                    success_count += 1
                    self.logger.info(f"✓ {db_name.upper()} database installed successfully")
                else:
                    self.logger.error(f"✗ Failed to install {db_name.upper()} database")
            except Exception as e:
                self.logger.error(f"✗ Error installing {db_name.upper()} database: {e}")
        
        self.logger.info(f"\nDatabase installation summary:")
        self.logger.info(f"Successfully installed: {success_count}/{total_count - optional_count} required databases")
        if optional_count > 0:
            self.logger.info(f"Skipped {optional_count} optional databases")
        
        if success_count == total_count - optional_count:
            self.logger.info("All required databases installed successfully!")
            return True
        else:
            self.logger.warning("Some databases failed to install. Check the logs above.")
            return False
    
    def _install_database(self, db_name: str, db_info: Dict[str, Any]) -> bool:
        """
        Install a specific database.
        
        Args:
            db_name: Database name
            db_info: Database configuration
            
        Returns:
            True if installation was successful
        """
        # For krona_taxonomy and refseq_viral, do not create or use a directory
        if db_name == 'krona_taxonomy':
            return self._install_krona_taxonomy(db_info)
        if db_name == 'refseq_viral':
            return self._install_refseq_viral(db_info)
        db_dir = self.databases_dir / db_name
        db_dir.mkdir(exist_ok=True)
        
        # Handle conda-installed tools (no download needed)
        if db_info.get('conda_installed', False):
            return self._install_conda_tool(db_name, db_info)
        
        # Download database
        if not self._download_file(db_info['url'], db_dir / db_info['filename']):
            return False
        
        # Extract and process database
        if db_name == 'rfam':
            return self._install_rfam(db_dir, db_info)
        elif db_name == 'kraken2_viral':
            return self._install_kraken2_viral(db_dir, db_info)
        elif db_name == 'rvdb':
            return self._install_rvdb(db_dir, db_info)
        elif db_name == 'vogs':
            return self._install_vogs(db_dir, db_info)
        elif db_name == 'vfam':
            return self._install_vfam(db_dir, db_info)
        elif db_name == 'phrogs':
            return self._install_phrogs(db_dir, db_info)
        else:
            self.logger.error(f"Unknown database: {db_name}")
            return False
    
    def _install_conda_tool(self, db_name: str, db_info: Dict[str, Any]) -> bool:
        """
        Install database for conda-installed tools.
        
        Args:
            db_name: Database name
            db_info: Database configuration
            
        Returns:
            True if installation was successful
        """
        if db_name == 'krona_taxonomy':
            return self._install_krona_taxonomy(db_info)
        else:
            self.logger.error(f"Unknown conda-installed tool: {db_name}")
            return False
    
    def _download_file(self, url: str, output_path: Path) -> bool:
        """
        Download a file from URL with a progress bar.
        """
        try:
            self.logger.info(f"Downloading {url}...")
            with tqdm(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=output_path.name) as t:
                def reporthook(block_num, block_size, total_size):
                    if total_size > 0:
                        t.total = total_size
                    t.update(block_size)
                urllib.request.urlretrieve(url, output_path, reporthook)
            return True
        except Exception as e:
            self.logger.error(f"Failed to download {url}: {e}")
            return False
    
    def _install_rfam(self, db_dir: Path, db_info: Dict[str, Any]) -> bool:
        """Install RFAM database."""
        try:
            # Extract gzipped file
            with gzip.open(db_dir / db_info['filename'], 'rb') as f_in:
                with open(db_dir / db_info['extracted'], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Run cmpress
            self.logger.info(f"Running: cmpress {db_dir / db_info['extracted']}")
            subprocess.run(['cmpress', str(db_dir / db_info['extracted'])], 
                         cwd=db_dir, check=True, capture_output=True)
            
            # Clean up
            (db_dir / db_info['filename']).unlink()
            return True
        except Exception as e:
            self.logger.error(f"Error installing RFAM: {e}")
            return False
    
    def _install_refseq_viral(self, db_info: Dict[str, Any]) -> bool:
        """Install RefSeq Viral database (no directory needed)."""
        import tempfile
        try:
            # Download and extract to a temporary file
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir_path = Path(tmpdir)
                gz_path = tmpdir_path / db_info['filename']
                extracted_path = tmpdir_path / db_info['extracted']
                # Download
                self.logger.info(f"Downloading {db_info['url']}...")
                if not self._download_file(db_info['url'], gz_path):
                    return False
                self.logger.info(f"Extracting {db_info['filename']}...")
                # Extract gzipped file with progress bar
                total_size = os.path.getsize(gz_path)
                with gzip.open(gz_path, 'rb') as f_in, open(extracted_path, 'wb') as f_out:
                    with tqdm(total=total_size, unit='B', unit_scale=True, desc='Extracting') as pbar:
                        while True:
                            chunk = f_in.read(1024 * 1024)
                            if not chunk:
                                break
                            f_out.write(chunk)
                            pbar.update(len(chunk))
                # Create DIAMOND database
                diamond_dir = self.databases_dir / "RefSeq_Viral_DIAMOND"
                diamond_dir.mkdir(exist_ok=True)
                self.logger.info("Building DIAMOND database...")
                shutil.copy(extracted_path, diamond_dir / db_info['final'])
                subprocess.run(['diamond', 'makedb', '--in', str(diamond_dir / db_info['final']), 
                              '-d', str(diamond_dir / 'refseq_viral_proteins')], 
                             cwd=diamond_dir, check=True, capture_output=True)
                (diamond_dir / db_info['final']).unlink()
                # Create BLAST database
                blast_dir = self.databases_dir / "RefSeq_Viral_BLAST"
                blast_dir.mkdir(exist_ok=True)
                self.logger.info("Building BLAST database...")
                shutil.copy(extracted_path, blast_dir / db_info['final'])
                subprocess.run(['makeblastdb', '-in', str(blast_dir / db_info['final']), 
                              '-dbtype', 'prot', '-out', str(blast_dir / 'refseq_viral_proteins')], 
                             cwd=blast_dir, check=True, capture_output=True)
                (blast_dir / db_info['final']).unlink()
            return True
        except Exception as e:
            self.logger.error(f"Error installing RefSeq Viral: {e}")
            return False
    
    def _install_kraken2_viral(self, db_dir: Path, db_info: Dict[str, Any]) -> bool:
        """Install Kraken2 viral database."""
        try:
            # Always create the expected directory
            extracted_dir = db_dir / db_info['extracted']
            extracted_dir.mkdir(exist_ok=True)
            tar_path = extracted_dir / db_info['filename']
            # Download tarball inside the directory
            if not self._download_file(db_info['url'], tar_path):
                return False
            self.logger.info(f"Extracting: {db_info['filename']} inside {extracted_dir}")
            with tarfile.open(tar_path, 'r:gz') as tar:
                tar.extractall(extracted_dir)
            tar_path.unlink(missing_ok=True)
            # Move any files in db_dir (but not in extracted_dir) that match key files into extracted_dir
            key_files = [
                'hash.k2d', 'library_report.tsv', 'names.dmp', 'nodes.dmp', 'seqid2taxid.map', 'taxo.k2d', 'ktaxonomy.tsv'
            ]
            for f in key_files:
                src = db_dir / f
                dst = extracted_dir / f
                if src.exists() and not dst.exists():
                    src.rename(dst)
            return True
        except Exception as e:
            self.logger.error(f"Error installing Kraken2 viral: {e}")
            return False
    
    def _install_krona_taxonomy(self, db_info: Dict[str, Any]) -> bool:
        """Install Krona taxonomy database (no directory needed)."""
        try:
            # Run ktUpdateTaxonomy to download taxonomy data
            self.logger.info("Downloading Krona taxonomy data...")
            
            # Run ktUpdateTaxonomy.sh from anywhere (it will download to Krona's default location)
            self.logger.info("Running: ktUpdateTaxonomy.sh")
            subprocess.run(['ktUpdateTaxonomy.sh'], 
                         check=True, capture_output=True)
            
            self.logger.info("Krona taxonomy data downloaded successfully")
            return True
        except Exception as e:
            self.logger.error(f"Error installing Krona taxonomy: {e}")
            return False
    
    def _install_rvdb(self, db_dir: Path, db_info: Dict[str, Any]) -> bool:
        """Install RVDB database."""
        try:
            # Extract xz file
            self.logger.info(f"Running: unxz {db_info['filename']}")
            subprocess.run(['unxz', str(db_dir / db_info['filename'])], 
                         cwd=db_dir, check=True, capture_output=True)
            
            # Rename to final name
            (db_dir / db_info['extracted']).rename(db_dir / db_info['final'])
            
            # Run hmmpress
            self.logger.info(f"Running: hmmpress {db_info['final']}")
            subprocess.run(['hmmpress', str(db_dir / db_info['final'])], 
                         cwd=db_dir, check=True, capture_output=True)
            
            # Clean up
            (db_dir / db_info['filename']).unlink(missing_ok=True)
            return True
        except Exception as e:
            self.logger.error(f"Error installing RVDB: {e}")
            return False
    
    def _install_vogs(self, db_dir: Path, db_info: Dict[str, Any]) -> bool:
        """Install VOGS database."""
        try:
            # Extract tar.gz
            self.logger.info(f"Extracting: {db_info['filename']}")
            with tarfile.open(db_dir / db_info['filename'], 'r:gz') as tar:
                tar.extractall(db_dir)
            
            # Combine HMM files
            hmm_files = list(db_dir.glob('hmm/*.hmm'))
            self.logger.info(f"Combining {len(hmm_files)} HMM files for VOGS")
            with open(db_dir / db_info['final'], 'w') as outfile:
                for hmm_file in hmm_files:
                    with open(hmm_file, 'r') as infile:
                        outfile.write(infile.read())
            
            # Run hmmpress
            self.logger.info(f"Running: hmmpress {db_info['final']}")
            subprocess.run(['hmmpress', str(db_dir / db_info['final'])], 
                         cwd=db_dir, check=True, capture_output=True)
            
            # Clean up
            (db_dir / db_info['filename']).unlink()
            shutil.rmtree(db_dir / 'hmm', ignore_errors=True)
            return True
        except Exception as e:
            self.logger.error(f"Error installing VOGS: {e}")
            return False
    
    def _install_vfam(self, db_dir: Path, db_info: Dict[str, Any]) -> bool:
        """Install VFAM database."""
        try:
            # Extract tar.gz
            self.logger.info(f"Extracting: {db_info['filename']}")
            with tarfile.open(db_dir / db_info['filename'], 'r:gz') as tar:
                tar.extractall(db_dir)
            
            # Combine HMM files
            hmm_files = list(db_dir.glob('hmm/*.hmm'))
            self.logger.info(f"Combining {len(hmm_files)} HMM files for VFAM")
            with open(db_dir / db_info['final'], 'w') as outfile:
                for hmm_file in hmm_files:
                    with open(hmm_file, 'r') as infile:
                        outfile.write(infile.read())
            
            # Run hmmpress
            self.logger.info(f"Running: hmmpress {db_info['final']}")
            subprocess.run(['hmmpress', str(db_dir / db_info['final'])], 
                         cwd=db_dir, check=True, capture_output=True)
            
            # Clean up
            (db_dir / db_info['filename']).unlink()
            shutil.rmtree(db_dir / 'hmm', ignore_errors=True)
            return True
        except Exception as e:
            self.logger.error(f"Error installing VFAM: {e}")
            return False
    
    def _install_phrogs(self, db_dir: Path, db_info: Dict[str, Any]) -> bool:
        """Install PHROGS database."""
        try:
            # Extract tar.gz
            self.logger.info(f"Extracting: {db_info['filename']}")
            with tarfile.open(db_dir / db_info['filename'], 'r:gz') as tar:
                tar.extractall(db_dir)
            
            # Process MSA files
            msa_dir = db_dir / 'MSA_Phrogs_M50_FASTA'
            if msa_dir.exists():
                hmm_files_to_build = list(msa_dir.glob('*.fma'))
                self.logger.info(f"Building HMMs for PHROGS: {len(hmm_files_to_build)} segments")
                for fma_file in hmm_files_to_build:
                    base_name = fma_file.stem
                    sto_file = msa_dir / f"{base_name}.sto"
                    hmm_file = msa_dir / f"{base_name}.hmm"
                    
                    # Convert to Stockholm format
                    self.logger.info(f"Running: esl-reformat -o {sto_file} stockholm {fma_file}")
                    subprocess.run(['esl-reformat', '-o', str(sto_file), 'stockholm', str(fma_file)], 
                                 cwd=msa_dir, check=True, capture_output=True)
                    
                    # Build HMM
                    # Only log once before the loop, not per file
                    subprocess.run(['hmmbuild', '--cpu', '8', str(hmm_file), str(sto_file)], 
                                 cwd=msa_dir, check=True, capture_output=True)
                
                # Combine HMM files
                hmm_files = list(msa_dir.glob('*.hmm'))
                with open(db_dir / db_info['final'], 'w') as outfile:
                    for hmm_file in hmm_files:
                        with open(hmm_file, 'r') as infile:
                            outfile.write(infile.read())
                
                # Run hmmpress
                self.logger.info(f"Running: hmmpress {db_info['final']}")
                subprocess.run(['hmmpress', str(db_dir / db_info['final'])], 
                             cwd=db_dir, check=True, capture_output=True)
                
                # Clean up
                shutil.rmtree(msa_dir, ignore_errors=True)
            
            # Clean up
            (db_dir / db_info['filename']).unlink()
            return True
        except Exception as e:
            self.logger.error(f"Error installing PHROGS: {e}")
            return False
    
    def check_all_databases(self) -> Dict[str, Any]:
        """
        Check if all databases are installed and ready.
        
        Returns:
            Dictionary with database status information
        """
        status = {
            'all_ready': True,
            'databases': {}
        }
        
        for db_name, db_info in self.databases.items():
            db_status = self._check_database(db_name, db_info)
            status['databases'][db_name] = db_status
            
            if not db_status['ready']:
                status['all_ready'] = False
        
        return status
    
    def _check_database(self, db_name: str, db_info: Dict[str, Any], fix_missing: bool = False) -> Dict[str, Any]:
        """
        Check if a specific database is installed and ready.
        
        Args:
            db_name: Database name
            db_info: Database configuration
            fix_missing: If True, attempt to repair missing files
        Returns:
            Dictionary with database status
        """
        db_dir = self.databases_dir / db_name
        
        if not db_dir.exists():
            return {
                'ready': False,
                'message': 'Database directory not found'
            }
        
        # Check specific database files
        if db_name == 'rfam':
            return self._check_rfam(db_dir, db_info, fix_missing=fix_missing)
        elif db_name == 'refseq_viral':
            return self._check_refseq_viral(db_dir, db_info, fix_missing=fix_missing)
        elif db_name == 'kraken2_viral':
            return self._check_kraken2_viral(db_dir, db_info, fix_missing=fix_missing)
        elif db_name == 'krona_taxonomy':
            return self._check_krona_taxonomy(db_info)
        elif db_name in ['rvdb', 'vogs', 'vfam', 'phrogs']:
            return self._check_hmmer_db(db_dir, db_info, fix_missing=fix_missing)
        else:
            return {
                'ready': False,
                'message': f'Unknown database type: {db_name}'
            }
    
    def _check_rfam(self, db_dir: Path, db_info: Dict[str, Any], fix_missing: bool = False) -> Dict[str, Any]:
        """Check RFAM database. Attempt to create missing index files if possible."""
        required_files = [
            db_info['final'],
            f"{db_info['final']}.i1m",
            f"{db_info['final']}.i1p",
            f"{db_info['final']}.i1f",
            f"{db_info['final']}.i1i"
        ]
        present_files = [f for f in required_files if (db_dir / f).exists()]
        missing_files = [f for f in required_files if not (db_dir / f).exists()]
        if (db_dir / db_info['final']).exists() and any((db_dir / f).exists() for f in required_files[1:]):
            if missing_files:
                self.logger.warning(f"RFAM: Some index files are missing: {', '.join(missing_files)}. Attempting to create them with cmpress." if fix_missing else f"RFAM: Some index files are missing: {', '.join(missing_files)}. Functionality may be limited.")
                if fix_missing:
                    # Delete any existing Rfam.cm.i1* files before running cmpress
                    for ext in ['i1m', 'i1p', 'i1f', 'i1i']:
                        idx_file = db_dir / f"{db_info['final']}.{ext}"
                        if idx_file.exists():
                            self.logger.info(f"RFAM: Deleting old index file {idx_file}")
                            idx_file.unlink()
                    try:
                        subprocess.run(['cmpress', db_info['final']], cwd=db_dir, check=True, capture_output=True)
                        # Re-check for missing files
                        missing_files = [f for f in required_files if not (db_dir / f).exists()]
                        if not missing_files:
                            self.logger.info("RFAM: All index files created successfully.")
                            return {
                                'ready': True,
                                'message': 'Database ready (index files created)'
                            }
                        else:
                            self.logger.warning(f"RFAM: Still missing index files after cmpress: {', '.join(missing_files)}. Functionality may be limited.")
                    except Exception as e:
                        self.logger.error(f"RFAM: Failed to run cmpress to create index files: {e}")
            return {
                'ready': True,
                'message': 'Database ready (some index files missing)' if missing_files else 'Database ready'
            }
        else:
            return {
                'ready': False,
                'message': f"Missing files: {', '.join(missing_files)}"
            }
    
    def _check_refseq_viral(self, db_dir: Path, db_info: Dict[str, Any], fix_missing: bool = False) -> Dict[str, Any]:
        """Check RefSeq Viral database."""
        diamond_dir = self.databases_dir / "RefSeq_Viral_DIAMOND"
        blast_dir = self.databases_dir / "RefSeq_Viral_BLAST"
        # Check if both directories exist
        if not diamond_dir.exists() and not blast_dir.exists():
            return {
                'ready': False,
                'message': 'Both RefSeq_Viral_DIAMOND and RefSeq_Viral_BLAST directories not found'
            }
        if not diamond_dir.exists():
            return {
                'ready': False,
                'message': 'RefSeq_Viral_DIAMOND directory not found'
            }
        if not blast_dir.exists():
            return {
                'ready': False,
                'message': 'RefSeq_Viral_BLAST directory not found'
            }
        diamond_db = diamond_dir / "refseq_viral_proteins.dmnd"
        essential_blast_files = [
            "refseq_viral_proteins.psq",
            "refseq_viral_proteins.phr",
            "refseq_viral_proteins.pin"
        ]
        optional_blast_files = [
            "refseq_viral_proteins.pdb",
            "refseq_viral_proteins.pog",
            "refseq_viral_proteins.psd",
            "refseq_viral_proteins.psi"
        ]
        missing_essential = [f for f in essential_blast_files if not (blast_dir / f).exists()]
        missing_optional = [f for f in optional_blast_files if not (blast_dir / f).exists()]
        if not diamond_db.exists():
            if fix_missing and (blast_dir / "refseq_viral_proteins.faa").exists():
                self.logger.info("RefSeq Viral: Attempting to create DIAMOND database with diamond makedb.")
                try:
                    subprocess.run(['diamond', 'makedb', '--in', 'refseq_viral_proteins.faa', '-d', 'refseq_viral_proteins'], cwd=diamond_dir, check=True, capture_output=True)
                except Exception as e:
                    self.logger.error(f"RefSeq Viral: Failed to create DIAMOND database: {e}")
                if diamond_db.exists():
                    self.logger.info("RefSeq Viral: DIAMOND database created successfully.")
                else:
                    return {
                        'ready': False,
                        'message': 'DIAMOND database not found and could not be created.'
                    }
            else:
                return {
                    'ready': False,
                    'message': 'DIAMOND database not found'
                }
        if missing_essential:
            if fix_missing and (blast_dir / "refseq_viral_proteins.faa").exists():
                self.logger.info("RefSeq Viral: Attempting to create BLAST database with makeblastdb.")
                try:
                    subprocess.run(['makeblastdb', '-in', 'refseq_viral_proteins.faa', '-dbtype', 'prot', '-out', 'refseq_viral_proteins'], cwd=blast_dir, check=True, capture_output=True)
                except Exception as e:
                    self.logger.error(f"RefSeq Viral: Failed to create BLAST database: {e}")
                missing_essential = [f for f in essential_blast_files if not (blast_dir / f).exists()]
                if not missing_essential:
                    self.logger.info("RefSeq Viral: BLAST database created successfully.")
                else:
                    return {
                        'ready': False,
                        'message': f"Missing essential BLAST files: {', '.join(missing_essential)}"
                    }
            else:
                return {
                    'ready': False,
                    'message': f"Missing essential BLAST files: {', '.join(missing_essential)}"
                }
        if missing_optional and fix_missing and (blast_dir / "refseq_viral_proteins.faa").exists():
            self.logger.info("RefSeq Viral: Attempting to regenerate optional BLAST files with makeblastdb.")
            try:
                subprocess.run(['makeblastdb', '-in', 'refseq_viral_proteins.faa', '-dbtype', 'prot', '-out', 'refseq_viral_proteins'], cwd=blast_dir, check=True, capture_output=True)
            except Exception as e:
                self.logger.error(f"RefSeq Viral: Failed to regenerate optional BLAST files: {e}")
            missing_optional = [f for f in optional_blast_files if not (blast_dir / f).exists()]
            if not missing_optional:
                self.logger.info("RefSeq Viral: All optional BLAST files created successfully.")
            else:
                self.logger.warning(f"RefSeq Viral BLAST: Still missing optional files after makeblastdb: {', '.join(missing_optional)}. Most functionality will work.")
        return {
            'ready': True,
            'message': 'Database ready (some optional BLAST files missing)' if missing_optional else 'Database ready'
        }
    
    def _check_kraken2_viral(self, db_dir: Path, db_info: Dict[str, Any], fix_missing: bool = False) -> Dict[str, Any]:
        """Check Kraken2 viral database."""
        expected_dir = db_dir / db_info['extracted']
        key_files = [
            'hash.k2d', 'library_report.tsv', 'names.dmp', 'nodes.dmp', 'seqid2taxid.map', 'taxo.k2d', 'ktaxonomy.tsv'
        ]
        if expected_dir.exists():
            return {
                'ready': True,
                'message': 'Database ready'
            }
        # If directory is missing, check for key files in db_dir
        if all((db_dir / f).exists() for f in key_files):
            if fix_missing:
                self.logger.info(f"Kraken2 viral: Moving key files into {expected_dir} to fix missing directory.")
                expected_dir.mkdir(exist_ok=True)
                for f in key_files:
                    src = db_dir / f
                    dst = expected_dir / f
                    if src.exists() and not dst.exists():
                        src.rename(dst)
                if all((expected_dir / f).exists() for f in key_files):
                    self.logger.info(f"Kraken2 viral: Directory {expected_dir} created and key files moved.")
                    return {
                        'ready': True,
                        'message': f"Database ready (directory '{db_info['extracted']}' created and key files moved)"
                    }
            self.logger.warning(f"Kraken2 viral: Expected directory '{db_info['extracted']}' not found, but key files are present in {db_dir}. Functionality should work.")
            return {
                'ready': True,
                'message': f"Database ready (directory '{db_info['extracted']}' missing, but key files found)"
            }
        return {
            'ready': False,
            'message': f"Kraken2 viral database not found (missing directory '{db_info['extracted']}' and/or key files)"
        }
    
    def _check_krona_taxonomy(self, db_info: Dict[str, Any]) -> Dict[str, Any]:
        """Check Krona taxonomy database."""
        # Check if ktUpdateTaxonomy.sh is available
        if not shutil.which('ktUpdateTaxonomy.sh'):
            return {
                'ready': False,
                'message': 'ktUpdateTaxonomy.sh not found (Krona not installed)'
            }
        
        # Check if taxonomy data exists (Krona stores it in its installation directory)
        try:
            # Try to run ktUpdateTaxonomy.sh to check if taxonomy data is available
            result = subprocess.run(['ktUpdateTaxonomy.sh', '--help'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                return {
                    'ready': True,
                    'message': 'Krona taxonomy available'
                }
            else:
                return {
                    'ready': False,
                    'message': 'Krona taxonomy not properly configured'
                }
        except Exception as e:
            return {
                'ready': False,
                'message': f'Error checking Krona taxonomy: {e}'
            }
    
    def _check_hmmer_db(self, db_dir: Path, db_info: Dict[str, Any], fix_missing: bool = False) -> Dict[str, Any]:
        """Check HMMER database. Attempt to create missing index files if possible."""
        required_files = [
            db_info['final'],
            f"{db_info['final']}.h3f",
            f"{db_info['final']}.h3i",
            f"{db_info['final']}.h3m",
            f"{db_info['final']}.h3p"
        ]
        missing_files = [f for f in required_files if not (db_dir / f).exists()]
        if (db_dir / db_info['final']).exists() and any((db_dir / f).exists() for f in required_files[1:]):
            if missing_files:
                self.logger.warning(f"HMMER DB: Some index files are missing: {', '.join(missing_files)}. Attempting to create them with hmmpress." if fix_missing else f"HMMER DB: Some index files are missing: {', '.join(missing_files)}. Functionality may be limited.")
                if fix_missing:
                    try:
                        subprocess.run(['hmmpress', db_info['final']], cwd=db_dir, check=True, capture_output=True)
                        missing_files = [f for f in required_files if not (db_dir / f).exists()]
                        if not missing_files:
                            self.logger.info("HMMER DB: All index files created successfully.")
                            return {
                                'ready': True,
                                'message': 'Database ready (index files created)'
                            }
                        else:
                            self.logger.warning(f"HMMER DB: Still missing index files after hmmpress: {', '.join(missing_files)}. Functionality may be limited.")
                    except Exception as e:
                        self.logger.error(f"HMMER DB: Failed to run hmmpress to create index files: {e}")
            return {
                'ready': True,
                'message': 'Database ready (some index files missing)' if missing_files else 'Database ready'
            }
        else:
            return {
                'ready': False,
                'message': f"Missing files: {', '.join(missing_files)}"
            } 