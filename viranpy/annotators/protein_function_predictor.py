#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan
"""
Protein function annotation in viral genomes using BLAST, DIAMOND, and HMMER.
"""

import os
import re
import csv
import subprocess
from typing import Dict, Any, List
from pathlib import Path
from collections import defaultdict, Counter, namedtuple
import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from ..core.base import BaseAnnotator
from ..core.results import AnnotationResult
from ..utils.file_utils import cmd_exists, safe_run_cmd

Result = namedtuple('Result', ['protein', 'query_name', 'database', 'pcoverage', 'evalue'])

class ViralProteinAnnotator(BaseAnnotator):
    """
    Annotate viral protein function using BLAST, DIAMOND, and HMMER.
    Unique ViRAnPy implementation.
    """
    def __init__(self, config, logger=None):
        super().__init__(config, logger)
        self.protein_annotations = {}
        self.protein_results = {}
        self.hmmer_annotations = {}
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        return cmd_exists("diamond") or cmd_exists("blastp")
    
    def validate_input(self, input_file: str) -> bool:
        """Validate input FASTA file."""
        return Path(input_file).exists() and Path(input_file).suffix.lower() in ['.fasta', '.fa', '.faa']
    
    def run(self, input_file: str, **kwargs) -> Dict[str, Any]:
        """
        Run protein function annotation for all sequences in the input file.
        Returns a dictionary with annotation results.
        """
        result = AnnotationResult(
            annotator_name=self.name,
            input_file=input_file
        )
        try:
            # Check if we have a protein file from gene prediction
            protein_file = self._get_protein_file(input_file)
            if not protein_file:
                self.logger.warning("No protein file found from gene prediction. Skipping protein function prediction.")
                result.success = False
                result.error_message = "No protein sequences available for annotation"
                return result.to_dict()
            
            self._run_homology_search(protein_file)
            self._process_proteins(protein_file)
            if not getattr(self.config, 'no_hmmer', False):
                self._run_hmmer_annotation(protein_file)
            result.add_annotation("protein_results", self.protein_results)
            result.add_annotation("homology_annotations", self.protein_annotations)
            result.add_annotation("hmmer_annotations", self.hmmer_annotations)
            result.add_metadata("total_proteins", sum(len(v) for v in self.protein_results.values()))
            result.add_metadata("total_hmmer_domains", len(self.hmmer_annotations))
            self.logger.info(f"Protein function annotation completed: {result.metadata['total_proteins']} proteins annotated, {result.metadata['total_hmmer_domains']} HMMER domains found")
        except Exception as e:
            result.success = False
            result.error_message = str(e)
            self.logger.error(f"Protein function annotation failed: {e}")
        return result.to_dict()
    
    def _get_protein_file(self, input_file: str) -> str:
        """Get the protein file from gene prediction if it exists."""
        # Check for protein file in gene prediction directory
        gene_pred_dir = Path(self.config.root_output) / "gene_prediction"
        protein_file = gene_pred_dir / f"{Path(input_file).stem}_proteins.faa"
        
        if protein_file.exists():
            self.logger.info(f"Using protein file from gene prediction: {protein_file}")
            return str(protein_file)
        
        # If no protein file found, check if input file is already protein
        if self._is_protein_file(input_file):
            self.logger.info(f"Input file appears to be protein sequences: {input_file}")
            return input_file
        
        self.logger.warning(f"No protein file found at {protein_file}")
        return None
    
    def _is_protein_file(self, file_path: str) -> bool:
        """Check if a FASTA file contains protein sequences."""
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        continue
                    if line.strip():
                        # Check if sequence contains DNA characters
                        seq = line.strip().upper()
                        dna_chars = set('ATCGN')
                        protein_chars = set('ACDEFGHIKLMNPQRSTVWY*')  # Include stop codons
                        
                        seq_chars = set(seq)
                        if seq_chars.issubset(dna_chars):
                            return False  # DNA sequence
                        elif seq_chars.issubset(protein_chars):
                            return True   # Protein sequence
                        break
        except Exception:
            pass
        return False
    
    def _run_homology_search(self, input_file: str) -> None:
        """Run BLAST or DIAMOND for protein homology search."""
        if getattr(self.config, 'blast_switch', False):
            self._run_blast_search(input_file)
        else:
            self._run_diamond_search(input_file)
    
    def _run_blast_search(self, input_file: str) -> None:
        """Run BLAST homology search."""
        blast_db = getattr(self.config, 'blast_database', None)
        blast_evalue = getattr(self.config, 'blast_evalue', 0.00001)
        ncpus = getattr(self.config, 'ncpus', 1)
        exhaustive = getattr(self.config, 'blast_exh', False)
        
        # Use config's output directory
        output_dir = Path(self.config.root_output) / "protein_function"
        os.makedirs(output_dir, exist_ok=True)
        output_file = output_dir / f"{Path(input_file).stem}.blast.csv"
        
        cmd = [
            'blastp', '-query', input_file, '-db', blast_db,
            '-evalue', str(blast_evalue),
            '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle',
            '-out', str(output_file), '-num_threads', str(ncpus)
        ]
        if exhaustive:
            cmd += ['-word_size', '2', '-gapopen', '8', '-gapextend', '2', '-matrix', 'PAM70', '-comp_based_stats', '0']
        safe_run_cmd(cmd, self.logger)
        self.protein_annotations = self._parse_homology_results(
            str(output_file), 
            getattr(self.config, 'blast_width_threshold', 50.0),
            getattr(self.config, 'blast_cov_threshold', 50.0),
            blast_evalue,
            "blast"
        )
        if output_file.exists():
            output_file.unlink()
    
    def _run_diamond_search(self, input_file: str) -> None:
        """Run DIAMOND homology search."""
        diamond_db = getattr(self.config, 'diamond_database', None)
        diamond_evalue = getattr(self.config, 'diamond_evalue', 0.00001)
        ncpus = getattr(self.config, 'ncpus', 1)
        
        # Use config's output directory
        output_dir = Path(self.config.root_output) / "protein_function"
        os.makedirs(output_dir, exist_ok=True)
        output_file = output_dir / f"{Path(input_file).stem}.diamond.csv"
        
        cmd = [
            'diamond', 'blastp', '-q', input_file, '-d', diamond_db,
            '-e', str(diamond_evalue), '-f', '6', 'qseqid', 'sseqid',
            'pident', 'length', 'qlen', 'slen', 'qstart', 'qend',
            'evalue', 'bitscore', 'stitle', '-o', str(output_file),
            '-p', str(ncpus), '--quiet'
        ]
        safe_run_cmd(cmd, self.logger)
        self.protein_annotations = self._parse_homology_results(
            str(output_file),
            getattr(self.config, 'diamond_width_threshold', 50.0),
            getattr(self.config, 'diamond_cov_threshold', 50.0),
            diamond_evalue,
            "diamond"
        )
        # Keep the DIAMOND results file for user access
        # if output_file.exists():
        #     output_file.unlink()
    
    def _parse_homology_results(self, file_path: str, threshold_width: float, threshold_cov: float, threshold_evalue: float, program: str) -> Dict[str, Any]:
        """
        Parse BLAST or DIAMOND results for protein annotation.
        """
        hypotheticalpat = re.compile(r'(?i)(((hypothetical|uncharacteri[z|s]ed|predicted)( phage)?( membrane)? protein)|(ORF|(unnamed protein product|gp\d+|protein of unknown function|phage protein)))')
        annotations = {}
        with open(file_path, "r") as results:
            reader = csv.DictReader(results, delimiter='\t', 
                                  fieldnames=['qseqid','sseqid','pident','length','qlen','slen','qstart','qend','evalue','bitscore','stitle'])
            for row in reader:
                perc_cover = round(100.00 * (float(row['length']) / float(row['qlen'])), 2)
                perc_id = float(row['pident'])
                ann = {
                    'sseqid': row['sseqid'],
                    'pident': perc_id,
                    'pcover': perc_cover,
                    'evalue': row['evalue'],
                    'descr': row['stitle']
                }
                if (not re.search(hypotheticalpat, ann['descr']) and perc_id >= threshold_width and perc_cover >= threshold_cov and float(row['evalue']) <= threshold_evalue):
                    if row['qseqid'] not in annotations or perc_id > annotations[row['qseqid']]['pident']:
                        annotations[row['qseqid']] = ann
        return annotations

    def _process_proteins(self, input_file: str) -> None:
        """Process protein sequences and create annotation dictionary."""
        records = list(SeqIO.parse(open(input_file, "r"), "fasta"))
        self.protein_results = self._init_protein_results(records)
        self.protein_results = self._update_protein_results(records, self.protein_annotations)
    
    def _init_protein_results(self, records) -> Dict[str, Any]:
        """Initialize protein results structure."""
        results = {}
        for record in records:
            contig_id = record.description.split(' # ')[0].split('_')[0]
            if contig_id not in results:
                results[contig_id] = {}
        return results
    
    def _update_protein_results(self, records, annotations: Dict[str, Any]) -> Dict[str, Any]:
        """Update protein results with annotation information."""
        for record in records:
            contig_id, prot_id, protinfo = self._extract_protein_features(record, annotations)
            self.protein_results[contig_id][prot_id] = protinfo
        return self.protein_results
    
    def _extract_protein_features(self, record, annotations: Dict[str, Any]) -> tuple:
        """Extract protein features and annotation from record."""
        dataprot = record.description.split(' # ')
        contig_id = dataprot[0].split('_')[0]
        prot_id = dataprot[0]
        
        # Remove stop codons (*) for protein analysis
        modseq = str(record.seq).replace("X", "").replace("*", "")
        
        # Only analyze if we have a valid protein sequence
        if len(modseq) > 0:
            try:
                analysed_seq = ProteinAnalysis(modseq)
                protinfo = {
                    'length': len(record.seq),
                    'isoelectric_point': analysed_seq.isoelectric_point(),
                    'molecular_weight_kda': analysed_seq.molecular_weight()/1000.00,
                    'instability_index': analysed_seq.instability_index(),
                    'protein_id': prot_id,
                    'translation': record.seq
                }
            except Exception as e:
                # Fallback for problematic sequences
                self.logger.warning(f"Protein analysis failed for {prot_id}: {e}. Using basic info only.")
                protinfo = {
                    'length': len(record.seq),
                    'isoelectric_point': 0.0,
                    'molecular_weight_kda': 0.0,
                    'instability_index': 0.0,
                    'protein_id': prot_id,
                    'translation': record.seq
                }
        else:
            # Handle empty sequences
            protinfo = {
                'length': len(record.seq),
                'isoelectric_point': 0.0,
                'molecular_weight_kda': 0.0,
                'instability_index': 0.0,
                'protein_id': prot_id,
                'translation': record.seq
            }
        
        hypotheticalpat = re.compile(r'(?i)(((hypothetical|uncharacteri(z|s)ed|predicted))( phage)?( membrane)? protein|(ORF|(unnamed protein product|gp\d+|protein of unknown function|phage protein)))')
        ann = annotations.get(prot_id)
        
        # Add HMMER domain information
        hmmer_ann = getattr(self, 'hmmer_annotations', {}).get(prot_id)
        if hmmer_ann:
            protinfo.update({
                'hmmer_domain': hmmer_ann['target_name'],
                'hmmer_accession': hmmer_ann['target_accession'],
                'hmmer_database': hmmer_ann['database'],
                'hmmer_score': hmmer_ann['domain_score'],
                'hmmer_evalue': hmmer_ann['domain_evalue'],
                'hmmer_coverage': hmmer_ann['coverage']
            })
        else:
            protinfo.update({
                'hmmer_domain': 'No domain',
                'hmmer_accession': 'NA',
                'hmmer_database': 'NA',
                'hmmer_score': 'NA',
                'hmmer_evalue': 'NA',
                'hmmer_coverage': 'NA'
            })
        
        if ann is None:
            protinfo.update({
                'descr': 'Hypothetical protein',
                'source': "NO_HIT",
                'pident': "NA",
                'pcover': "NA",
                'evalue': "NA"
            })
        elif re.search(hypotheticalpat, ann['descr']):
            protinfo.update({
                'descr': 'Conserved hypothetical protein',
                'source': "Homology",
                'pident': ann['pident'],
                'pcover': ann['pcover'],
                'evalue': ann['evalue']
            })
        else:
            descr = " ".join(ann['descr'].split()[1:])
            protinfo.update({
                'descr': descr,
                'source': "Homology",
                'pident': ann['pident'],
                'pcover': ann['pcover'],
                'evalue': ann['evalue']
            })
        return contig_id, prot_id, protinfo
    
    def _run_hmmer_annotation(self, input_file: str) -> None:
        """Run HMMER for protein domain annotation using viral databases."""
        # Check if HMMER is available
        if not cmd_exists("hmmscan"):
            self.logger.warning("HMMER not found. Skipping HMMER domain annotation.")
            self.logger.info("To install HMMER: conda install -c bioconda hmmer")
            return
        
        # Get viral database paths from config
        viral_databases = {}
        if not getattr(self.config, 'no_vfam', False):
            vfam_db = getattr(self.config, 'vfam_database', None)
            if vfam_db and Path(vfam_db).exists():
                viral_databases['VFAM'] = vfam_db
            else:
                self.logger.warning("VFAM database not found. Skipping VFAM annotation.")
        
        if not getattr(self.config, 'no_rvdb', False):
            rvdb_db = getattr(self.config, 'rvdb_database', None)
            if rvdb_db and Path(rvdb_db).exists():
                viral_databases['RVDB'] = rvdb_db
            else:
                self.logger.warning("RVDB database not found. Skipping RVDB annotation.")
        
        if not getattr(self.config, 'no_vogs', False):
            vogs_db = getattr(self.config, 'vogs_database', None)
            if vogs_db and Path(vogs_db).exists():
                viral_databases['VOGs'] = vogs_db
            else:
                self.logger.warning("VOGs database not found. Skipping VOGs annotation.")
        
        if not getattr(self.config, 'no_phrogs', False):
            phrogs_db = getattr(self.config, 'phrogs_database', None)
            if phrogs_db and Path(phrogs_db).exists():
                viral_databases['PHROGS'] = phrogs_db
            else:
                self.logger.warning("PHROGS database not found. Skipping PHROGS annotation.")
        
        if not viral_databases:
            self.logger.warning("No viral databases found. Skipping HMMER domain annotation.")
            self.logger.info("To set up viral databases, run: viranpy --build-databases")
            return
        
        hmmer_evalue = getattr(self.config, 'hmmer_evalue', 0.001)
        hmmer_coverage = getattr(self.config, 'hmmer_cov_threshold', 50.0) / 100.0  # Convert percentage to decimal
        ncpus = getattr(self.config, 'ncpus', 1)
        
        # Use config's output directory
        output_dir = Path(self.config.root_output) / "protein_function"
        os.makedirs(output_dir, exist_ok=True)
        
        all_annotations = {}
        
        # Run HMMER against each viral database
        for db_name, db_path in viral_databases.items():
            self.logger.info(f"Running HMMER against {db_name} database...")
            output_file = output_dir / f"{Path(input_file).stem}.{db_name.lower()}.hmmer.txt"
            
            cmd = [
                'hmmscan', '--cpu', str(ncpus), '--domtblout', str(output_file),
                '-E', str(hmmer_evalue), '--domE', str(hmmer_evalue),
                db_path, input_file
            ]
            
            try:
                safe_run_cmd(cmd, self.logger)
                db_annotations = self._parse_hmmer_results(str(output_file), hmmer_coverage, db_name)
                all_annotations.update(db_annotations)
                self.logger.info(f"{db_name} annotation completed: {len(db_annotations)} domain annotations found")
            except Exception as e:
                self.logger.error(f"{db_name} annotation failed: {e}")
        
        self.hmmer_annotations = all_annotations
        self.logger.info(f"Total HMMER annotation completed: {len(self.hmmer_annotations)} domain annotations found across all databases")
    
    def _parse_hmmer_results(self, file_path: str, coverage_threshold: float, database_name: str) -> Dict[str, Any]:
        """Parse HMMER domain table output."""
        annotations = {}
        
        if not Path(file_path).exists():
            return annotations
        
        try:
            with open(file_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split()
                    if len(fields) < 22:
                        continue
                    
                    try:
                        # Parse HMMER domtblout format with proper error handling
                        query_name = fields[0]
                        target_name = fields[3]
                        target_accession = fields[4]
                        
                        # Convert numeric fields with proper error handling
                        query_len = int(float(fields[5]))  # Handle scientific notation
                        target_len = int(float(fields[6]))  # Handle scientific notation
                        full_evalue = float(fields[7])
                        full_score = float(fields[8])
                        full_bias = float(fields[9])
                        dom_num = int(float(fields[10]))  # Handle scientific notation
                        dom_total = int(float(fields[11]))  # Handle scientific notation
                        dom_cvalue = float(fields[12])
                        dom_ivalue = float(fields[13])
                        dom_score = float(fields[14])
                        dom_bias = float(fields[15])
                        hmm_start = int(float(fields[16]))  # Handle scientific notation
                        hmm_end = int(float(fields[17]))  # Handle scientific notation
                        ali_start = int(float(fields[18]))  # Handle scientific notation
                        ali_end = int(float(fields[19]))  # Handle scientific notation
                        env_start = int(float(fields[20]))  # Handle scientific notation
                        env_end = int(float(fields[21]))  # Handle scientific notation
                        
                        # Calculate coverage
                        domain_length = ali_end - ali_start + 1
                        coverage = domain_length / query_len
                        
                        # Apply coverage threshold
                        if coverage >= coverage_threshold:
                            domain_info = {
                                'target_name': target_name,
                                'target_accession': target_accession,
                                'database': database_name,
                                'domain_score': dom_score,
                                'domain_evalue': dom_ivalue,
                                'coverage': coverage,
                                'ali_start': ali_start,
                                'ali_end': ali_end,
                                'hmm_start': hmm_start,
                                'hmm_end': hmm_end
                            }
                            
                            # Keep best hit per protein (highest score)
                            if query_name not in annotations or dom_score > annotations[query_name]['domain_score']:
                                annotations[query_name] = domain_info
                    
                    except (ValueError, IndexError) as e:
                        # Log the problematic line but continue parsing
                        self.logger.warning(f"Skipping malformed line {line_num} in {database_name} results: {line.strip()[:100]}... Error: {e}")
                        continue
                            
        except Exception as e:
            self.logger.error(f"Error parsing HMMER results from {database_name}: {e}")
        
        return annotations
    
    def get_output_files(self) -> List[str]:
        return [f"protein_results_{contig_id}.csv" for contig_id in self.protein_results.keys()] 