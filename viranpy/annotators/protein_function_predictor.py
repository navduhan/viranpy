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
            self._run_homology_search(input_file)
            self._process_proteins(input_file)
            if not getattr(self.config, 'no_hmmer', False):
                self._run_hmmer_annotation(input_file)
            result.add_annotation("protein_results", self.protein_results)
            result.add_annotation("homology_annotations", self.protein_annotations)
            result.add_metadata("total_proteins", sum(len(v) for v in self.protein_results.values()))
            self.logger.info(f"Protein function annotation completed: {result.metadata['total_proteins']} proteins annotated")
        except Exception as e:
            result.success = False
            result.error_message = str(e)
            self.logger.error(f"Protein function annotation failed: {e}")
        return result.to_dict()
    
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
        output_file = f"{input_file}.blast.csv"
                    cmd = [
                        'blastp', '-query', input_file, '-db', blast_db,
                        '-evalue', str(blast_evalue),
                        '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend evalue bitscore stitle',
            '-out', output_file, '-num_threads', str(ncpus)
        ]
        if exhaustive:
            cmd += ['-word_size', '2', '-gapopen', '8', '-gapextend', '2', '-matrix', 'PAM70', '-comp_based_stats', '0']
        safe_run_cmd(cmd, self.logger)
        self.protein_annotations = self._parse_homology_results(
            output_file, 
            getattr(self.config, 'blast_width_threshold', 50.0),
            getattr(self.config, 'blast_cov_threshold', 50.0),
            blast_evalue,
            "blast"
        )
        if os.path.exists(output_file):
            os.remove(output_file)
    
    def _run_diamond_search(self, input_file: str) -> None:
        """Run DIAMOND homology search."""
        diamond_db = getattr(self.config, 'diamond_database', None)
        diamond_evalue = getattr(self.config, 'diamond_evalue', 0.00001)
        ncpus = getattr(self.config, 'ncpus', 1)
        output_file = f"{input_file}.diamond.csv"
                cmd = [
                    'diamond', 'blastp', '-q', input_file, '-d', diamond_db,
                    '-e', str(diamond_evalue), '-f', '6', 'qseqid', 'sseqid',
                    'pident', 'length', 'qlen', 'slen', 'qstart', 'qend',
                    'evalue', 'bitscore', 'stitle', '-o', output_file,
            '-p', str(ncpus), '--quiet'
                ]
        safe_run_cmd(cmd, self.logger)
        self.protein_annotations = self._parse_homology_results(
            output_file,
            getattr(self.config, 'diamond_width_threshold', 50.0),
            getattr(self.config, 'diamond_cov_threshold', 50.0),
            diamond_evalue,
            "diamond"
        )
        if os.path.exists(output_file):
            os.remove(output_file)
    
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
        modseq = str(record.seq).replace("X", "")
        analysed_seq = ProteinAnalysis(modseq)
        protinfo = {
            'length': len(record.seq),
            'isoelectric_point': analysed_seq.isoelectric_point(),
            'molecular_weight_kda': analysed_seq.molecular_weight()/1000.00,
            'instability_index': analysed_seq.instability_index(),
            'protein_id': prot_id,
            'translation': record.seq
        }
        hypotheticalpat = re.compile(r'(?i)(((hypothetical|uncharacteri(z|s)ed|predicted))( phage)?( membrane)? protein|(ORF|(unnamed protein product|gp\d+|protein of unknown function|phage protein)))')
        ann = annotations.get(prot_id)
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
        """Run HMMER for protein function refinement (placeholder)."""
        self.logger.info("HMMER annotation would be implemented here.")
    
    def get_output_files(self) -> List[str]:
        return [f"protein_results_{contig_id}.csv" for contig_id in self.protein_results.keys()] 