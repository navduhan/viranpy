#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan
"""
Viral gene prediction using Prodigal.
"""

import os
import subprocess
from typing import Dict, Any
from pathlib import Path
from Bio import SeqIO

from ..core.base import BaseAnnotator
from ..core.results import AnnotationResult
from ..utils.file_utils import safe_run_cmd

class ViralGeneFinder(BaseAnnotator):
    """
    Predict protein-coding genes in viral genomes using Prodigal or Prodigal-GV.
    This ViRAnPy implementation is modular and distinct from VIGA.
    """
    
    def __init__(self, config, logger=None):
        super().__init__(config, logger)
        self.predicted_proteins = {}
    
    def check_dependencies(self) -> bool:
        """Check if Prodigal or Prodigal-GV is available."""
        from ..utils.file_utils import cmd_exists
        return cmd_exists("prodigal-gv") or cmd_exists("prodigal")
    
    def validate_input(self, input_file: str) -> bool:
        """Validate input FASTA file."""
        return Path(input_file).exists() and Path(input_file).suffix.lower() in ['.fasta', '.fa', '.fna']
    
    def run(self, input_file: str, **kwargs) -> Dict[str, Any]:
        """
        Run viral gene prediction for all sequences in the input file.
        Returns a dictionary with annotation results.
        """
        result = AnnotationResult(
            annotator_name=self.name,
            input_file=input_file
        )
        try:
            for record in SeqIO.parse(input_file, "fasta"):
                self._find_genes_in_sequence(record, input_file)
            result.add_annotation("predicted_proteins", self.predicted_proteins)
            result.add_metadata("total_genes", sum(len(v) for v in self.predicted_proteins.values()))
            self.logger.info(f"Gene prediction completed: {result.metadata['total_genes']} genes found")
        except Exception as e:
            result.success = False
            result.error_message = str(e)
            self.logger.error(f"Gene prediction failed: {e}")
        return result.to_dict()

    def _find_genes_in_sequence(self, record, fasta_path: str) -> None:
        """
        Predict genes for a single viral sequence using Prodigal/Prodigal-GV.
        Results are stored in self.predicted_proteins.
        """
        temp_fasta = f"temp_{record.id}.fasta"
        with open(temp_fasta, "w") as f:
            SeqIO.write(record, f, "fasta")
        protein_faa = f"proteins_{record.id}.faa"
        nucleotide_fna = f"proteins_{record.id}.fna"
        use_gv = getattr(self.config, 'use_prodigal_gv', False)
        genetic_code = str(getattr(self.config, 'genetic_code_table', 11))
        genome_shape = getattr(self, 'topology_results', {}).get(record.id, {}).get("topology", "linear")
        if use_gv:
            cmd = ["prodigal-gv", "-p", "meta", "-i", temp_fasta, "-a", protein_faa, "-d", nucleotide_fna, "-o", "/dev/null", "-q"]
            if genome_shape == 'linear':
                cmd += ["-c"]
        else:
            cmd = ["prodigal", "-a", protein_faa, "-d", nucleotide_fna, "-i", temp_fasta, "-o", "/dev/null", "-g", genetic_code, "-q"]
            if len(record.seq) >= 100000:
                if genome_shape == 'linear':
                    cmd += ["-c"]
            else:
                cmd += ["-p", "meta"]
                if genome_shape == 'linear':
                    cmd += ["-c"]
        safe_run_cmd(cmd, self.logger)
        self.predicted_proteins[record.id] = self._parse_predicted_proteins(protein_faa)
        os.remove(temp_fasta)
        if os.path.exists(protein_faa):
            os.remove(protein_faa)
        if os.path.exists(nucleotide_fna):
            os.remove(nucleotide_fna)

    def _parse_predicted_proteins(self, faa_file: str) -> list:
        """
        Parse predicted protein sequences from a FASTA file.
        Returns a list of protein records (dicts).
        """
        proteins = []
        if not os.path.exists(faa_file):
            return proteins
        for seq_record in SeqIO.parse(faa_file, "fasta"):
            proteins.append({
                "id": seq_record.id,
                "description": seq_record.description,
                "sequence": str(seq_record.seq).rstrip("*")
            })
        return proteins

    def get_output_files(self) -> list:
        return [f"proteins_{contig_id}.faa" for contig_id in self.predicted_proteins.keys()] 