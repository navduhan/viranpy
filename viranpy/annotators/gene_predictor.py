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
            # Use config's output directory for output files
            output_dir = Path(self.config.root_output) / "gene_prediction"
            os.makedirs(output_dir, exist_ok=True)
            
            # Combined protein file for downstream analysis
            combined_protein_file = output_dir / f"{Path(input_file).stem}_proteins.faa"
            
            for record in SeqIO.parse(input_file, "fasta"):
                self._find_genes_in_sequence(record, input_file, combined_protein_file)
            
            result.add_annotation("predicted_proteins", self.predicted_proteins)
            result.add_metadata("total_genes", sum(len(v) for v in self.predicted_proteins.values()))
            result.add_metadata("protein_file", str(combined_protein_file))
            self.logger.info(f"Gene prediction completed: {result.metadata['total_genes']} genes found")
        except Exception as e:
            result.success = False
            result.error_message = str(e)
            self.logger.error(f"Gene prediction failed: {e}")
        return result.to_dict()

    def _find_genes_in_sequence(self, record, fasta_path: str, combined_protein_file: Path) -> None:
        """
        Predict genes for a single viral sequence using Prodigal/Prodigal-GV.
        Results are stored in self.predicted_proteins.
        """
        # Validate sequence before processing
        if not self._is_valid_sequence(record):
            self.logger.warning(f"Contig {record.id} is too short or invalid for gene prediction. Skipping.")
            self.predicted_proteins[record.id] = []
            return
        
        # Use config's output directory for temporary and output files
        output_dir = Path(self.config.root_output) / "gene_prediction"
        os.makedirs(output_dir, exist_ok=True)
        
        temp_fasta = output_dir / f"temp_{record.id}.fasta"
        protein_faa = output_dir / f"proteins_{record.id}.faa"
        nucleotide_fna = output_dir / f"proteins_{record.id}.fna"
        
        with open(temp_fasta, "w") as f:
            SeqIO.write(record, f, "fasta")
        
        use_gv = getattr(self.config, 'use_prodigal_gv', False)
        genetic_code = str(getattr(self.config, 'genetic_code_table', 11))
        genome_shape = getattr(self, 'topology_results', {}).get(record.id, {}).get("topology", "linear")
        
        if use_gv:
            cmd = ["prodigal-gv", "-p", "meta", "-i", str(temp_fasta), "-a", str(protein_faa), "-d", str(nucleotide_fna), "-o", "/dev/null", "-q"]
            if genome_shape == 'linear':
                cmd += ["-c"]
        else:
            cmd = ["prodigal", "-a", str(protein_faa), "-d", str(nucleotide_fna), "-i", str(temp_fasta), "-o", "/dev/null", "-g", genetic_code, "-q"]
            if len(record.seq) >= 100000:
                if genome_shape == 'linear':
                    cmd += ["-c"]
            else:
                cmd += ["-p", "meta"]
                if genome_shape == 'linear':
                    cmd += ["-c"]
        
        try:
            safe_run_cmd(cmd, self.logger)
            self.predicted_proteins[record.id] = self._parse_predicted_proteins(str(protein_faa))
            
            # Append proteins to combined file
            if protein_faa.exists():
                with open(protein_faa, 'r') as infile, open(combined_protein_file, 'a') as outfile:
                    outfile.write(infile.read())
            
        except Exception as e:
            # Handle segmentation faults and other Prodigal crashes
            self.logger.warning(f"Prodigal failed for contig {record.id}: {e}. Skipping this contig.")
            self.predicted_proteins[record.id] = []  # Empty list for failed contigs
        
        # Clean up temporary files but keep individual protein files for now
        if temp_fasta.exists():
            temp_fasta.unlink()
        if nucleotide_fna.exists():
            nucleotide_fna.unlink()

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

    def _is_valid_sequence(self, record) -> bool:
        """Check if a sequence is valid for gene prediction."""
        # Check sequence length (Prodigal needs at least 200 bp)
        if len(record.seq) < 200:
            return False
        
        # Check for valid DNA characters
        valid_chars = set('ATCGN')
        seq_chars = set(str(record.seq).upper())
        if not seq_chars.issubset(valid_chars):
            return False
        
        # Check for too many N's (more than 50% N's)
        n_count = str(record.seq).upper().count('N')
        if n_count / len(record.seq) > 0.5:
            return False
        
        return True 