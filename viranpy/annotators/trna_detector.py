#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan
"""
tRNA and tmRNA detection using ARAGORN.
"""

import os
import subprocess
import re
from typing import Dict, Any, List
from pathlib import Path
from Bio import SeqIO
from ..core.base import BaseAnnotator
from ..core.results import AnnotationResult
from ..utils.file_utils import safe_run_cmd
from ..utils.sequence_utils import get_sequence_length_from_fasta

class ViralTRNAFinder(BaseAnnotator):
    """
    Identify tRNA and tmRNA genes in viral genomes using ARAGORN.
    Unique ViRAnPy implementation.
    """
    def __init__(self, config, logger=None):
        super().__init__(config, logger)
        self.trna_results = {}
        self.tmrna_results = {}

    def check_dependencies(self) -> bool:
        """Check if ARAGORN is available."""
        from ..utils.file_utils import cmd_exists
        return cmd_exists("aragorn")

    def validate_input(self, input_file: str) -> bool:
        """Validate input FASTA file."""
        return Path(input_file).exists() and Path(input_file).suffix.lower() in ['.fasta', '.fa', '.fna']

    def run(self, input_file: str, **kwargs) -> Dict[str, Any]:
        """
        Run tRNA/tmRNA detection for all sequences in the input file.
        Returns a dictionary with annotation results.
        """
        result = AnnotationResult(
            annotator_name=self.name,
            input_file=input_file
        )
        try:
            for record in SeqIO.parse(input_file, "fasta"):
                self._find_trna_in_sequence(record)
            result.add_annotation("trna_results", self.trna_results)
            result.add_annotation("tmrna_results", self.tmrna_results)
            result.add_metadata("total_trna", sum(len(v) for v in self.trna_results.values()))
            result.add_metadata("total_tmrna", sum(len(v) for v in self.tmrna_results.values()))
            self.logger.info(f"tRNA/tmRNA detection completed: {result.metadata['total_trna']} tRNAs, {result.metadata['total_tmrna']} tmRNAs found")
        except Exception as e:
            result.success = False
            result.error_message = str(e)
            self.logger.error(f"tRNA/tmRNA detection failed: {e}")
        return result.to_dict()

    def _find_trna_in_sequence(self, record) -> None:
        """
        Detect tRNA/tmRNA for a single viral sequence using ARAGORN.
        Results are stored in self.trna_results and self.tmrna_results.
        """
        temp_fasta = f"temp_{record.id}.fasta"
        with open(temp_fasta, "w") as f:
            SeqIO.write(record, f, "fasta")
        output_file = f"trna_results_{record.id}.fasta"
        cmd = [
            "aragorn", "-fon", f"-gc{self.config.genetic_code_table}", temp_fasta
        ]
        genome_shape = getattr(self, 'topology_results', {}).get(record.id, {}).get("topology", "linear")
        if genome_shape == "circular":
            cmd.append("-c")
        else:
            cmd.append("-l")
        safe_run_cmd(cmd, self.logger, output_file=output_file)
        self._parse_aragorn_output(record.id, output_file, temp_fasta)
        os.remove(temp_fasta)
        if os.path.exists(output_file):
            os.remove(output_file)

    def _parse_aragorn_output(self, contig_id: str, output_file: str, fasta_file: str) -> None:
        """
        Parse ARAGORN output to extract tRNA/tmRNA information.
        Results are stored in self.trna_results and self.tmrna_results.
        """
        self.trna_results[contig_id] = []
        self.tmrna_results[contig_id] = []
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            return
        for tRNA_seq in SeqIO.parse(output_file, "fasta"):
            tRNA_information = tRNA_seq.description.split(" ")
            if tRNA_information[1] == "tmRNA":
                tmrna_info = self._extract_tmrna_info(tRNA_information, fasta_file)
                self.tmrna_results[contig_id].append(tmrna_info)
            elif re.match("^tRNA-", tRNA_information[1]):
                trna_info = self._extract_trna_info(tRNA_information, fasta_file)
                self.trna_results[contig_id].append(trna_info)

    def _extract_trna_info(self, tRNA_data: List[str], fasta_file: str) -> Dict[str, Any]:
        """
        Extract tRNA information from ARAGORN output.
        """
        contig_id = Path(fasta_file).stem.replace("temp_", "")
        seq_length = get_sequence_length_from_fasta(contig_id, fasta_file)
        info = {}
        info['product'] = re.sub(r"\(\w{3}\)", "", tRNA_data[1])
        coords = tRNA_data[2] if tRNA_data[2] != "(Permuted)" else tRNA_data[3]
        info['strand'] = -1 if coords.startswith("c") else 1
        coords = coords.replace("c[", "").replace("[", "").replace("]", "").split(",")
        begin = int(coords[0])
        end = int(coords[1])
        info['begin'] = max(begin, 1)
        info['end'] = min(end, seq_length)
        if info['begin'] > info['end']:
            info['begin'], info['end'] = info['end'], info['begin']
        info['locus_tag'] = f"{contig_id}_t{tRNA_data[0]}"
        return info

    def _extract_tmrna_info(self, tRNA_data: List[str], fasta_file: str) -> Dict[str, Any]:
        """
        Extract tmRNA information from ARAGORN output.
        """
        contig_id = Path(fasta_file).stem.replace("temp_", "")
        seq_length = get_sequence_length_from_fasta(contig_id, fasta_file)
        info = {}
        info['product'] = tRNA_data[1]
        coords = tRNA_data[2] if tRNA_data[2] != "(Permuted)" else tRNA_data[3]
        info['strand'] = -1 if coords.startswith("c") else 1
        coords = coords.replace("c[", "").replace("[", "").replace("]", "").split(",")
        begin = int(coords[0])
        end = int(coords[1])
        info['begin'] = max(begin, 1)
        info['end'] = min(end, seq_length)
        if info['begin'] > info['end']:
            info['begin'], info['end'] = info['end'], info['begin']
        info['locus_tag'] = f"{contig_id}_tm{tRNA_data[0]}"
        return info

    def get_output_files(self) -> List[str]:
        return [f"trna_results_{contig_id}.fasta" for contig_id in self.trna_results.keys()] 