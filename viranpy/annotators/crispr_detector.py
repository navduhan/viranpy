#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan
"""
CRISPR array detection in viral genomes using PILER-CR.
"""

import os
import re
import subprocess
from typing import Dict, Any, List
from pathlib import Path
from ..core.base import BaseAnnotator
from ..core.results import AnnotationResult
from ..utils.file_utils import cmd_exists, safe_run_cmd

class ViralCRISPRFinder(BaseAnnotator):
    """
    Detect CRISPR arrays in viral genomes using PILER-CR.
    Unique ViRAnPy implementation.
    """
    def __init__(self, config, logger=None):
        super().__init__(config, logger)
        self.crispr_results = {}
    
    def check_dependencies(self) -> bool:
        """Check if PILER-CR is available."""
        return cmd_exists("pilercr")
    
    def validate_input(self, input_file: str) -> bool:
        """Validate input FASTA file."""
        return Path(input_file).exists() and Path(input_file).suffix.lower() in ['.fasta', '.fa', '.fna']
    
    def run(self, input_file: str, **kwargs) -> Dict[str, Any]:
        """
        Run CRISPR array detection for all sequences in the input file.
        Returns a dictionary with annotation results.
        """
        result = AnnotationResult(
            annotator_name=self.name,
            input_file=input_file
        )
        try:
            min_repeat = getattr(self.config, 'min_crispr_repeat', 16)
            max_repeat = getattr(self.config, 'max_crispr_repeat', 64)
            min_spacer = getattr(self.config, 'min_crispr_spacer', 8)
            max_spacer = getattr(self.config, 'max_crispr_spacer', 64)
            output_file = f"crispr_results_{Path(input_file).stem}.txt"
            cmd = [
                "pilercr", "-in", input_file, "-out", output_file, "-noinfo",
                "-minrepeat", str(min_repeat), "-maxrepeat", str(max_repeat),
                "-minspacer", str(min_spacer), "-maxspacer", str(max_spacer)
            ]
            safe_run_cmd(cmd, self.logger)
            self.crispr_results = self._parse_crispr_output(output_file)
            result.add_annotation("crispr_results", self.crispr_results)
            result.add_metadata("total_crispr_arrays", sum(len(v) for v in self.crispr_results.values()))
            self.logger.info(f"CRISPR detection completed: {result.metadata['total_crispr_arrays']} arrays found")
            if os.path.exists(output_file):
                os.remove(output_file)
        except Exception as e:
            result.success = False
            result.error_message = str(e)
            self.logger.error(f"CRISPR detection failed: {e}")
        return result.to_dict()

    def _parse_crispr_output(self, filename: str) -> Dict[str, Any]:
        """
        Parse PILER-CR output to extract CRISPR array information.
        Returns a dictionary keyed by sequence ID.
        """
        crispr_arrays = {}
        with open(filename, "r") as crisprfile:
            in_summary = False
            for line in crisprfile:
                if "SUMMARY BY POSITION" in line:
                    in_summary = True
                            continue
                if in_summary:
                    match = re.match(r"^\s+(\d+)\s+(.{16})\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+\d?\s+(\w+)", line)
                    if match:
                        key, seq_id, start, length, repeat_seq = match.groups()
                        seq_id = seq_id.strip()
                        if seq_id not in crispr_arrays:
                            crispr_arrays[seq_id] = []
                        crispr_arrays[seq_id].append({
                            'start': int(start),
                            'end': int(start) + int(length),
                            'repeat_sequence': repeat_seq,
                            'repeat_end': int(start) + len(repeat_seq)
                        })
        return crispr_arrays

    def get_output_files(self) -> List[str]:
        return [f"crispr_results_{contig_id}.txt" for contig_id in self.crispr_results.keys()] 