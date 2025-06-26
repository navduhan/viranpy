#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan
"""
ncRNA detection in viral genomes using RFAM and cmscan.
"""

import os
import subprocess
import re
from typing import Dict, Any, List
from pathlib import Path
from ..core.base import BaseAnnotator
from ..core.results import AnnotationResult
from ..utils.file_utils import cmd_exists, safe_run_cmd

class ViralNcRNAFinder(BaseAnnotator):
    """
    Detect non-coding RNAs (ncRNAs) in viral genomes using RFAM and cmscan.
    Unique ViRAnPy implementation.
    """
    def __init__(self, config, logger=None):
        super().__init__(config, logger)
        self.ncrna_results = {}
    
    def check_dependencies(self) -> bool:
        """Check if cmscan is available."""
        return cmd_exists("cmscan")
    
    def validate_input(self, input_file: str) -> bool:
        """Validate input FASTA file."""
        return Path(input_file).exists() and Path(input_file).suffix.lower() in ['.fasta', '.fa', '.fna']
    
    def run(self, input_file: str, **kwargs) -> Dict[str, Any]:
        """
        Run ncRNA detection for all sequences in the input file.
        Returns a dictionary with annotation results.
        """
        result = AnnotationResult(
            annotator_name=self.name,
            input_file=input_file
        )
        try:
            rfam_db = getattr(self.config, 'rfam_database', None)
            ncpus = getattr(self.config, 'ncpus', 1)
            norfam = getattr(self.config, 'nor_fam', False)
            hmmonly = getattr(self.config, 'hmmonly', False)
            
            # Use config's output directory
            output_dir = Path(self.config.root_output) / "ncrna_detection"
            os.makedirs(output_dir, exist_ok=True)
            output_file = output_dir / f"ncrna_results_{Path(input_file).stem}.csv"
            
            cmd = ["cmscan"]
            if not norfam:
                cmd.append("--rfam")
            elif hmmonly:
                cmd.append("--hmmonly")
            cmd.extend(["--cut_ga", "--tblout", str(output_file), "--cpu", str(ncpus), rfam_db, input_file])
            safe_run_cmd(cmd, self.logger)
            self.ncrna_results = self._parse_ncrna_output(str(output_file))
            
            # Write results to log file
            self._write_ncrna_log()
            
            result.add_annotation("ncrna_results", self.ncrna_results)
            result.add_metadata("total_ncrna", sum(len(v) for v in self.ncrna_results.values()))
            self.logger.info(f"ncRNA detection completed: {result.metadata['total_ncrna']} ncRNAs found")
            
            # Clean up output file
            if output_file.exists():
                output_file.unlink()
        except Exception as e:
            result.success = False
            result.error_message = str(e)
            self.logger.error(f"ncRNA detection failed: {e}")
        return result.to_dict()

    def _parse_ncrna_output(self, filename: str) -> Dict[str, List[Dict[str, Any]]]:
        """
        Parse cmscan output to extract ncRNA information.
        Returns a dictionary keyed by contig ID.
        """
        ncrna_dict = {}
        with open(filename, "r") as ncrnafile:
            for line in ncrnafile:
                if line.startswith("#"):
                    continue
                fields = re.sub("\s{2,}", ",", line).split(",")
                if len(fields) < 15:
                        continue
                contig_id = fields[1] if len(fields) == 15 else fields[2]
                if contig_id not in ncrna_dict:
                    ncrna_dict[contig_id] = []
                ncrna_info = {
                    'type': fields[0],
                    'product': fields[-1].strip(),
                    'begin': int(fields[6]) if len(fields) == 15 else int(fields[7]),
                    'end': int(fields[7]) if len(fields) == 15 else int(fields[8]),
                    'score': float(fields[13].replace(" !", "")) if len(fields) == 15 else float(fields[14].replace(" !", "")),
                    'rfam_code': fields[0].split()[-1] if len(fields) == 15 else fields[1],
                    'strand': 1 if fields[8] == "+" else -1 if len(fields) == 15 else 1 if fields[9] == "+" else -1,
                    'locus_tag': f"{contig_id}_nc{len(ncrna_dict[contig_id])+1}"
                }
                ncrna_dict[contig_id].append(ncrna_info)
        return ncrna_dict

    def get_output_files(self) -> List[str]:
        return [f"ncrna_results_{contig_id}.csv" for contig_id in self.ncrna_results.keys()]

    def _write_ncrna_log(self) -> None:
        """Write ncRNA results to a log file."""
        output_dir = Path(self.config.root_output) / "ncrna_detection"
        log_file = output_dir / "ncrna_results.txt"
        
        with open(log_file, 'w') as f:
            f.write("# ncRNA Detection Results\n")
            f.write("# Contig_ID\tLocus_Tag\tType\tProduct\tStart\tEnd\tStrand\tScore\tRFAM_Code\n")
            
            for contig_id, ncrnas in self.ncrna_results.items():
                for ncrna in ncrnas:
                    f.write(f"{contig_id}\t{ncrna['locus_tag']}\t{ncrna['type']}\t{ncrna['product']}\t{ncrna['begin']}\t{ncrna['end']}\t{ncrna['strand']}\t{ncrna['score']}\t{ncrna['rfam_code']}\n")
        
        self.logger.info(f"ncRNA results written to: {log_file}") 