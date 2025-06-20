#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan
"""
Genome topology prediction (circular vs linear) using LASTZ.
"""

import os
import subprocess
import fractions
from typing import Dict, Any, List, Tuple
from pathlib import Path
import numpy as np
import scipy.signal as signal
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ..core.base import BaseAnnotator
from ..core.results import AnnotationResult
from ..utils.file_utils import cmd_exists, safe_run_cmd
from ..utils.sequence_utils import get_combined_seqs, write_temp_file

try:
    from Bio.Alphabet import IUPAC
except ImportError:
    IUPAC = None

class ViralGenomeTopology(BaseAnnotator):
    """
    Predict viral genome topology (circular or linear) using LASTZ.
    Unique ViRAnPy implementation.
    """
    def __init__(self, config, logger=None):
        super().__init__(config, logger)
        self.topology_results = {}

    def check_dependencies(self) -> bool:
        """Check if LASTZ is available."""
        return cmd_exists("lastz")

    def validate_input(self, input_file: str) -> bool:
        """Validate input FASTA file."""
        return Path(input_file).exists() and Path(input_file).suffix.lower() in ['.fasta', '.fa', '.fna']

    def run(self, input_file: str, **kwargs) -> Dict[str, Any]:
        """
        Run genome topology prediction for all sequences in the input file.
        Returns a dictionary with annotation results.
        """
        result = AnnotationResult(
            annotator_name=self.name,
            input_file=input_file
        )
        try:
            for record in SeqIO.parse(input_file, "fasta"):
                self._analyze_sequence_topology(record)
            result.add_annotation("topology_results", self.topology_results)
            result.add_metadata("total_sequences", len(self.topology_results))
            result.add_metadata("circular_count", sum(1 for t in self.topology_results.values() if t.get("topology") == "circular"))
            result.add_metadata("linear_count", sum(1 for t in self.topology_results.values() if t.get("topology") == "linear"))
            self.logger.info(f"Genome topology prediction completed for {len(self.topology_results)} sequences")
        except Exception as e:
            result.success = False
            result.error_message = str(e)
            self.logger.error(f"Genome topology prediction failed: {e}")
        return result.to_dict()

    def _analyze_sequence_topology(self, record: SeqRecord) -> None:
        """
        Analyze topology for a single viral sequence using LASTZ.
        Results are stored in self.topology_results.
        """
        self.topology_results[record.id] = {"topology": "linear"}
        combined_seqs = get_combined_seqs(record, self.config, IUPAC)
        write_temp_file(combined_seqs)
        try:
            output = subprocess.check_output([
                "lastz", "temporal_circular.fasta", "--self", "--notrivial",
                "--nomirror", "--ambiguous=iupac",
                "--format=general-:start1,end1,start2,end2,score,strand1,strand2,identity,length1"
            ], stderr=subprocess.PIPE, text=True)
            self._parse_lastz_output(output, record)
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"LASTZ failed for {record.id}: {e}")
        finally:
            if os.path.exists("temporal_circular.fasta"):
                os.remove("temporal_circular.fasta")

    def _parse_lastz_output(self, output: str, record: SeqRecord) -> None:
        """
        Parse LASTZ output to determine genome topology.
        Results are stored in self.topology_results.
        """
        for line in output.strip().split('\n'):
            if not line:
                continue
            try:
                start1, end1, start2, end2, _, strand1, strand2, identity, _, length = line.split()
                length = int(length)
                identity = float(fractions.Fraction(identity))
                read_len = getattr(self.config, 'gc_skew_read_length', 101)
                if (strand1 == strand2 and
                    length > 0.4 * read_len and
                    identity > 0.95 and
                    int(start1) < 5 and
                    int(start2) > read_len and
                    int(end1) < read_len and
                    int(end2) > read_len * 2 * 0.9):
                    self.topology_results[record.id]["topology"] = "circular"
                    self.topology_results[record.id]["identity"] = min(identity, self.topology_results[record.id].get("identity", identity))
                    self.topology_results[record.id]["length"] = max(length, self.topology_results[record.id].get("length", length))
            except (ValueError, IndexError) as e:
                self.logger.warning(f"Failed to parse LASTZ line: {line}, error: {e}")

    def get_output_files(self) -> List[str]:
        return ["genome_topology_log.txt"]

def calculate_gc_skew(name: str, length: int, seq: str, gc_skew_window: int, gc_skew_slide: int) -> Tuple[List, List, List]:
    """
    Calculate GC skew and determine origin and terminus regions.
    Args:
        name: Sequence name
        length: Sequence length
        seq: Sequence string
        gc_skew_window: Window size
        gc_skew_slide: Sliding window size
    Returns:
        Tuple of (origin_terminus, skew, cumulative_skew)
    """
    replacements = {'G': 1, 'C': -1, 'A': 0, 'T': 0, 'N': 0}
    gmc = [replacements.get(base, 0) for base in seq]
    gpc = [abs(i) for i in gmc]
    weights = np.ones(gc_skew_window) / gc_skew_window
    gmc = [[i, c] for i, c in enumerate(signal.fftconvolve(gmc, weights, 'same').tolist())]
    gpc = [[i, c] for i, c in enumerate(signal.fftconvolve(gpc, weights, 'same').tolist())]
    skew = [[], []]
    c_skew = [[], []]
    cs = 0
    for i, m in gmc[0::gc_skew_slide]:
        p = gpc[i][1]
        gcs = m/p if p != 0 else 0
        cs += gcs
        skew[0].append(i)
        c_skew[0].append(i)
        skew[1].append(gcs)
        c_skew[1].append(cs)
    c_skew_min = signal.argrelextrema(np.asarray(c_skew[1]), np.less, order=1)[0].tolist()
    c_skew_max = signal.argrelextrema(np.asarray(c_skew[1]), np.greater, order=1)[0].tolist()
    if len(c_skew_min) == 0 or len(c_skew_max) == 0:
        return [False, False], skew, c_skew
    c_skew_min = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_min]
    c_skew_max = [[c_skew[0][i], c_skew[1][i]] for i in c_skew_max]
    closest, farthest = int(length * 0.45), int(length * 0.55)
    pairs = []
    from itertools import product
    for pair in product(*[[c_skew_min, c_skew_max]]):
        tr, pk = sorted(list(pair), key=lambda x: x[1], reverse=False)
        a = min(tr[0] - pk[0], pk[0] - tr[0])
        pt = abs(tr[1] - pk[1])
        if closest <= a <= farthest:
            pairs.append([pt, tr, pk])
    if not pairs:
        return [False, False], skew, c_skew
    pt, tr, pk = max(pairs, key=lambda x: x[0])
    return [tr[0], pk[0]], skew, c_skew 