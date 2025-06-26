#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Annotation modules for viral genome annotation.
"""

from .genome_shape import ViralGenomeTopology
from .trna_detector import ViralTRNAFinder
from .ncrna_detector import ViralNcRNAFinder
from .crispr_detector import ViralCRISPRFinder
from .gene_predictor import ViralGeneFinder
from .protein_function_predictor import ViralProteinAnnotator
from .taxonomic_classification import TaxonomicClassificationAnnotator

__all__ = [
    'ViralGenomeTopology',
    'ViralTRNAFinder',
    'ViralNcRNAFinder',
    'ViralCRISPRFinder',
    'ViralGeneFinder',
    'ViralProteinAnnotator',
    'TaxonomicClassificationAnnotator'
] 