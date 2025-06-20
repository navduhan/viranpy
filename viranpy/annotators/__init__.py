#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Annotation modules for viral genome annotation.
"""

from .genome_shape import GenomeShapeAnnotator
from .trna_detection import TRNADetector
from .ncrna_detection import NCRNADetector
from .crispr_detection import CRISPRDetector
from .gene_prediction import GenePredictor
from .protein_function import ProteinFunctionPredictor
from .taxonomic_classification import TaxonomicClassificationAnnotator

__all__ = [
    'GenomeShapeAnnotator',
    'TRNADetector',
    'NCRNADetector',
    'CRISPRDetector',
    'GenePredictor',
    'ProteinFunctionPredictor',
    'TaxonomicClassificationAnnotator'
] 