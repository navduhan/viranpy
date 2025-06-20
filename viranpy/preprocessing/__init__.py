#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Preprocessing module for viral genome annotation.

This module provides functionality for:
- Quality control and trimming
- Host removal and taxonomic classification
- Assembly and contig processing
"""

from .quality_control import QualityController
from .host_removal import HostRemover
from .assembly import Assembler
from .taxonomic_classification import TaxonomicClassifier

__all__ = [
    'QualityController',
    'HostRemover', 
    'Assembler',
    'TaxonomicClassifier'
] 