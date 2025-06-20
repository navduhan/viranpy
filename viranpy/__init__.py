#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
ViRAnPy - Viral Metagenomic Analysis Pipeline

A comprehensive pipeline for viral metagenomic analysis including quality control,
host removal, assembly, taxonomic classification, coverage analysis, and annotation.

This pipeline integrates multiple bioinformatics tools into a unified platform
for viral metagenomic research, featuring automated database management and
flexible input options with automatic paired-end detection.

Note: This pipeline was inspired by the concept of viral genome annotation
but has been completely redesigned and implemented as a unique metagenomic
analysis platform with distinct features and capabilities.
"""

__version__ = "0.1.0"
__author__ = "ViRAnPy Development Team"
__email__ = "viranpy@example.com"

# Import main components
from .config import PipelineConfig
from .pipeline import ViralAnnotationPipeline

# Import utilities
from .utils.metagenomic_metadata import create_viral_metagenomic_metadata
from .utils.database_manager import DatabaseManager

__all__ = [
    'PipelineConfig',
    'ViralAnnotationPipeline', 
    'create_viral_metagenomic_metadata',
    'DatabaseManager'
] 