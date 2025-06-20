#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Core components for ViRAnPy pipeline.
"""

from .base import BaseAnnotator, BasePredictor
from .pipeline import ViralAnnotationPipeline
from .results import AnnotationResult, PredictionResult

__all__ = [
    "BaseAnnotator",
    "BasePredictor", 
    "ViralAnnotationPipeline",
    "AnnotationResult",
    "PredictionResult",
] 