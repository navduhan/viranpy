#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Main pipeline module for ViRAnPy.

This module provides the main entry points for the ViRAnPy pipeline.
"""

from .core.pipeline import ViralAnnotationPipeline
from .config import PipelineConfig

__all__ = [
    "ViralAnnotationPipeline",
    "PipelineConfig",
] 