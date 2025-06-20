#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Logging utilities for ViRAnPy pipeline.
"""

import logging
import sys
from pathlib import Path
from typing import Optional


def setup_logger(
    name: str = "viranpy",
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    console_output: bool = True
) -> logging.Logger:
    """
    Set up a logger for the ViRAnPy pipeline.
    
    Args:
        name: Logger name
        level: Logging level
        log_file: Optional log file path
        console_output: Whether to output to console
        
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Clear existing handlers
    logger.handlers.clear()
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler
    if console_output:
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    
    # File handler
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


# Create default logger instance
logger = setup_logger()


def eprint(*args, **kwargs):
    """
    Print to stderr (equivalent to the original eprint function).
    
    Args:
        *args: Arguments to print
        **kwargs: Keyword arguments for print
    """
    print(*args, file=sys.stderr, **kwargs)


class PipelineLogger:
    """
    Context manager for pipeline logging with timing information.
    """
    
    def __init__(self, step_name: str, logger_instance: Optional[logging.Logger] = None):
        """
        Initialize the pipeline logger.
        
        Args:
            step_name: Name of the pipeline step
            logger_instance: Logger instance to use
        """
        self.step_name = step_name
        self.logger = logger_instance or logger
        self.start_time = None
        
    def __enter__(self):
        """Enter the logging context."""
        import time
        self.start_time = time.time()
        self.logger.info(f"Starting: {self.step_name}")
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit the logging context."""
        import time
        if self.start_time:
            duration = time.time() - self.start_time
            if exc_type is None:
                self.logger.info(f"Completed: {self.step_name} (took {duration:.2f} seconds)")
            else:
                self.logger.error(f"Failed: {self.step_name} after {duration:.2f} seconds")
                self.logger.error(f"Error: {exc_val}") 