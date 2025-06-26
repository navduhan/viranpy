#!/usr/bin/env python3
"""
Test script to verify viral database configuration and HMMER implementation.
"""

import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

from viranpy.config import PipelineConfig

def test_viral_databases():
    """Test viral database configuration."""
    print("Testing viral database configuration...")
    
    # Create config with minimal requirements
    config = PipelineConfig(
        non_crna=True,  # Skip RFAM
        blast_switch=True,  # Use BLAST instead of DIAMOND
        no_hmmer=False,  # Enable HMMER
        no_vfam=False,   # Enable VFAM
        no_rvdb=False,   # Enable RVDB
        no_vogs=False,   # Enable VOGs
        no_phrogs=False  # Enable PHROGS
    )
    
    print(f"VFAM database: {config.vfam_database}")
    print(f"RVDB database: {config.rvdb_database}")
    print(f"VOGs database: {config.vogs_database}")
    print(f"PHROGS database: {config.phrogs_database}")
    
    # Check if databases exist
    databases = {
        'VFAM': config.vfam_database,
        'RVDB': config.rvdb_database,
        'VOGs': config.vogs_database,
        'PHROGS': config.phrogs_database
    }
    
    for name, path in databases.items():
        if path and Path(path).exists():
            print(f"✓ {name} database found: {path}")
        else:
            print(f"✗ {name} database not found or not accessible")
    
    return config

if __name__ == "__main__":
    test_viral_databases() 