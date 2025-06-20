#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Utilities for generating metadata for viral metagenomic contigs.
"""

import os
from datetime import datetime
from typing import Dict, Any, Optional, List
from pathlib import Path


class ViralMetagenomicMetadataGenerator:
    """
    Generate appropriate metadata for viral metagenomic contigs.
    
    This class helps users create metadata files when working with
    viral metagenomic data, including segmented viruses.
    """
    
    def __init__(self):
        """Initialize the viral metadata generator."""
        self.default_metadata = {
            'moltype': 'DNA',
            'tech': 'wgs',
            'genetic_code_table': '11',
            'topology': 'linear'
        }
    
    def generate_viral_metadata(self,
                                sample_name: str,
                                viral_family: Optional[str] = None,
                                viral_genus: Optional[str] = None,
                                viral_species: Optional[str] = None,
                                collection_date: Optional[str] = None,
                                location: Optional[str] = None,
                                source: Optional[str] = None,
                                host: Optional[str] = None,
                                segment_info: Optional[str] = None) -> str:
        """
        Generate metadata specifically for viral contigs from metagenomes.
        
        Args:
            sample_name: Name/ID of the viral metagenomic sample
            viral_family: Viral family if known
            viral_genus: Viral genus if known
            viral_species: Viral species if known
            collection_date: Collection date
            location: Geographic location
            source: Sample source
            host: Host organism if known
            segment_info: Segment information for segmented viruses
            
        Returns:
            String containing viral metadata in SeqIn format
        """
        metadata = []
        
        # Viral organism name construction
        organism_parts = []
        if viral_species:
            organism_parts.append(viral_species)
        elif viral_genus:
            organism_parts.append(f"{viral_genus} virus")
        elif viral_family:
            organism_parts.append(f"{viral_family} virus")
        else:
            organism_parts.append("virus")
        
        organism_name = f"uncultured {' '.join(organism_parts)}"
        metadata.append(f"[organism={organism_name}]")
        
        # Strain name with segment info for segmented viruses
        if segment_info:
            strain_name = f"{sample_name}_{segment_info}"
        else:
            strain_name = sample_name
        metadata.append(f"[strain={strain_name}]")
        
        # Basic metadata
        metadata.append(f"[moltype={self.default_metadata['moltype']}]")
        metadata.append(f"[tech={self.default_metadata['tech']}]")
        metadata.append(f"[genetic_code_table={self.default_metadata['genetic_code_table']}]")
        metadata.append(f"[topology={self.default_metadata['topology']}]")
        
        # Optional metadata
        if collection_date:
            metadata.append(f"[collection-date={collection_date}]")
        else:
            current_date = datetime.now().strftime("%Y-%m-%d")
            metadata.append(f"[collection-date={current_date}]")
        
        if location:
            metadata.append(f"[country={location}]")
        
        if source:
            metadata.append(f"[isolation-source={source}]")
        else:
            metadata.append("[isolation-source=metagenome]")
        
        if host:
            metadata.append(f"[host={host}]")
        
        # Viral-specific notes
        metadata.append("[note=metagenomic viral contig]")
        
        if viral_family:
            metadata.append(f"[note=viral family: {viral_family}]")
        if viral_genus:
            metadata.append(f"[note=viral genus: {viral_genus}]")
        if viral_species:
            metadata.append(f"[note=viral species: {viral_species}]")
        if segment_info:
            metadata.append(f"[note=viral segment: {segment_info}]")
        
        return " ".join(metadata)
    
    def generate_taxonomy_based_viral_metadata(self,
                                              sample_name: str,
                                              taxonomy_results: Dict[str, Any],
                                              collection_date: Optional[str] = None,
                                              location: Optional[str] = None,
                                              source: Optional[str] = None,
                                              host: Optional[str] = None,
                                              segment_info: Optional[str] = None) -> str:
        """
        Generate viral metadata based on taxonomic classification results.
        
        Args:
            sample_name: Name/ID of the viral metagenomic sample
            taxonomy_results: Results from Kraken2 or other taxonomic classifier
            collection_date: Collection date
            location: Geographic location
            source: Sample source
            host: Host organism if known
            segment_info: Segment information for segmented viruses
            
        Returns:
            String containing taxonomy-based viral metadata
        """
        metadata = []
        
        # Extract best viral taxonomic classification
        best_viral_taxon = self._get_best_viral_taxonomic_classification(taxonomy_results)
        
        if best_viral_taxon:
            organism_name = f"uncultured {best_viral_taxon['name']}"
            metadata.append(f"[organism={organism_name}]")
            
            # Add taxonomic information as notes
            if best_viral_taxon.get('rank'):
                metadata.append(f"[note=taxonomic rank: {best_viral_taxon['rank']}]")
            if best_viral_taxon.get('taxid'):
                metadata.append(f"[note=taxonomic ID: {best_viral_taxon['taxid']}]")
            if best_viral_taxon.get('percentage'):
                metadata.append(f"[note=classification confidence: {best_viral_taxon['percentage']:.2f}%]")
        else:
            metadata.append("[organism=uncultured virus]")
        
        # Strain name with segment info
        if segment_info:
            strain_name = f"{sample_name}_{segment_info}"
        else:
            strain_name = sample_name
        metadata.append(f"[strain={strain_name}]")
        
        # Basic metadata
        metadata.append(f"[moltype={self.default_metadata['moltype']}]")
        metadata.append(f"[tech={self.default_metadata['tech']}]")
        metadata.append(f"[genetic_code_table={self.default_metadata['genetic_code_table']}]")
        metadata.append(f"[topology={self.default_metadata['topology']}]")
        
        # Optional metadata
        if collection_date:
            metadata.append(f"[collection-date={collection_date}]")
        else:
            current_date = datetime.now().strftime("%Y-%m-%d")
            metadata.append(f"[collection-date={current_date}]")
        
        if location:
            metadata.append(f"[country={location}]")
        
        if source:
            metadata.append(f"[isolation-source={source}]")
        else:
            metadata.append("[isolation-source=metagenome]")
        
        if host:
            metadata.append(f"[host={host}]")
        
        metadata.append("[note=metagenomic viral contig]")
        
        if segment_info:
            metadata.append(f"[note=viral segment: {segment_info}]")
        
        return " ".join(metadata)
    
    def _get_best_viral_taxonomic_classification(self, taxonomy_results: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Extract the best viral taxonomic classification from results."""
        if not taxonomy_results.get("success"):
            return None
        
        best_viral_taxon = None
        highest_percentage = 0
        
        for file_results in taxonomy_results.get("results", {}).values():
            for taxid, info in file_results.get("taxonomy", {}).items():
                # Check if this is a viral classification
                if self._is_viral_taxon(info):
                    if info.get("percentage", 0) > highest_percentage:
                        highest_percentage = info.get("percentage", 0)
                        best_viral_taxon = {
                            "name": info.get("name", ""),
                            "rank": info.get("rank", ""),
                            "taxid": taxid,
                            "percentage": info.get("percentage", 0)
                        }
        
        return best_viral_taxon
    
    def _is_viral_taxon(self, taxon_info: Dict[str, Any]) -> bool:
        """Check if a taxonomic classification is viral."""
        name = taxon_info.get("name", "").lower()
        rank = taxon_info.get("rank", "").lower()
        
        # Check for viral indicators
        viral_indicators = [
            "virus", "viridae", "virales", "virinae", "phage", "bacteriophage",
            "viruses", "viral", "viroid", "satellite"
        ]
        
        for indicator in viral_indicators:
            if indicator in name:
                return True
        
        # Check viral ranks
        viral_ranks = ["viruses", "viral", "virus"]
        if rank in viral_ranks:
            return True
        
        return False
    
    def create_viral_metadata_file(self, 
                                  output_file: str,
                                  sample_name: str,
                                  viral_family: Optional[str] = None,
                                  viral_genus: Optional[str] = None,
                                  viral_species: Optional[str] = None,
                                  **kwargs) -> str:
        """
        Create a viral metadata file for metagenomic contigs.
        
        Args:
            output_file: Path to output metadata file
            sample_name: Name/ID of the viral metagenomic sample
            viral_family: Viral family if known
            viral_genus: Viral genus if known
            viral_species: Viral species if known
            **kwargs: Additional arguments for metadata generation
            
        Returns:
            Path to created metadata file
        """
        metadata = self.generate_viral_metadata(
            sample_name, viral_family, viral_genus, viral_species, **kwargs
        )
        
        with open(output_file, 'w') as f:
            f.write(metadata + '\n')
        
        return output_file
    
    def create_segmented_virus_metadata(self,
                                       sample_name: str,
                                       segments: List[str],
                                       viral_family: Optional[str] = None,
                                       viral_genus: Optional[str] = None,
                                       viral_species: Optional[str] = None,
                                       **kwargs) -> Dict[str, str]:
        """
        Create metadata for segmented viruses.
        
        Args:
            sample_name: Name/ID of the viral metagenomic sample
            segments: List of segment identifiers (e.g., ['segment1', 'segment2'])
            viral_family: Viral family if known
            viral_genus: Viral genus if known
            viral_species: Viral species if known
            **kwargs: Additional arguments for metadata generation
            
        Returns:
            Dictionary mapping segment names to metadata strings
        """
        segment_metadata = {}
        
        for segment in segments:
            segment_metadata[segment] = self.generate_viral_metadata(
                sample_name, viral_family, viral_genus, viral_species,
                segment_info=segment, **kwargs
            )
        
        return segment_metadata
    
    def create_batch_viral_metadata(self,
                                   contigs_file: str,
                                   output_dir: str,
                                   sample_name: str,
                                   taxonomy_results: Optional[Dict[str, Any]] = None,
                                   **kwargs) -> str:
        """
        Create viral metadata for a batch of contigs with potential taxonomic information.
        
        Args:
            contigs_file: Input contigs FASTA file
            output_dir: Output directory
            sample_name: Sample name
            taxonomy_results: Optional taxonomic classification results
            **kwargs: Additional arguments
            
        Returns:
            Path to created metadata file
        """
        output_file = Path(output_dir) / f"{sample_name}_viral_metadata.txt"
        
        if taxonomy_results:
            metadata = self.generate_taxonomy_based_viral_metadata(
                sample_name, taxonomy_results, **kwargs
            )
        else:
            metadata = self.generate_viral_metadata(sample_name, **kwargs)
        
        with open(output_file, 'w') as f:
            f.write(metadata + '\n')
        
        return str(output_file)


def create_viral_metagenomic_metadata(sample_name: str,
                                     output_file: str = "viral_metadata.txt",
                                     viral_family: Optional[str] = None,
                                     viral_genus: Optional[str] = None,
                                     viral_species: Optional[str] = None,
                                     **kwargs) -> str:
    """
    Convenience function to create viral metadata for metagenomic contigs.
    
    Args:
        sample_name: Name/ID of the viral metagenomic sample
        output_file: Output metadata file path
        viral_family: Viral family if known
        viral_genus: Viral genus if known
        viral_species: Viral species if known
        **kwargs: Additional arguments
        
    Returns:
        Path to created metadata file
    """
    generator = ViralMetagenomicMetadataGenerator()
    return generator.create_viral_metadata_file(
        output_file, sample_name, viral_family, viral_genus, viral_species, **kwargs
    ) 