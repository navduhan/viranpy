#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) Naveen Duhan

"""
Assembly and contig processing for viral genome annotation.
"""

import os
import subprocess
import shutil
from typing import Dict, Any, List, Optional
from pathlib import Path
import psutil

from Bio import SeqIO
from ..utils.file_utils import cmd_exists


class Assembler:
    """
    Assembly and contig processing for viral genomes.
    
    This module handles:
    - SPAdes assembly
    - MEGAHIT assembly
    - CD-HIT clustering
    - Hybrid assembly creation
    """
    
    def __init__(self, config, logger=None):
        """Initialize the assembler."""
        self.config = config
        self.logger = logger or config.logger
        self.assembly_results = {}
    
    def check_dependencies(self) -> bool:
        """Check if required tools are available."""
        required_tools = ["spades.py", "megahit", "cd-hit-est"]
        missing_tools = [tool for tool in required_tools if not cmd_exists(tool)]
        if missing_tools:
            self.logger.error(f"Missing required tools: {', '.join(missing_tools)}")
            return False
        return True
    
    def _get_system_memory_gb(self) -> float:
        """Get system memory in GB."""
        try:
            memory_bytes = psutil.virtual_memory().total
            memory_gb = memory_bytes / (1024**3)
            return memory_gb
        except Exception:
            # Fallback to a conservative default
            return 8.0
    
    def _get_megahit_memory_bytes(self) -> int:
        """Get MEGAHIT memory in bytes based on system capacity."""
        system_memory_gb = self._get_system_memory_gb()
        
        # If user specified memory, use that
        if hasattr(self.config, 'memory') and self.config.memory:
            try:
                # Parse memory string (e.g., "16G", "8GB", "16384M")
                memory_str = str(self.config.memory).upper().strip()
                if memory_str.endswith('G') or memory_str.endswith('GB'):
                    memory_gb = float(memory_str.rstrip('GB'))
                elif memory_str.endswith('M') or memory_str.endswith('MB'):
                    memory_gb = float(memory_str.rstrip('MB')) / 1024
                else:
                    # Assume it's already in GB
                    memory_gb = float(memory_str)
                return int(memory_gb * (1024**3))
            except (ValueError, AttributeError) as e:
                self.logger.warning(f"Invalid memory format: {self.config.memory}. Using system-based default. Error: {e}")
        
        # Default logic based on system memory
        if system_memory_gb >= 16:
            # For systems with 16GB+, use 16GB for MEGAHIT
            default_memory_gb = 16
        else:
            # For smaller systems, use 75% of available memory
            default_memory_gb = max(4, system_memory_gb * 0.75)
        
        memory_bytes = int(default_memory_gb * (1024**3))
        self.logger.info(f"System memory: {system_memory_gb:.1f}GB, MEGAHIT memory: {default_memory_gb:.1f}GB ({memory_bytes:,} bytes)")
        return memory_bytes
    
    def run_spades(self, input_files: List[str], paired: bool = False, 
                   output_dir: str = "spades_assembly") -> Dict[str, Any]:
        """
        Run SPAdes assembly.
        
        Args:
            input_files: List of input FASTQ files
            paired: Whether the reads are paired-end
            output_dir: Output directory for SPAdes results
            
        Returns:
            Dictionary containing SPAdes results
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        # Debug logging
        self.logger.info(f"SPAdes input files: {input_files}")
        self.logger.info(f"SPAdes paired mode: {paired}")
        self.logger.info(f"SPAdes number of input files: {len(input_files)}")
        
        cmd = [
            "spades.py",
            "--meta",
            "-o", output_dir,
            "--threads", str(self.config.ncpus or 1),
            "--memory", str(getattr(self.config, 'spades_memory', 16))
        ]
        
        if paired:
            if len(input_files) >= 2:
                cmd.extend(["-1", input_files[0], "-2", input_files[1]])
                self.logger.info(f"SPAdes paired-end command: R1={input_files[0]}, R2={input_files[1]}")
            else:
                cmd.extend(["-s", input_files[0]])
                self.logger.warning(f"SPAdes paired mode but only 1 file provided: {input_files[0]}")
        else:
            cmd.extend(["-s", input_files[0]])
            self.logger.info(f"SPAdes single-end command: {input_files[0]}")
        
        self.logger.info(f"SPAdes full command: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"SPAdes assembly completed: {output_dir}")
            
            # Get scaffolds file (preferred over contigs for viral assembly)
            scaffolds_file = Path(output_dir) / "scaffolds.fasta"
            contigs_file = Path(output_dir) / "contigs.fasta"
            
            # Prefer scaffolds if available, fallback to contigs
            if scaffolds_file.exists():
                self.logger.info("Using SPAdes scaffolds for assembly")
                assembly_file = scaffolds_file
            elif contigs_file.exists():
                self.logger.info("Using SPAdes contigs (scaffolds not available)")
                assembly_file = contigs_file
            else:
                return {"success": False, "error": "No scaffolds or contigs file found"}
            
            # Rename contigs with descriptive headers
            renamed_file = Path(output_dir) / ("scaffolds_renamed.fasta" if scaffolds_file.exists() else "contigs_renamed.fasta")
            self._rename_contigs_with_descriptive_headers(
                str(assembly_file), 
                str(renamed_file),
                assembly_method="spades"
            )
            return {"success": True, "contigs_file": str(renamed_file)}
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"SPAdes assembly failed: {e}")
            return {"success": False, "error": str(e)}
    
    def run_megahit(self, input_files: List[str], paired: bool = False, 
                    output_dir: str = "megahit_assembly") -> Dict[str, Any]:
        """
        Run MEGAHIT assembly.
        
        Args:
            input_files: List of input FASTQ files
            paired: Whether the reads are paired-end
            output_dir: Output directory for MEGAHIT results
            
        Returns:
            Dictionary containing MEGAHIT results
        """
        # Always write to output_dir inside self.config.root_output
        output_dir = os.path.join(self.config.root_output, output_dir)
        
        # Check if MEGAHIT has already completed successfully
        final_contigs_file = os.path.join(output_dir, "final.contigs.fa")
        renamed_file = os.path.join(output_dir, "contigs_renamed.fasta")
        
        # If we're resuming and the final contigs file exists, use it
        if hasattr(self.config, 'resume') and self.config.resume and os.path.exists(final_contigs_file):
            self.logger.info(f"[RESUME] Found existing MEGAHIT results: {final_contigs_file}")
            # Check if renamed file exists, if not create it
            if not os.path.exists(renamed_file):
                self.logger.info("[RESUME] Creating renamed contigs file from existing MEGAHIT results")
                self._rename_contigs_with_descriptive_headers(
                    final_contigs_file, 
                    renamed_file,
                    assembly_method="megahit"
                )
            return {"success": True, "contigs_file": renamed_file}
        
        # MEGAHIT does not allow pre-existing output dir, so we need to handle this carefully
        if os.path.exists(output_dir):
            # Check if MEGAHIT is currently running (has intermediate files but no final output)
            if os.path.exists(final_contigs_file):
                self.logger.warning(f"MEGAHIT output directory {output_dir} already exists with final results. Using existing results.")
                # Rename contigs if needed
                if not os.path.exists(renamed_file):
                    self._rename_contigs_with_descriptive_headers(
                        final_contigs_file, 
                        renamed_file,
                        assembly_method="megahit"
                    )
                return {"success": True, "contigs_file": renamed_file}
            else:
                # Check if MEGAHIT is still running (has intermediate files)
                intermediate_dir = os.path.join(output_dir, "intermediate_contigs")
                if os.path.exists(intermediate_dir) and any(os.listdir(intermediate_dir)):
                    self.logger.warning(f"MEGAHIT appears to be running in {output_dir}. Waiting for completion...")
                    # Wait a bit and check again
                    import time
                    time.sleep(10)
                    if os.path.exists(final_contigs_file):
                        self.logger.info("MEGAHIT completed while waiting")
                        if not os.path.exists(renamed_file):
                            self._rename_contigs_with_descriptive_headers(
                                final_contigs_file, 
                                renamed_file,
                                assembly_method="megahit"
                            )
                        return {"success": True, "contigs_file": renamed_file}
                    else:
                        self.logger.warning("MEGAHIT still running, will start fresh")
                        shutil.rmtree(output_dir)
                else:
                    self.logger.warning(f"MEGAHIT output directory {output_dir} already exists. Deleting it to avoid MEGAHIT error.")
                    shutil.rmtree(output_dir)
        
        # Do NOT pre-create output_dir for MEGAHIT
        
        # Debug logging
        self.logger.info(f"MEGAHIT input files: {input_files}")
        self.logger.info(f"MEGAHIT input files type: {type(input_files)}")
        self.logger.info(f"MEGAHIT input files length: {len(input_files)}")
        self.logger.info(f"MEGAHIT paired mode: {paired}")
        self.logger.info(f"MEGAHIT number of input files: {len(input_files)}")
        
        # Check if input_files is actually a list and has content
        if not isinstance(input_files, list):
            self.logger.error(f"MEGAHIT input_files is not a list: {type(input_files)}")
            return {"success": False, "error": "input_files must be a list"}
        
        if len(input_files) == 0:
            self.logger.error("MEGAHIT input_files is empty")
            return {"success": False, "error": "No input files provided"}
        
        if len(input_files) == 1:
            self.logger.warning(f"MEGAHIT only one file provided: {input_files[0]}")
            if paired:
                self.logger.warning("MEGAHIT paired mode enabled but only one file provided")
        
        cmd = [
            "megahit",
            "-o", output_dir,
            "-t", str(self.config.ncpus or 1),
            "-m", str(self._get_megahit_memory_bytes())
        ]
        
        # Debug logging for memory
        memory_bytes = self._get_megahit_memory_bytes()
        memory_gb = memory_bytes / (1024**3)
        self.logger.info(f"MEGAHIT memory setting: {memory_gb:.1f}GB ({memory_bytes:,} bytes)")
        
        if paired:
            if len(input_files) >= 2:
                cmd.extend(["-1", input_files[0], "-2", input_files[1]])
                self.logger.info(f"MEGAHIT paired-end command: R1={input_files[0]}, R2={input_files[1]}")
            else:
                cmd.extend(["-r", input_files[0]])
                self.logger.warning(f"MEGAHIT paired mode but only 1 file provided: {input_files[0]}")
        else:
            cmd.extend(["-r", input_files[0]])
            self.logger.info(f"MEGAHIT single-end command: {input_files[0]}")
        
        self.logger.info(f"MEGAHIT full command: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, check=True)
            self.logger.info(f"MEGAHIT assembly completed: {output_dir}")
            
            # Get contigs file
            contigs_file = Path(output_dir) / "final.contigs.fa"
            if contigs_file.exists():
                # Rename contigs with descriptive headers
                renamed_file = Path(output_dir) / "contigs_renamed.fasta"
                self._rename_contigs_with_descriptive_headers(
                    str(contigs_file), 
                    str(renamed_file),
                    assembly_method="megahit"
                )
                return {"success": True, "contigs_file": str(renamed_file)}
            else:
                return {"success": False, "error": "No contigs file found"}
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"MEGAHIT assembly failed: {e}")
            return {"success": False, "error": str(e)}
    
    def run_cdhit(self, input_files: List[str], output_file: str = "cdhit_contigs.fasta",
                  identity: float = 0.95) -> Dict[str, Any]:
        """
        Run CD-HIT clustering to remove redundant contigs.
        
        Args:
            input_files: List of input FASTA files
            output_file: Output file for clustered contigs
            identity: Identity threshold for clustering
            
        Returns:
            Dictionary containing CD-HIT results
        """
        # If only one input file and it exists, use it directly
        if len(input_files) == 1 and os.path.exists(input_files[0]):
            combined_input = input_files[0]
        else:
            # Combine input files if multiple
            combined_input = os.path.join(self.config.root_output, "combined_contigs.fasta")
            self._combine_fasta_files(input_files, combined_input)
        
        try:
            # Change to the output directory before running CD-HIT
            original_cwd = os.getcwd()
            output_dir = os.path.dirname(os.path.abspath(output_file))
            os.chdir(output_dir)
            
            # Always use the absolute path for the input file
            combined_input_abs = os.path.abspath(combined_input)
            output_filename = os.path.basename(output_file)
            
            cmd = [
                "cd-hit-est",
                "-i", combined_input_abs,
                "-o", output_filename,
                "-c", str(identity),
                "-T", str(self.config.ncpus or 1),
                "-M", "0"  # Disable memory limit
            ]
            
            subprocess.run(cmd, check=True)
            
            # Change back to original directory
            os.chdir(original_cwd)
            
            # CD-HIT should create the file with the exact name we specified
            final_output_file = os.path.abspath(output_file)
            if os.path.exists(final_output_file):
                self.logger.info(f"CD-HIT clustering completed: {final_output_file}")
            else:
                self.logger.error(f"CD-HIT output file not found: {final_output_file}")
                return {"success": False, "error": "CD-HIT output file not found"}
            
            # Get statistics
            stats = self._get_cdhit_statistics(final_output_file)
            return {"success": True, "output_file": final_output_file, "statistics": stats}
            
        except subprocess.CalledProcessError as e:
            if 'original_cwd' in locals():
                os.chdir(original_cwd)
            self.logger.error(f"CD-HIT clustering failed: {e}")
            return {"success": False, "error": str(e)}
        finally:
            # Do NOT delete the combined input file
            pass
    
    def create_hybrid_assembly(self, spades_contigs: str, megahit_contigs: str,
                              output_file: str = "hybrid_assembly.fasta") -> Dict[str, Any]:
        """
        Create hybrid assembly from SPAdes and MEGAHIT contigs.
        
        Args:
            spades_contigs: SPAdes contigs file (should be the renamed file)
            megahit_contigs: MEGAHIT contigs file (should be the renamed file)
            output_file: Output file for hybrid assembly
            
        Returns:
            Dictionary containing hybrid assembly results
        """
        # Create hybrid_assembly directory
        if hasattr(self.config, 'root_output') and self.config.root_output:
            output_dir = os.path.join(self.config.root_output, "hybrid_assembly")
        else:
            output_dir = "hybrid_assembly"
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Update output file path to be inside the directory
        output_basename = os.path.basename(output_file)
        final_output_file = os.path.join(output_dir, output_basename)
        
        try:
            # Combine contigs from both assemblers
            all_contigs = []
            
            # Read SPAdes contigs (should be the renamed file)
            if os.path.exists(spades_contigs):
                self.logger.info(f"Reading SPAdes contigs from: {spades_contigs}")
                for record in SeqIO.parse(spades_contigs, "fasta"):
                    # The contigs are already renamed, so we don't need to add "spades_" prefix
                    all_contigs.append(record)
            else:
                self.logger.warning(f"SPAdes contigs file not found: {spades_contigs}")
            
            # Read MEGAHIT contigs (should be the renamed file)
            if os.path.exists(megahit_contigs):
                self.logger.info(f"Reading MEGAHIT contigs from: {megahit_contigs}")
                for record in SeqIO.parse(megahit_contigs, "fasta"):
                    # The contigs are already renamed, so we don't need to add "megahit_" prefix
                    all_contigs.append(record)
            else:
                self.logger.warning(f"MEGAHIT contigs file not found: {megahit_contigs}")
            
            if not all_contigs:
                return {"success": False, "error": "No contigs found from either assembler"}
            
            # Write combined contigs to the hybrid_assembly directory
            temp_combined_file = os.path.join(output_dir, "combined_contigs.fasta")
            SeqIO.write(all_contigs, temp_combined_file, "fasta")
            self.logger.info(f"Combined {len(all_contigs)} contigs from both assemblers: {temp_combined_file}")
            
            # Run CD-HIT to remove redundancy
            # Use the correct output path: output_dir/hybrid_assembly (without .fasta extension for CD-HIT)
            cdhit_output_prefix = os.path.join(output_dir, "hybrid_assembly")
            cdhit_result = self.run_cdhit([temp_combined_file], output_file, identity=0.95)
            
            # Keep the temporary file for debugging
            # if os.path.exists(temp_combined_file):
            #     os.remove(temp_combined_file)
            
            if cdhit_result["success"]:
                # Rename contigs with descriptive headers
                renamed_file = self._rename_contigs_with_descriptive_headers(
                    cdhit_result["output_file"], 
                    final_output_file,
                    assembly_method="hybrid"
                )
                
                self.logger.info(f"Hybrid assembly created: {renamed_file}")
                return {"success": True, "output_file": renamed_file}
            else:
                return cdhit_result
                
        except Exception as e:
            self.logger.error(f"Hybrid assembly failed: {e}")
            return {"success": False, "error": str(e)}
    
    def _rename_contigs_with_descriptive_headers(self, input_file: str, output_file: str, 
                                               assembly_method: str = "unknown") -> str:
        """
        Rename contigs with descriptive headers including contig number, length, GC content, etc.
        
        Args:
            input_file: Input FASTA file
            output_file: Output FASTA file with renamed contigs
            assembly_method: Assembly method used (spades, megahit, hybrid)
            
        Returns:
            Path to the renamed output file
        """
        try:
            renamed_contigs = []
            contig_number = 1
            
            # Check if input file exists and is readable
            if not os.path.exists(input_file):
                self.logger.error(f"Input file does not exist: {input_file}")
                return input_file
            
            for record in SeqIO.parse(input_file, "fasta"):
                try:
                    # Validate record
                    if record is None or not hasattr(record, 'seq') or record.seq is None:
                        self.logger.warning(f"Skipping invalid record: {record.id if hasattr(record, 'id') else 'unknown'}")
                        continue
                    
                    # Calculate contig properties
                    try:
                        length = len(record.seq)
                        if length <= 0:
                            self.logger.warning(f"Skipping contig with zero length: {record.id}")
                            continue
                    except Exception as len_error:
                        self.logger.warning(f"Could not calculate length for contig {record.id}: {len_error}")
                        length = 0
                    
                    # Calculate GC content manually (avoid Bio.SeqUtils.GC for now)
                    gc_content = 0.0
                    try:
                        seq_str = str(record.seq).upper()
                        gc_count = seq_str.count('G') + seq_str.count('C')
                        total_count = len(seq_str)
                        if total_count > 0:
                            gc_content = (gc_count / total_count) * 100
                    except Exception as gc_error:
                        self.logger.warning(f"Could not calculate GC content for contig {record.id}: {gc_error}")
                        gc_content = 0.0
                    
                    # Create descriptive header
                    new_id = f"viranpy_contig_{contig_number:06d}"
                    new_description = f"len={length}bp gc={gc_content:.1f}% method={assembly_method} orig_id={record.id}"
                    
                    # Create new record
                    new_record = record
                    new_record.id = new_id
                    new_record.description = new_description
                    renamed_contigs.append(new_record)
                    
                    contig_number += 1
                    
                except Exception as record_error:
                    self.logger.warning(f"Could not process contig {record.id if hasattr(record, 'id') else 'unknown'}: {record_error}")
                    # Skip this contig and continue with the next one
                    continue
            
            # Write renamed contigs
            if renamed_contigs:
                SeqIO.write(renamed_contigs, output_file, "fasta")
                self.logger.info(f"Renamed {len(renamed_contigs)} contigs with descriptive headers")
                return output_file
            else:
                self.logger.warning("No contigs were successfully processed, copying original file")
                import shutil
                shutil.copy2(input_file, output_file)
                return output_file
            
        except Exception as e:
            self.logger.error(f"Failed to rename contigs: {e}")
            # Fallback: just copy the original file
            import shutil
            shutil.copy2(input_file, output_file)
            return output_file
    
    def filter_contigs_by_length(self, input_file: str, min_length: int = 200,
                                output_file: str = None) -> Dict[str, Any]:
        """
        Filter contigs by minimum length.
        
        Args:
            input_file: Input FASTA file (with descriptive headers)
            min_length: Minimum contig length
            output_file: Output file (if None, auto-generated)
            
        Returns:
            Dictionary containing filtering results
            
        Note:
            This method preserves the descriptive headers created by the assembly methods.
            Headers follow the format: viranpy_contig_NNNNNN len=Xbp gc=Y% method=Z orig_id=W
        """
        # Create output file path
        if output_file is None:
            input_basename = os.path.basename(input_file)
            input_name = os.path.splitext(input_basename)[0]
            output_basename = f"{input_name}_filtered_{min_length}bp.fasta"
        else:
            output_basename = os.path.basename(output_file)
        
        # Save filtered file in the same directory as the input file
        input_dir = os.path.dirname(input_file)
        if not input_dir:
            # If input_file has no directory, use root output directory
            if hasattr(self.config, 'root_output') and self.config.root_output:
                input_dir = self.config.root_output
            else:
                input_dir = "."
        
        final_output_file = os.path.join(input_dir, output_basename)
        
        try:
            filtered_contigs = []
            total_contigs = 0
            
            for record in SeqIO.parse(input_file, "fasta"):
                total_contigs += 1
                if len(record.seq) >= min_length:
                    # Preserve the descriptive header
                    filtered_contigs.append(record)
            
            SeqIO.write(filtered_contigs, final_output_file, "fasta")
            
            stats = {
                "total_contigs": total_contigs,
                "filtered_contigs": len(filtered_contigs),
                "total_bases": sum(len(contig.seq) for contig in filtered_contigs)
            }
            
            self.logger.info(f"Contig filtering completed: {len(filtered_contigs)}/{total_contigs} contigs retained")
            return {"success": True, "output_file": final_output_file, "statistics": stats}
            
        except Exception as e:
            self.logger.error(f"Contig filtering failed: {e}")
            return {"success": False, "error": str(e)}
    
    def _combine_fasta_files(self, input_files: List[str], output_file: str) -> None:
        """Combine multiple FASTA files into one."""
        with open(output_file, 'w') as outfile:
            for input_file in input_files:
                if os.path.exists(input_file):
                    with open(input_file, 'r') as infile:
                        shutil.copyfileobj(infile, outfile)
    
    def _get_cdhit_statistics(self, output_file: str) -> Dict[str, Any]:
        """Get statistics from CD-HIT output."""
        stats = {"total_contigs": 0, "total_bases": 0}
        
        try:
            for record in SeqIO.parse(output_file, "fasta"):
                stats["total_contigs"] += 1
                stats["total_bases"] += len(record.seq)
        except Exception as e:
            self.logger.warning(f"Could not parse CD-HIT output for statistics: {e}")
        
        return stats 