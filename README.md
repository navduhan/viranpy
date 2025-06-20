# ViRAnPy - Viral Metagenomic Analysis Pipeline

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/navduhan/viranpy/releases)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://github.com/navduhan/viranpy)

**ViRAnPy** is a comprehensive viral metagenomic analysis pipeline designed for viral genome assembly, annotation, taxonomic classification, coverage analysis, and quality assessment. Built for researchers working with viral metagenomic data, including both DNA and RNA viruses, as well as segmented viruses.

**Primary Workflow**: ViRAnPy is optimized for viral metagenomic analysis starting from raw sequencing reads, providing a complete pipeline from quality control through assembly to viral annotation. For users with pre-assembled contigs, direct annotation is also supported.

**Key Innovations**: ViRAnPy features automated database management, flexible input options with automatic paired-end detection, comprehensive reporting, and seamless integration of multiple bioinformatics tools into a unified viral metagenomic analysis platform.

## 🧬 Key Features

### 🦠 Viral-Specific Analysis
- **Viral taxonomy integration**: Support for viral family, genus, and species classification
- **DNA and RNA virus support**: Handle both DNA and RNA viral genomes
- **Segmented virus support**: Handle multi-segment viruses (influenza, rotavirus, etc.)
- **Host information**: Include host organism metadata in annotations
- **Viral-specific gene prediction**: Optimized for viral genome characteristics

### 📊 Complete Metagenomic Pipeline
- **Quality Control**: FastQC and Trim Galore for read quality assessment
- **Host Removal**: Bowtie2-based host genome filtering
- **Assembly**: SPAdes/MEGAHIT hybrid assembly optimized for viral genomes
- **Taxonomic Classification**: Kraken2 integration for viral taxonomy
- **Coverage Analysis**: BWA/samtools-based coverage calculation
- **Assembly Statistics**: QUAST integration for quality assessment
- **Gene Annotation**: Viral-specific gene prediction and functional annotation

### 🔧 Automated Database Management
- **One-command installation**: `viranpy --build-databases` installs all required databases
- **Automatic path detection**: No need to specify database paths manually
- **Database verification**: `viranpy --check-databases` verifies installation
- **Comprehensive databases**: RFAM, DIAMOND, BLAST, VOGS, RVDB, PHROGS, VFAM

### 📈 Advanced Reporting
- **Interactive HTML reports**: Comprehensive analysis reports with visualizations
- **Coverage statistics**: Detailed coverage analysis per contig
- **Taxonomic classification**: Viral taxonomy with confidence scores
- **Assembly quality**: QUAST metrics and viral-specific statistics
- **TSV outputs**: Machine-readable results for downstream analysis

## 🚀 Quick Start

### Installation

#### Option 1: Conda Environment (Recommended)
```bash
# Clone the repository
git clone https://github.com/viranpy/viranpy.git
cd viranpy

# Create conda environment with all dependencies
conda env create -f environment.yml

# Activate the environment
conda activate viranpy

# Install ViRAnPy in development mode
pip install -e .
```

#### Option 2: Manual Installation
```bash
# Clone the repository
git clone https://github.com/viranpy/viranpy.git
cd viranpy

# Install Python dependencies
pip install -r requirements.txt

# Install bioinformatics tools (using conda or your system package manager)
conda install -c bioconda bwa samtools bowtie2 spades megahit cd-hit kraken2 quast prodigal blast diamond hmmer infernal fastqc trim-galore

# Install ViRAnPy
pip install -e .
```

### Simplified Command-Line Interface

ViRAnPy features a streamlined command-line interface designed for viral metagenomic analysis:

- **Automated database management**: No need to specify database paths manually
- **Read-based workflow**: Optimized for viral metagenomic analysis from raw reads
- **Flexible input options**: Use `--pe1/--pe2` for explicit paired-end, `--single` for single-end, or `--reads` for automatic detection
- **Automatic paired-end detection**: When using `--reads`, automatically detects paired-end (2 files) vs single-end (1 file)
- **One-command database installation**: `viranpy --build-databases`
- **Automatic path detection**: All database paths are automatically configured
- **Clear workflow separation**: Distinct options for metagenomic analysis vs. direct annotation

### Basic Usage

#### Database Management
```bash
# Install all required databases (first time setup)
viranpy --build-databases

# Check if all databases are installed and ready
viranpy --check-databases
```

#### Generate Viral Metadata
```bash
# Basic viral contigs
viranpy --generate-metadata \
    --sample-name VIROME001 \
    --viral-family Myoviridae \
    --viral-genus T4virus \
    --source gut \
    --location USA \
    --host "Escherichia coli"

# Segmented viruses (e.g., influenza)
viranpy --generate-metadata \
    --sample-name INFLUENZA001 \
    --viral-family Orthomyxoviridae \
    --viral-genus Influenzavirus \
    --segments segment1 segment2 segment3 segment4 segment5 segment6 segment7 segment8 \
    --source human \
    --host "Homo sapiens"
```

#### Run Complete Viral Metagenomic Pipeline (RECOMMENDED)
```bash
# Full pipeline with paired-end reads (explicit)
viranpy --pe1 R1.fastq --pe2 R2.fastq \
    --generate-metadata \
    --sample-name VIROME001 \
    --assemble \
    --host-genome host.fasta \
    --taxonomy-contigs \
    --coverage-analysis \
    --quast-analysis

# Full pipeline with paired-end reads (automatic detection)
viranpy --reads R1.fastq R2.fastq \
    --generate-metadata \
    --sample-name VIROME001 \
    --assemble \
    --host-genome host.fasta \
    --taxonomy-contigs \
    --coverage-analysis \
    --quast-analysis

# Full pipeline with single-end reads (automatic detection)
viranpy --reads reads.fastq \
    --viral-metadata viral_metadata.txt

# Full pipeline with comprehensive reporting
viranpy --pe1 R1.fastq --pe2 R2.fastq \
    --viral-metadata viral_metadata.txt

# By default, ViRAnPy will:
# - Run quality control, host removal (if host genome/index provided), assembly, taxonomy, coverage, QUAST, and reporting
# - Manage all required databases automatically (no need to specify database paths)
# - Only skip a step if you use --skip-host-removal, --skip-taxonomy, --skip-qc, etc.

# To skip host removal (for example):
viranpy --reads R1.fastq R2.fastq --viral-metadata viral_metadata.txt --skip-host-removal

# To run only quality control:
viranpy --reads R1.fastq R2.fastq --viral-metadata viral_metadata.txt --qc-only

# To run only assembly:
viranpy --reads R1.fastq R2.fastq --viral-metadata viral_metadata.txt --assemble-only
```

#### Advanced: Direct Annotation of Pre-assembled Contigs
```bash
# For users who already have assembled contigs
viranpy --input genome.fasta --viral-metadata viral_metadata.txt
```

## 📋 Analysis Components

### Complete Metagenomic Pipeline (Default Workflow)

#### 1. Quality Control
- **FastQC**: Read quality assessment and reports
- **Trim Galore**: Adapter trimming and quality filtering
- **Statistics**: Raw vs cleaned read counts, quality metrics

#### 2. Host Removal
- **Bowtie2**: Fast host genome alignment
- **Statistics**: Host vs viral read percentages
- **Output**: Host-free reads for assembly

#### 3. Assembly
- **SPAdes**: High-quality assembly
- **MEGAHIT**: Fast metagenomic assembly
- **CD-HIT**: Hybrid assembly with redundancy removal
- **Statistics**: Contig counts, N50, GC content

#### 4. Taxonomic Classification
- **Kraken2**: Viral taxonomy classification
- **Confidence scores**: Classification reliability
- **Viral filtering**: Focus on viral contigs
- **Output**: Taxonomic hierarchy (phylum to species)

#### 5. Coverage Analysis
- **BWA**: Read mapping to contigs
- **Samtools**: Coverage depth calculation
- **Statistics**: Mean coverage, coverage breadth, depth distribution
- **Output**: Coverage per contig with read counts

#### 6. Assembly Quality (QUAST)
- **QUAST**: Assembly quality assessment
- **Viral metrics**: Contig length distribution, GC content
- **Quality assessment**: Assembly completeness and fragmentation
- **Output**: Detailed assembly statistics

#### 7. Gene Annotation
- **Prodigal**: Gene prediction optimized for viral genomes
- **Functional annotation**: BLAST/DIAMOND/HMMER integration with automated database management
- **Viral databases**: VOGS, RVDB, PHROGS, VFAM (automatically configured)
- **Output**: GenBank, GFF, CSV files

### Advanced: Direct Annotation Workflow
- **Pre-assembled contigs**: Direct annotation of already assembled viral contigs
- **Same annotation pipeline**: Uses identical gene prediction and functional annotation
- **Automated database management**: All databases automatically detected and configured

## 📁 Output Files

### Analysis Results
- `coverage_report.tsv`: Coverage statistics per contig
- `quast_results/`: QUAST assembly quality reports
- `taxonomy_coverage.tsv`: Combined taxonomy and coverage data

### Comprehensive Reports
- `VIROME001_comprehensive_report.html`: Interactive HTML report
- `assembly_report.txt`: Detailed assembly statistics
- `viral_annotation_results/`: Annotation files (GenBank, GFF, CSV)

### Pipeline Logs
- `viranpy.log`: Pipeline execution log
- `quality_control/`: FastQC reports
- `host_removal/`: Bowtie2 alignment statistics
- `assembly/`: Assembly logs and statistics

## 🐍 Python API

### Basic Usage

```python
from viranpy.utils import create_viral_metagenomic_metadata

# Generate viral metadata
metadata_file = create_viral_metagenomic_metadata(
    sample_name="VIROME001",
    viral_family="Myoviridae",
    viral_genus="T4virus",
    source="gut",
    location="USA",
    host="Escherichia coli"
)
```

### Metagenomic Pipeline

```python
from viranpy.core.pipeline import ViralAnnotationPipeline
from viranpy.config import PipelineConfig

# Create configuration for metagenomic analysis
config = PipelineConfig(
    input_file="dummy.fasta",  # Not used for read-based workflow
    viral_metadata_file="viral_metadata.txt"
)

# Create pipeline
pipeline = ViralAnnotationPipeline(config, logger)

# Run full metagenomic pipeline (paired-end automatically detected)
pipeline.run_preprocessing_pipeline(
    read_files=["R1.fastq", "R2.fastq"]  # 2 files = paired-end
)

# Run with single-end reads (automatically detected)
pipeline.run_preprocessing_pipeline(
    read_files=["reads.fastq"]  # 1 file = single-end
)
```

### Coverage Analysis

```python
from viranpy.utils.coverage_calculator import CoverageCalculator

calculator = CoverageCalculator(config, logger)
coverage_results = calculator.calculate_coverage(
    contigs_file="assembled_contigs.fasta",
    read_files=["R1.fastq", "R2.fastq"]  # Paired-end automatically detected
)
```

### Assembly Statistics

```python
from viranpy.utils.assembly_stats import AssemblyStatsCalculator

stats_calculator = AssemblyStatsCalculator(config, logger)
quast_results = stats_calculator.run_quast_analysis(
    contigs_file="assembled_contigs.fasta",
    min_contig_length=200
)
```

### Comprehensive Reporting

```python
from viranpy.utils.comprehensive_reporter import ComprehensiveReporter

reporter = ComprehensiveReporter(config, logger)
html_report = reporter.generate_comprehensive_report(
    output_dir="reports",
    pipeline_results=results,
    sample_name="VIROME001"
)
```

## 🔧 Dependencies

### Required Tools
- **BWA**: Read mapping and alignment
- **Samtools**: BAM file processing and coverage calculation
- **QUAST**: Assembly quality assessment
- **Bowtie2**: Host removal
- **SPAdes/MEGAHIT**: Assembly
- **Kraken2**: Taxonomic classification
- **Prodigal**: Gene prediction
- **BLAST/DIAMOND**: Functional annotation
- **HMMER**: Protein family analysis

### Python Dependencies
- **Biopython**: Sequence processing
- **Pandas**: Data manipulation
- **NumPy**: Numerical computing
- **Matplotlib/Seaborn**: Visualization (for reports)

### Databases
ViRAnPy automatically manages all required databases. Use the database management commands to install and verify databases:

```bash
# Install all required databases (first time setup)
viranpy --build-databases

# Check if all databases are installed and ready
viranpy --check-databases
```

The `--build-databases` command installs the following databases:

#### Annotation Databases
- **RFAM**: ncRNA detection database (Infernal format)
- **RefSeq Viral Proteins**: Protein database for BLAST/DIAMOND searches
- **RVDB**: Viral protein families database (HMMER format)
- **VOGS**: Viral orthologous groups database (HMMER format)
- **VFAM**: Viral protein families database (HMMER format)
- **PHROGS**: Phage protein families database (HMMER format)

#### Taxonomic Classification
- **Kraken2 Viral Database**: Viral-only taxonomic classification database (0.6 GB)
- **Krona Taxonomy**: Taxonomy data for Krona visualization (optional)

#### Database Sources
- **Kraken2 Viral**: [AWS Indexes](https://benlangmead.github.io/aws-indexes/k2) - April 2025 update
- **RFAM**: EBI RFAM database
- **RefSeq Viral**: NCBI RefSeq viral proteins
- **RVDB**: RVDB viral protein database
- **VOGS/VFAM**: VOG database project
- **PHROGS**: PHROGS phage protein database

**Note**: Database paths are automatically detected and managed by ViRAnPy. Users no longer need to specify database paths manually.

## 📖 Documentation

- [Viral Metagenomic Workflow](https://github.com/navduhan/viranpy/blob/main/docs/metagenomic_workflow.md): Comprehensive guide for viral metagenomic analysis
- [Installation Guide](https://github.com/navduhan/viranpy/blob/main/docs/installation.md): Detailed installation instructions
- [API Reference](https://github.com/navduhan/viranpy/blob/main/docs/api.md): Python API documentation
- [Examples](https://github.com/navduhan/viranpy/blob/main/docs/examples.md): Usage examples and tutorials

## 🎯 Use Cases

### Primary: Viral Metagenomics from Raw Reads
- **Viral discovery**: Identify novel viruses in metagenomic samples starting from raw sequencing data
- **Viral diversity**: Characterize viral communities through complete metagenomic analysis
- **Viral abundance**: Quantify viral populations using coverage analysis
- **Functional annotation**: Annotate viral genes and proteins from assembled contigs

### DNA and RNA Viruses
- **DNA viruses**: Bacteriophages, herpesviruses, adenoviruses, etc.
- **RNA viruses**: Influenza, coronaviruses, rotaviruses, etc.
- **Segmented viruses**: Influenza (8 segments), rotavirus (11 segments), etc.
- **Custom segmented viruses**: Support for any number of segments

### Research Applications
- **Environmental virology**: Soil, water, and air viral communities
- **Human virome**: Gut, respiratory, and skin viral communities
- **Animal virome**: Livestock and wildlife viral communities
- **Clinical virology**: Disease-associated viral communities

### Advanced: Direct Annotation
- **Pre-assembled contigs**: Direct annotation of already assembled viral contigs
- **Reference genomes**: Annotation of viral reference genomes
- **Custom analysis**: Specialized annotation workflows

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](https://github.com/navduhan/viranpy/blob/main/CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone the repository
git clone https://github.com/viranpy/viranpy.git
cd viranpy

# Create conda environment with development dependencies
conda env create -f environment.yml
conda activate viranpy

# Install development dependencies
pip install -r requirements-dev.txt
pip install -e .

# Run tests
pytest tests/
```

## 📄 License

This project is licensed under the GPLv3 License - see the [LICENSE](https://github.com/navduhan/viranpy/blob/main/LICENSE) file for details.

## 🔧 External Tools and Databases

ViRAnPy integrates and builds upon several excellent bioinformatics tools and databases:

### Assembly and Quality Control
- **SPAdes**: High-quality genome assembly
- **MEGAHIT**: Fast metagenomic assembly
- **CD-HIT**: Sequence clustering and redundancy removal
- **FastQC**: Read quality assessment
- **Trim Galore**: Adapter trimming and quality filtering
- **QUAST**: Assembly quality assessment

### Read Processing and Alignment
- **BWA**: Read mapping and alignment
- **Bowtie2**: Fast host genome alignment
- **Samtools**: BAM file processing and coverage calculation

### Gene Prediction and Annotation
- **Prodigal**: Gene prediction optimized for viral genomes
- **BLAST**: Sequence similarity search
- **DIAMOND**: Fast protein sequence alignment
- **HMMER**: Protein family analysis
- **Infernal**: RNA structure and alignment

### Taxonomic Classification
- **Kraken2**: Taxonomic classification using k-mers

### Databases
- **RefSeq Viral**: NCBI viral protein database
- **RFAM**: RNA families database
- **VOGS**: Viral orthologous groups
- **RVDB**: Viral protein families
- **VFAM**: Viral protein families
- **PHROGS**: Phage protein families

### Visualization and Reporting
- **Matplotlib/Seaborn**: Data visualization
- **Pandas**: Data manipulation and analysis
- **NumPy**: Numerical computing

We thank the developers and maintainers of these tools for their excellent work in advancing bioinformatics research.

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/navduhan/viranpy/issues)
- **Discussions**: [GitHub Discussions](https://github.com/navduhan/viranpy/discussions)
- **Documentation**: [GitHub Wiki](https://github.com/navduhan/viranpy/wiki)

## 🔬 Citation

If you use ViRAnPy in your research, please cite:

