# ViRAnPy: Viral Metagenomic Analysis Pipeline

---

## Executive Summary

**ViRAnPy** is a modern, publication-ready pipeline for comprehensive viral metagenomic analysis. It is designed for researchers who need robust, reproducible, and richly annotated viral genome analysis from raw sequencing reads to interactive reports. ViRAnPy uniquely supports segmented viruses, viral-specific metadata, and seamless integration with public databases.

---

## Table of Contents
- [ViRAnPy: Viral Metagenomic Analysis Pipeline](#viranpy-viral-metagenomic-analysis-pipeline)
  - [Executive Summary](#executive-summary)
  - [Table of Contents](#table-of-contents)
  - [Quickstart](#quickstart)
  - [Workflow Overview](#workflow-overview)
  - [Stepwise Pipeline Guide](#stepwise-pipeline-guide)
    - [Step 1: Generate Viral Metadata](#step-1-generate-viral-metadata)
      - [Example: Simple Viral Sample](#example-simple-viral-sample)
      - [Example: Segmented Viruses](#example-segmented-viruses)
    - [Step 2: Run the Complete Pipeline](#step-2-run-the-complete-pipeline)
    - [Step 3: Review Results \& Reports](#step-3-review-results--reports)
  - [Python API Usage](#python-api-usage)
    - [Generate Metadata (Function)](#generate-metadata-function)
    - [Generate Metadata for Segmented Viruses (Class)](#generate-metadata-for-segmented-viruses-class)
    - [Coverage Analysis](#coverage-analysis)
    - [Assembly Statistics](#assembly-statistics)
    - [Comprehensive Reporting](#comprehensive-reporting)
  - [Best Practices \& Tips](#best-practices--tips)
  - [Troubleshooting \& FAQ](#troubleshooting--faq)
    - [Common Issues](#common-issues)
    - [Error Messages](#error-messages)
  - [How to Cite ViRAnPy](#how-to-cite-viranpy)
  - [Contact \& Support](#contact--support)

---

## Quickstart

> **Minimal Example:** Analyze a viral metagenome from raw reads to annotated report in two commands.

```bash
# 1. Generate viral metadata (metadata)
viranpy --generate-metadata --sample-name VIROME001 --viral-family Myoviridae --source gut

# 2. Run the full pipeline
viranpy --reads R1.fastq R2.fastq --viral-metadata viral_metadata.txt --assemble --comprehensive-report
```

---

## Workflow Overview

> **ViRAnPy Workflow Diagram:**
>
> [Insert pipeline diagram here: Raw Reads → QC → Host Removal → Assembly → Taxonomy → Coverage → Annotation → Report]

ViRAnPy automates every step of viral metagenomic analysis:
- **Quality Control** (FastQC, Trim Galore)
- **Host Removal** (Bowtie2)
- **Assembly** (MEGAHIT, SPAdes, Hybrid)
- **Taxonomic Classification** (Kraken2)
- **Coverage Analysis** (BWA, samtools)
- **Assembly Quality** (QUAST)
- **Gene Annotation** (Prodigal, viral DBs)
- **Comprehensive Reporting** (HTML, TSV, GenBank)

---

## Stepwise Pipeline Guide

### Step 1: Generate Viral Metadata

ViRAnPy uses "metadata" to attach rich, viral-specific metadata to your contigs. This ensures your results are publication-ready and compatible with public databases.

> **Tip:** Metadata are required for downstream annotation and reporting. You can generate them automatically or provide your own.

#### Example: Simple Viral Sample
```bash
viranpy --generate-metadata \
    --sample-name VIROME001 \
    --viral-family Myoviridae \
    --viral-genus T4virus \
    --source gut \
    --location USA \
    --collection-date 2024-01-15 \
    --host "Escherichia coli"
```
**What this does:**
- Creates a `viral_metadata.txt` file with all the metadata needed for annotation.

#### Example: Segmented Viruses
```bash
viranpy --generate-metadata \
    --sample-name INFLUENZA001 \
    --viral-family Orthomyxoviridae \
    --viral-genus Influenzavirus \
    --viral-species "Influenza A virus" \
    --segments segment1 segment2 segment3 segment4 segment5 segment6 segment7 segment8 \
    --source human \
    --location Canada \
    --host "Homo sapiens"
```
**What this does:**
- Creates a separate metadata file for each segment, e.g., `viral_metadata_segment1.txt`.

---

### Step 2: Run the Complete Pipeline

Run the full analysis pipeline from raw reads to annotated, publication-ready results.

```bash
viranpy --reads R1.fastq R2.fastq \
    --viral-metadata viral_metadata.txt \
    --assemble \
    --host-genome host.fasta \
    --taxonomy-contigs \
    --coverage-analysis \
    --quast-analysis \
    --comprehensive-report
```

**What this does:**
- Performs QC, host removal, assembly, taxonomy, coverage, annotation, and generates interactive reports.

> **ViRAnPy Advantage:**
> - Automatic detection of paired/single-end reads
> - Viral-specific annotation and reporting
> - Segmented virus support

---

### Step 3: Review Results & Reports

After the pipeline completes, you'll find:
- `viral_metadata.txt` (or per-segment files)
- `coverage_report.tsv`, `quast_results/`, `taxonomy_coverage.tsv`
- `VIROME001_comprehensive_report.html` (interactive report)
- Annotation files: GenBank, GFF, CSV
- Pipeline logs: `viranpy.log`

> **Tip:** Open the HTML report in your browser for an interactive summary of your results.

---

## Python API Usage

ViRAnPy provides a Python API for advanced users and batch processing.

### Generate Metadata (Function)
```python
from viranpy.utils.metagenomic_metadata import create_viral_metagenomic_metadata

metadata_file = create_viral_metagenomic_metadata(
    sample_name="VIROME001",
    viral_family="Myoviridae",
    viral_genus="T4virus",
    source="gut",
    location="USA",
    host="Escherichia coli"
)
```

### Generate Metadata for Segmented Viruses (Class)
```python
from viranpy.utils.metagenomic_metadata import ViralMetagenomicMetadataGenerator

generator = ViralMetagenomicMetadataGenerator()
segment_metadata = generator.create_segmented_virus_metadata(
    sample_name="INFLUENZA001",
    segments=["segment1", "segment2", ...],
    viral_family="Orthomyxoviridae",
    viral_genus="Influenzavirus",
    source="human"
)
for segment, metadata in segment_metadata.items():
    with open(f"viral_metadata_{segment}.txt", "w") as f:
        f.write(metadata + "\n")
```

### Coverage Analysis
```python
from viranpy.utils import CoverageCalculator

calculator = CoverageCalculator(config, logger)
coverage_results = calculator.calculate_coverage(
    contigs_file="assembled_contigs.fasta",
    read_files=["R1.fastq", "R2.fastq"],
    paired=True
)
```

### Assembly Statistics
```python
from viranpy.utils import AssemblyStatsCalculator

stats_calculator = AssemblyStatsCalculator(config, logger)
quast_results = stats_calculator.run_quast_analysis(
    contigs_file="assembled_contigs.fasta",
    min_contig_length=200
)
```

### Comprehensive Reporting
```python
from viranpy.utils import ComprehensiveReporter

reporter = ComprehensiveReporter(config, logger)
html_report = reporter.generate_comprehensive_report(
    output_dir="reports",
    pipeline_results=results,
    sample_name="VIROME001"
)
```

---

## Best Practices & Tips

> **Tip:** Always include host and geographic metadata for best downstream compatibility.

- Use viral-specific databases for taxonomy and annotation.
- For segmented viruses, generate and use per-segment metadata.
- Check read quality and host removal efficiency before assembly.
- Use QUAST and coverage analysis to assess assembly quality.
- Review the HTML report for a comprehensive summary.

---

## Troubleshooting & FAQ

### Common Issues

**Q: Low viral coverage?**
- Check host removal efficiency and input read quality.

**Q: Poor assembly?**
- Verify read quality, coverage, and assembly parameters.

**Q: Taxonomic misclassification?**
- Use up-to-date, viral-specific databases.

**Q: Segmented virus handling?**
- Ensure all segments are processed and have metadata.

### Error Messages

```
ERROR: --sample-name is required when using --generate-metadata
```
**Solution:** Add `--sample-name YOUR_SAMPLE_NAME`

```
ERROR: Missing required tools: bwa, samtools
```
**Solution:** Install BWA and samtools for coverage analysis

```
ERROR: QUAST analysis failed
```
**Solution:** Check QUAST installation and contig file format

```
No host genome or Bowtie2 index provided. Skipping host removal. This may result in host contamination in the assembly.
```
**Solution:** No need to add a solution, as this is a warning message and not an error

---

## How to Cite ViRAnPy

If you use ViRAnPy in your research, please cite:

> **ViRAnPy: A Modern Pipeline for Viral Metagenomic Analysis**
> Naveen et al., 2024. [GitHub: https://github.com/naveenbioinfo/viranpy](https://github.com/naveenbioinfo/viranpy)

Please cite: Duhan N. et al., ViRAnPy, https://github.com/naveenduhan/viranpy

---

## Contact & Support

For questions, feature requests, or bug reports, please open an issue on [GitHub](https://github.com/naveenbioinfo/viranpy/issues) or contact the maintainer.

---

**ViRAnPy: Making viral metagenomics easy, reproducible, and publication-ready.** 