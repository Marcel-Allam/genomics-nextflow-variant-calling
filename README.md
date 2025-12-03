# Genomics Nextflow Variant Calling Pipeline

Repository: `genomics-nextflow-variant-calling`  
Workflow engine: **Nextflow DSL2**  
Environment: **Conda** (`genomics_env`)

> ⚠️ **IMPORTANT**  
> This pipeline is designed for training and portfolio purposes only.  
> It is **not** validated for clinical decision-making and must **not** be used to guide patient care.

---

## 1. Overview

This repository contains a **Nextflow DSL2 pipeline** that performs:

- Read quality control
- (Optional) read trimming
- Short-read alignment to a reference genome
- Alignment and coverage-based QC
- Variant calling (bcftools)
- Variant-level QC
- Functional annotation (VEP, offline cache, GRCh38)
- Generation of:
  - A **clinical-style QC & variant summary** in Markdown **and HTML**
  - A **MultiQC** dashboard aggregating all QC metrics

The design and structure mimic real-world pipelines used in clinical and translational genomics (e.g. tumor / exome workflows), but packaged as a **teaching and portfolio project**.

---

## 2. Key Features

- 🧬 **Variant calling pipeline**
  - `bwa-mem2` alignment
  - `samtools` for BAM processing and QC
  - `bcftools mpileup + call` for variant calling (germline-style)
- 📊 **Comprehensive QC stack**
  - Raw read QC: **FastQC**
  - Alignment QC: `samtools flagstat`, `samtools stats`, `samtools idxstats`
  - Coverage QC: `samtools coverage`, `samtools depth`
  - Variant QC: `bcftools stats` (on *raw* VCF)
- 🧠 **Clinical-style reporting**
  - Single-sample **Markdown + HTML** report
  - Summary overview (“Section 0”) with:
    - Overall QC status: **PASS / WARN / FAIL**
    - Total reads, mapping %, properly paired %
    - Mean coverage and breadth
    - Total variants & Ts/Tv ratio
    - High-impact variant count
  - Detailed sections for alignment, coverage, variant QC, VEP impact, and tool versions
- 📈 **MultiQC integration**
  - Aggregated QC across:
    - FastQC
    - samtools alignment metrics
    - bcftools stats
    - coverage QC
- ✅ **Version logging**
  - Records versions of Nextflow, samtools, bcftools, and VEP in the report
- 📦 **Reproducible environment**
  - Conda environment (e.g. `envs/genomics_env.yml`) with all required tools
- 🧪 **Offline-friendly**
  - Uses VEP in **offline cache** mode (GRCh38)

---

## 3. Pipeline Workflow

### 3.1 High-level flow

1. **Input FASTQs** are discovered (paired-end) under `data/fastq/*_{R1,R2}.fastq.gz`
2. **FASTQC** runs on raw reads
3. **(Optional) TRIM_READS** using `fastp` (quality-based trimming)
4. **ALIGN_AND_SORT** reads to the reference with `bwa-mem2`, then sort + index BAM
5. **ALIGN_QC** runs `samtools flagstat`, `stats`, and `idxstats`
6. **COVERAGE_QC** runs `samtools coverage` and `depth`
7. **CALL_VARIANTS** uses `bcftools mpileup + call` to produce a raw VCF
8. **VARIANT_QC** runs `bcftools stats` on the **raw VCF**
9. **ANNOTATE_VEP** annotates the VCF using Ensembl VEP (offline cache, GRCh38)
10. **CLINICAL_SUMMARY** combines:
    - alignment QC
    - coverage QC
    - variant QC
    - VEP annotation
    into a **clinical-style report** (`*.md` + `*.html`)
11. **MULTIQC** scans the final `results/` folder and produces an HTML QC dashboard

### 3.2 Mermaid diagram

### 3.2 Mermaid diagram

```mermaid
flowchart LR
    R[Paired FASTQ files (data/fastq/*.fastq.gz)] --> FQC[FASTQC]

    R --> TRIM{Trim reads?}
    TRIM -- "yes" --> T[TRIM_READS (fastp)]
    TRIM -- "no"  --> U[Use raw reads]

    T --> ALN[ALIGN_AND_SORT (bwa-mem2 + samtools sort)]
    U --> ALN

    ALN --> AQC[ALIGN_QC (flagstat / stats / idxstats)]
    ALN --> COV[COVERAGE_QC (coverage / depth)]
    ALN --> VC[CALL_VARIANTS (bcftools mpileup + call)]

    VC --> VQC[VARIANT_QC (bcftools stats on raw VCF)]
    VC --> VEP[ANNOTATE_VEP (VEP offline, GRCh38)]

    AQC --> CS[CLINICAL_SUMMARY (Markdown + HTML)]
    COV --> CS
    VQC --> CS
    VEP --> CS

    CS --> MQC[MULTIQC (multiqc_report.html)]

    FQC --> MQC
    AQC --> MQC
    COV --> MQC
    VQC --> MQC
```

