#!/usr/bin/env bash
#
# download_test_data.sh
#
# Purpose:
#   Download small example FASTQ files **and** a lightweight reference genome
#   (chr22 only) so users can run the full pipeline quickly.
#
# Output:
#   data/fastq/        – example FASTQs
#   data/reference/    – chr22 FASTA + indexes (bwa-mem2, samtools, picard)
#
# Usage:
#   bash scripts/download_test_data.sh
#

set -euo pipefail

echo "======================================"
echo " Creating directory structure"
echo "======================================"
mkdir -p data/fastq
mkdir -p data/reference

#############################################
# 1) DOWNLOAD EXAMPLE FASTQ FILES
#############################################

echo
echo "======================================"
echo " Downloading example paired-end FASTQs"
echo "======================================"

# These are small Illumina 2×150 bp FASTQs (public demo data).
# Replace with any public FASTQ if desired.
wget -O data/fastq/EXAMPLE_T_R1.fastq.gz \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/004/SRR3907284/SRR3907284_1.fastq.gz"

wget -O data/fastq/EXAMPLE_T_R2.fastq.gz \
  "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/004/SRR3907284/SRR3907284_2.fastq.gz"

echo "Example FASTQs downloaded:"
ls -1 data/fastq

#############################################
# 2) DOWNLOAD MINI REFERENCE GENOME (CHR22)
#############################################
#
# Why chr22?
#   - Very small → indexing is fast
#   - Works perfectly for demo pipelines
#   - Users do NOT need to download 3+ GB GRCh38

echo
echo "======================================"
echo " Downloading chr22 reference FASTA"
echo "======================================"

wget -O data/reference/chr22.fa.gz \
  "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz"

gunzip -f data/reference/chr22.fa.gz

# Standardize the name to genome.fa (expected by your pipeline)
mv data/reference/chr22.fa data/reference/genome.fa

echo "FASTA downloaded:"
ls -lh data/reference/genome.fa

#############################################
# 3) INDEX THE REFERENCE
#############################################

echo
echo "======================================"
echo " Indexing reference for samtools (.fai)"
echo "======================================"
samtools faidx data/reference/genome.fa

echo
echo "======================================"
echo " Creating Picard sequence dictionary (.dict)"
echo "======================================"
picard CreateSequenceDictionary \
    REFERENCE=data/reference/genome.fa \
    OUTPUT=data/reference/genome.fa.dict

echo
echo "======================================"
echo " Indexing reference for bwa-mem2"
echo "======================================"
bwa-mem2 index data/reference/genome.fa

echo
echo "======================================"
echo " Reference genome ready"
echo "======================================"
ls -1 data/reference

#############################################
# 4) COMPLETION MESSAGE
#############################################

echo
echo "======================================"
echo " All example data downloaded!"
echo "======================================"

echo
echo "Run the pipeline using:"
echo
echo "  nextflow run main.nf -profile conda \\"
echo "      --reads \"data/fastq/*_{R1,R2}.fastq.gz\" \\"
echo "      --reference \"data/reference/genome.fa\""
echo
