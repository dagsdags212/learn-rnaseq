#!/usr/bin/env bash

# Print commands and exit on error.
set -uex

## Store reference genomw path.
REF='data/refs/22.fa'
GFF='data/refs/22.gtf'

# Setup project directory structure.
mkdir -vp data
mkdir -vp output
mkdir -vp plots

# Store URL where the data is located.
DATA='http://data.biostarhandbook.com/rnaseq/projects/griffith/griffith-data.tar.gz'

# Download the data from source.
wget -nc ${DATA} -O data/griffith-data.tar.gz

# Decompress the data.
tar xvzf data/griffith-data.tar.gz -C data/

# Inspect structure of untarred directory.
tree data/

# Extract transcripts from the 22.gtf file
gffread -w data/refs/transcripts.fa -g data/refs/22.fa data/refs/22.gtf

# Extract unique read ids.
bash scripts/generate_ids.sh

# Build reference genome index.
hisat2-build ${REF} ${REF}

# Align reads to reference index.

## Create folder for storing BAM files.
mkdir -p data/bam

## Align reads to the reference genome.
cat data/ids.txt | \
parallel "hisat2 -x ${REF} -1 data/reads/{}_R1.fq -2 data/reads/{}_R2.fq | samtools view -Sbh - | samtools sort - data/bam/{}"

## Generate BAM indices.
cat data/ids.txt | \
parallel "samtools index data/bam/{}.bam"

# Count alignments from overlaps.
cat data/ids.txt | \
parallel -j 1 echo "data/bam/{}.bam" | \
xargs featureCounts -p --countReadPairs \
    -a ${GFF} \
    -o output/counts.txt

# Standardize values of count matrix.
Rscript scripts/parse_featurecounts.r

# Run statistical method on matrix using DESeq2.
Rscript scripts/deseq2.r

# Generate heatmap from DESeq2 output.
Rscript scripts/create_heatmap.r output/results_deseq2.csv plots/result_deseq2.pdf