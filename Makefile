# Set default shell.
SHELL := bash

# Run recipe on a single shell.
.ONESHELL:

# Enable bash strict mode.
.SHELLFLAGS := -eu -o pipefail -c

# Delete generated files upon encountering an error.
.DELETE_ON_ERROR:

# Set GNU Make options.
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# Override default prefix for recipes.
.RECIPEPREFIX = >

# Check if all required exectuables are installed.
EXE = gffread hisat2 hisat2-build parallel Rscript
CHECK := $(foreach exec, $(EXE),\
				 $(if $(shell which $(exec)),some string,$(error 'Program: $(exec) not found.')))

# Define path variables.

## This points to the project root directory.
ROOT = ${HOME}/pipelines/rnaseq/practice

## This points to the data directory.
DATADIR = ${ROOT}/data

## This points to the output directory.
OUTDIR = ${ROOT}/output

## This points to the plots directory.
PLOTDIR = ${ROOT}/plots

## URL to complete data
URL ?= http://data.biostarhandbook.com/rnaseq/projects/griffith/griffith-data.tar.gz

## This points to the reference genome (FASTA).
REF ?= ${DATADIR}/refs/22.fa

## This ponts to the reference genome index (ht2).
REFIDX ?= ${REF}%.ht2

## This points to the annotation file (GTF).
GTF ?= ${DATADIR}/refs/22.gtf

## This points to the ids text file.
IDS ?= ${DATADIR}/ids.txt

## This points to the overlap counts.
COUNTS ?= ${OUTDIR}/counts.txt

## This points to the STANDARDIZED overlap counts.
SCOUNTS ?= ${OUTDIR}/counts.csv

## This points to filepath pattern for BAM files.
BAM ?= ${DATA}/bam/%.bam

## This points to filepath pattern for BAM indices.
BAI ?= ${DATA}/bam/%.bai

vars:
> @echo 'Working directory: $(PWD)'
> @echo 'REG: ${REF}'
> @echo 'GTF: ${GTF}'
> @echo 'IDS: ${IDS}'

data:
> # Create data directory.
> mkdir -p ${DATADIR}
>
> # Download the data from source.
> wget -nc ${URL} -O '${DATADIR}/griffith-data.tar.gz'
> 
> # Decompress tar file.
> tar xvzf '${DATADIR}/griffith-data.tar.gz' -C ${DATADIR}

data/refs/transcripts.fa: data
> # Extract transcripts from the reference and annotation files.
> gffread -w '${DATADIR}/refs/transcripts.fa' -g ${REF} ${GTF}

${IDS}:
> # Generate sample ids from either read filenames or design file.
> bash '${ROOT}/scripts/generate_ids.sh'

${REFIDX}: ${REF}
> # Build reference genome index.
> hisat2-build ${REF} ${REF}

${BAM}: ${IDS} ${REFIDX}
> # Create bam directory.
> mkdir -p ${DATADIR}/bam
> # Align reads to the reference genome
> cat ${IDS} | \
>	parallel "hisat2 -x ${REF} -1 ${DATADIR}/reads/{}_R1.fq -2 ${DATADIR}/reads/{}_R2.fq | samtools view -Sbh - | samtools sort - ${DATADIR}/bam/{}"

${BAI}: ${BAM}
> # Generate BAM indices
> cat ${IDS} | parallel 'samtools index ${DATADIR}/bam/{}.bam'

${COUNTS}: ${BAI}
> cat ${IDS} | \
> parallel -j 1 echo '${DATADIR}/bam/{}.bam' | \
> xargs featureCounts -p --countReadPairs \
>		-a ${GTF} \
>		-o ${OUTDIR}/counts.txt

${SCOUNTS}: ${COUNTS}
> # Standardize values in the count matrix.
> Rscript ${ROOT}/scripts/parse_featurecounts.r

deseq2: ${SCOUNTS}
> # Run statistical method on standardized count matrix using DESeq2.
> Rscript ${ROOT}/scripts/deseq2.r ${OUTDIR}/counts.csv ${OUTDIR}/results_deseq2.csv

heatmap: ${OUTDIR}/results_deseq2.csv
> # Create a plots directory.
> mkdir -p ${PLOTDIR}
> # Generate heatmat from DESeq2 output.
> Rscript ${ROOT}/scripts/create_heatmap.r ${OUTDIR}/results_deseq2.csv ${PLOTDIR}/heatmap_deseq2.pdf

open-plot:
> # Open generated heatmap.
> open ${PLOTDIR}/heatmap_deseq2.pdf

# Run RNA-seq pipeline.
all: deseq2

.PHONY: vars open-plot
