#!/bin/bash
#BSUB -J 2026-03-23_results
#BSUB -P acc_oscarlr
#BSUB -q premium
#BSUB -n 1
#BSUB -W 24:00
#BSUB -R "rusage[mem=50000]"
#BSUB -R "span[hosts=1]"
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

set -euo pipefail
set -x

# ================================ MODULES ================================
module load micromamba/1.5.3-0 gcc/14.2.0

export JAVA_HOME="/hpc/packages/minerva-rocky9/java/21.0.4/jdk-21.0.4" # set JAVA paths for conda/mamba
export JAVA_LD_LIBRARY_PATH="${JAVA_HOME}/lib"

module load samtools/1.21 blast/2.16.0+

# ================================ DIRECTORIES ================================
oscar_pro="/sc/arion/projects/oscarlr"
scratch="${oscar_pro}/woodrk03/projects/nanopore_sequencing/results/2026-03-23_sequence_matching"
exp_data="${oscar_pro}/woodrk03/projects/nanopore_sequencing/data/2026-03-23_sequence_matching"
work="/sc/arion/work/woodrk03/nanopore_sequencing/results/2026-03-23_sequence_matching"
data_dir="/sc/arion/work/woodrk03/nanopore_sequencing/data/2026-03-23_sequence_matching"
database="/sc/arion/work/woodrk03/databases"
bin="/sc/arion/work/woodrk03/bin"

# ================================ PREVIOUSLY GENERATED DATA USED IN THIS EXPERIMENT ================================
# w = work directory, d = scratch data directory, r = scratch results directory
# none

# ================================ EXPERIMENTAL CODE ================================
mkdir -p "${scratch}"
mkdir -p "${work}"

function generate_sequences_to_blast() {
    mkdir -p "${scratch}/sequences_to_blast"

    # STEP 1: IGHC vs CHR14
    # we are blasting the last 51363 bases of chr14 against the non-ighc region of the ighc chromosome
    # chr14:105428637-105480000
    # ighc:1-51363

    samtools faidx "${database}/reference/franken/reference.fasta" "ighc:1-51363" \
        > "${scratch}/sequences_to_blast/ighc.fa"
    samtools faidx "${database}/reference/franken/reference.fasta" "chr14:105428638-105480000" \
        > "${scratch}/sequences_to_blast/chr14.fa"
}

function blast_sequences_ighc() {
    mkdir -p "${scratch}/blast_sequences_ighc"
    mkdir -p "${scratch}/blast_sequences_ighc/databases"
    mkdir -p "${scratch}/blast_sequences_ighc/databases/ighc" "${scratch}/blast_sequences_ighc/databases/chr14"

    makeblastdb -in "${scratch}/sequences_to_blast/ighc.fa" -dbtype nucl -out "${scratch}/blast_sequences_ighc/databases/ighc/ighc_db"
    makeblastdb -in "${scratch}/sequences_to_blast/chr14.fa" -dbtype nucl -out "${scratch}/blast_sequences_ighc/databases/chr14/chr14_db"

    # chr14 tail as query vs IGHC (tests whether tail aligns to ighc:1–51363)
    blastn -task blastn \
        -query "${scratch}/sequences_to_blast/chr14.fa" -db "${scratch}/blast_sequences_ighc/databases/ighc/ighc_db" \
        -out "${scratch}/blast_sequences_ighc/chr14_vs_ighc.blastn" \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
}

# ==============================================================================
# Run pipeline
# ==============================================================================
generate_sequences_to_blast
blast_sequences_ighc