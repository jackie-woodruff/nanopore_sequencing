#!/bin/bash
#BSUB -J 2026-03-23_build_reference
#BSUB -P acc_oscarlr
#BSUB -q premium
#BSUB -n 1
#BSUB -W 3:00
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

set -euo pipefail
set -x

# ================================ MODULES ================================
module load micromamba/1.5.3-0

export JAVA_HOME="/hpc/packages/minerva-rocky9/java/21.0.4/jdk-21.0.4" # set JAVA paths for conda/mamba
export JAVA_LD_LIBRARY_PATH="${JAVA_HOME}/lib"

module load samtools/1.21
module load bedtools/2.31.0

# ================================ DIRECTORIES ================================
scratch="/sc/arion/projects/oscarlr/woodrk03/projects/nanopore_adaptive_sampling/results/2026-03-23_build_reference"
exp_data="/sc/arion/projects/oscarlr/woodrk03/projects/nanopore_adaptive_sampling/data/2026-03-23_build_reference"
work="/sc/arion/work/woodrk03/nanopore_adaptive_sampling/results/2026-03-23_build_reference"
database="/sc/arion/work/woodrk03/databases"
bin="/sc/arion/work/woodrk03/bin"

# ================================ PREVIOUSLY GENERATED DATA ================================
# w = work directory, d = data directory, r = results directory
# franken reference:  ${database}/reference/franken/reference.fasta
# hg38 reference:     ${database}/reference/hg38/hg38.fa
# hg19 reference:     ${database}/reference/hg19/hg19.fa

# ================================ EXPERIMENTAL CODE ================================
mkdir -p "${scratch}"
mkdir -p "${work}"

# --- Step 1: Build adaptive_sampling_reference.fasta ---
# Assembles a modified chr14 from franken + ighc + igh + hg38 tail, then appends
# all remaining franken contigs and hg19 chr14 (renamed chr14_hg19).
#
# New chr14 layout (1-based coordinates):
#   chr14:1-105480000        franken chr14 body          (105,480,000 bp)
#   ighc:51364-401741        franken ighc, no overlap    (    350,378 bp)
#   igh:3-1193129            franken igh, skip 2 bp      (  1,193,127 bp)
#   hg38 chr14:106879845-107043718  post-IGH buffer      (    163,874 bp)
#   Total new chr14:                                      107,187,379 bp
#
# Note: ighc:1-51363 overlaps franken chr14:105428637-105480000; chr14 version is kept.
function build_reference {
    mkdir -p "${scratch}/build_reference"

    conda activate basic_tools

    echo "[build_reference] Assembling adaptive_sampling_reference.fasta ..."
    python3 "${work}/build_reference.py" \
        --franken "${database}/reference/franken/reference.fasta" \
        --hg38    "${database}/reference/hg38/hg38.fa" \
        --hg19    "${database}/reference/hg19/hg19.fa" \
        --output  "${scratch}/build_reference/adaptive_sampling_reference.fasta"

    echo "[build_reference] Indexing ..."
    samtools faidx "${scratch}/build_reference/adaptive_sampling_reference.fasta"

    cut -f1,2 "${scratch}/build_reference/adaptive_sampling_reference.fasta.fai" \
        > "${scratch}/build_reference/adaptive_sampling_reference.fasta.chrom.sizes"

    echo "[build_reference] Done. Reference is at ${scratch}/build_reference/adaptive_sampling_reference.fasta"
}

# --- Step 2: Build IG region BED files ---
# Three BED files are generated:
#   ig_regions.bed          exact IG loci
#   ig_regions_50kb.bed     IG loci ± 50 kb
#   ig_regions_62kb.bed     IG loci ± 62 kb

function build_bed_files {
    mkdir -p "${scratch}/build_bed_files"

    echo "[build_bed_files] Writing ig_regions.bed ..."
    cat > "${scratch}/build_bed_files/ig_regions.bed" << 'EOF'
chr14	105480000	107023505	IGH
chr14_hg19	106531320	106569343	IGH_1-8_1-9_SV
chr2	88837161	89340311	IGK-p
chr2	89841997	90280099	IGK-d
chr22	22378774	23423319	IGL
EOF

    chrom_sizes="${scratch}/build_reference/adaptive_sampling_reference.fasta.chrom.sizes"

    echo "[build_bed_files] Writing ig_regions_50kb.bed ..."
    bedtools slop \
        -i "${scratch}/build_bed_files/ig_regions.bed" \
        -g "${chrom_sizes}" \
        -b 50000 \
        > "${scratch}/build_bed_files/ig_regions_50kb.bed"

    echo "[build_bed_files] Writing ig_regions_62kb.bed ..."
    bedtools slop \
        -i "${scratch}/build_bed_files/ig_regions.bed" \
        -g "${chrom_sizes}" \
        -b 62000 \
        > "${scratch}/build_bed_files/ig_regions_62kb.bed"

    echo "[build_bed_files] Done. BED files are at ${scratch}/build_bed_files/"
}

# ================================ EXECUTION ================================
build_reference
build_bed_files
