#!/bin/bash
set -e

source ~/miniconda3/etc/profile.d/conda.sh

BASE_DIR="/data/Microcystis/0_rawdata"
ILLUMINA_DIR="$BASE_DIR/Illumina"
NANOPORE_DIR="$BASE_DIR/Nanopore"

THREADS=16

# Loop over all strains based on R1 fastq file names
for R1 in "$ILLUMINA_DIR"/*_R1.fastq.gz; do
    STRAIN=$(basename "$R1" _R1.fastq.gz)

    POD5_DIR="$NANOPORE_DIR/$STRAIN"
    ILLUMINA_R1="$ILLUMINA_DIR/${STRAIN}_R1.fastq.gz"
    ILLUMINA_R2="$ILLUMINA_DIR/${STRAIN}_R2.fastq.gz"
    FASTQ="${STRAIN}.fastq"
    WORK_DIR="./${STRAIN}_assembly"
    
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"

    echo "[DORADO] Basecalling for $STRAIN"
    conda activate dorado
    dorado basecaller --trim adapters --emit-fastq sup@v5 "$POD5_DIR" > "$FASTQ"
    conda deactivate

    echo "[FLYE] Assembling $STRAIN"
    conda activate flye
    if [[ "$STRAIN" == "FBCC-A68" ]]; then
        flye --threads $THREADS --meta --nano-hq "$FASTQ" --out-dir FLYE
    else
        flye --threads $THREADS --nano-hq "$FASTQ" --out-dir FLYE
    fi
    conda deactivate

    echo "[MEDAKA] Polishing $STRAIN"
    conda activate medaka
    medaka_consensus -i "$FASTQ" -d FLYE/assembly.fasta \
                     -o Medaka -t $THREADS -m r104_e82_400bps_sup_v4.0.0
    conda deactivate

    echo "[PILON] Polishing with Illumina for $STRAIN"
    conda activate pilon
    bowtie2-build Medaka/consensus.fasta pilon_ref
    bowtie2 --no-mixed -x pilon_ref -1 "$ILLUMINA_R1" -2 "$ILLUMINA_R2" -S pilon_aln.sam -p $THREADS
    samtools view -@ $THREADS -bS pilon_aln.sam | samtools sort -@ $THREADS -o pilon_aln.sorted.bam
    samtools index pilon_aln.sorted.bam
    pilon --genome Medaka/consensus.fasta --frags pilon_aln.sorted.bam --output ${STRAIN}_pilon --threads $THREADS
    conda deactivate

    echo "Completed: $STRAIN"
done
