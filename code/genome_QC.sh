#!/bin/bash
set -e

source ~/miniconda3/etc/profile.d/conda.sh

# === Path ===
INPUT_DIR="/data/Microcystis/0_rawdata"
OUTPUT_DIR="/data/Microcystis/QC"
THREADS=16
BUSCO_LINEAGE="chroococcales_odb10"

# === BUSCO ===
conda activate busco
for fasta in "$INPUT_DIR"/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    outdir="$OUTPUT_DIR/$sample"
    mkdir -p "$outdir"

    echo "Running BUSCO for $sample"
    busco -i "$fasta" \
          -l "$BUSCO_LINEAGE" \
          -o "${sample}_busco" \
          -m genome \
          -c $THREADS \
          --out_path "$outdir"
done
conda deactivate

# === CheckM2 ===
conda activate checkm2
for fasta in "$INPUT_DIR"/*.fasta; do
    sample=$(basename "$fasta" .fasta)
    outdir="$OUTPUT_DIR/$sample"
    echo "Running CheckM2 for $sample"

    checkm2 predict \
        --input "$fasta" \
        --output "$outdir/${sample}_checkm2.tsv" \
        --threads $THREADS \
        --models all
done
conda deactivate

echo "All QC complete. Results saved in $OUTPUT_DIR"
