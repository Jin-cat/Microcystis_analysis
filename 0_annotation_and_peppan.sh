#!/bin/bash
set -e

# === [Initialize Conda environment] ===
source ~/miniconda3/etc/profile.d/conda.sh

# === [Set up log output] ===
exec > >(tee -i PEPPAN_core_gene_analysis_log.txt)
exec 2>&1

# === [User-defined directory paths] ===
INPUT_DIR="/data2/250401_Microcystis/0_rawdata"       # Raw genome input (.fna)
OUTPUT_DIR="/data2/250401_Microcystis"                # Base output directory

# === [1. Genome annotation using Prokka] ===
echo "Running Prokka for genome annotation"
PROKKA_DIR="$OUTPUT_DIR/1_prokka"
mkdir -p $PROKKA_DIR

for GENOME in $INPUT_DIR/*.fna; do
  BASENAME=$(basename $GENOME .fna)
  mkdir -p $PROKKA_DIR/$BASENAME
  conda activate prokka
  prokka --outdir $PROKKA_DIR/$BASENAME --cpus 64 --force --prefix $BASENAME $GENOME 
  conda deactivate
done

mkdir -p $PROKKA_DIR/gff
cp $PROKKA_DIR/*/*.gff $PROKKA_DIR/gff/

# === [2. Run PEPPAN: Based on GFF3 and FASTA] ===
echo "Running PEPPAN for genome analysis"
PEPPAN_DIR="$OUTPUT_DIR/2_peppan"
mkdir -p $PEPPAN_DIR
cd $PEPPAN_DIR
conda activate peppan
PEPPAN -t 64 "$PROKKA_DIR/gff/"*.gff
conda deactivate

echo "[✓] Step 1 completed: Prokka annotation + PEPPAN pan-genome analysis"
