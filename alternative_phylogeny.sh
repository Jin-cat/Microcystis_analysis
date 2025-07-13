#!/bin/bash
set -e

# --- Configuration ---
BASE_DIR="/data2/250401_Microcystis"
# This script assumes a core gene alignment from Roary is available
ROARY_ALIGNMENT="/data2/250401_Microcystis/2_roary/result/core_gene_alignment.aln"
ROARY_PARTITION="/data2/250401_Microcystis/2_roary/result/partition.txt"
THREADS=32

# --- Initialize Log File ---
exec > >(tee -a alternative_phylogeny_log.txt)
exec 2>&1

echo "Starting Alternative Phylogenetic Analysis"
source /root/miniconda3/etc/profile.d/conda.sh # Adjust path to your conda installation

# --- 1. Best-fit Model Selection with ModelTest-NG ---
echo "[Step 1/2] Running ModelTest-NG for model selection..."
MODELTEST_DIR="$BASE_DIR/6_modeltest/"
mkdir -p "$MODELTEST_DIR"

conda activate modeltest-ng
modeltest-ng -p "$THREADS" \
             -i "$ROARY_ALIGNMENT" \
             --partitions "$ROARY_PARTITION" \
             -o "$MODELTEST_DIR/roary_core_gene.model"
conda deactivate
echo "[Step 1/2] ModelTest-NG complete."


# --- 2. Phylogenetic Tree with RAxML-NG ---
echo "[Step 2/2] Building phylogenetic tree with RAxML-NG..."
RAXML_DIR="$BASE_DIR/7_raxml/roary_core_gene"
mkdir -p "$RAXML_DIR"

conda activate raxml-ng
# The model GTR+I+G4 is used here, but you should use the best model
# suggested by ModelTest-NG from the output file in MODELTEST_DIR.
raxml-ng --all \
         --msa "$ROARY_ALIGNMENT" \
         --model GTR+I+G4 \
         --prefix "$RAXML_DIR/tree" \
         --bs-trees 100 \
         --threads "$THREADS"
conda deactivate
echo "[Step 2/2] RAxML-NG tree construction complete."

echo "Alternative pipeline finished."
