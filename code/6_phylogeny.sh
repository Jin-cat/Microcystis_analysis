#!/bin/bash
set -e

# === [Log all output] ===
exec > >(tee -a phy_log.txt)
exec 2>&1

# === [Initialize conda] ===
source ~/miniconda3/etc/profile.d/conda.sh

# === [User Editable Paths] ===
ALIGNMENT="./results/core_gene_alignment.aln"        # Core gene alignment file
PARTITIONS="./results/partition.txt"                 # Partition file for ModelTest-NG
MODELTEST_DIR="./results/modeltest"
RAXML_DIR="./results/raxml"
THREADS=32
# ==========================

# === [Step 1: Model selection using ModelTest-NG] ===
echo "[+] Running ModelTest-NG for model selection"
mkdir -p "$MODELTEST_DIR"
conda activate modeltest-ng
modeltest-ng -p $THREADS --force \
  -i "$ALIGNMENT" \
  --partitions "$PARTITIONS" \
  -o "$MODELTEST_DIR"
conda deactivate

# === [Step 2: Tree construction using RAxML-NG] ===
echo "[+] Running RAxML-NG for maximum likelihood tree inference"
mkdir -p "$RAXML_DIR"
conda activate raxml-ng
raxml-ng --all \
  --msa "$ALIGNMENT" \
  --model GTR+G \
  --threads $THREADS \
  --bs-trees 100 \
  --prefix "$RAXML_DIR/tree"
conda deactivate

echo "[✓] Phylogenetic analysis completed."
