#!/bin/bash
set -e

# === [Set up log output] ===
exec > >(tee -i phylogenetic_analysis_log.txt)
exec 2>&1

# ===========================================================================
# Phylogenomic pipeline using PEPPAN strict single-copy core genes
# This script assumes that:
# - MAFFT alignment has already been completed (step 4)
# - FASconCAT-G_v1.06.1.pl has already been run (step 5)
#
# [FASconCAT-G Options Used:]
# - READ IN ALL files:         YES
# - READ IN SINGLE files:      NO
# - FILE HANDLING:             Supermatrix
# - SEQUENCE TRANSLATION:      NO
# - EXCLUDING 3rd POSITION:    Remain
# - RY CODING:                 NO
# - BUILD CONSENSUS:           NO
# - RENAME SEQUENCE NAMES:     NO
# - ABSENT SEQUENCE CODING:    Missing
# ===========================================================================

# Set base directories
BASE_DIR="/data/Microcystis"
FASCONCAT_DIR="$BASE_DIR/5_fasconcat"
MODELTEST_DIR="$BASE_DIR/6_modeltest"
RAXML_DIR="$BASE_DIR/7_raxml"

mkdir -p "$MODELTEST_DIR" "$RAXML_DIR"

# ------------------------------------------------
# Run ModelTest-NG based on FASconCAT-G result
# ------------------------------------------------
echo "Running ModelTest-NG using FASconCAT-G output"
conda activate modeltest-ng
modeltest-ng -p 32 --force \
    -i "$FASCONCAT_DIR/FcC_supermatrix.fas" \
    --partitions "$FASCONCAT_DIR/FcC_partitions.txt" \
    -o "$MODELTEST_DIR"
conda deactivate
echo "ModelTest-NG completed"

# ------------------------------------------------
# Infer maximum likelihood tree with RAxML-NG
# ------------------------------------------------
echo "Running RAxML-NG for phylogenetic tree construction"
conda activate raxml-ng
raxml-ng --all \
    --msa "$FASCONCAT_DIR/FcC_supermatrix.fas" \
    --model GTR+I+G4 \
    --prefix "$RAXML_DIR/PEPPAN_core_gene/tree" \
    --bs-tree 1000
conda deactivate
echo "RAxML-NG completed"
