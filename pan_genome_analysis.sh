#!/bin/bash
set -e

# --- Configuration ---
# Set the base directory for the project
BASE_DIR="/data2/250401_Microcystis"
INPUT_DIR="$BASE_DIR/0_rawdata"
OUTPUT_DIR="$BASE_DIR"
THREADS=64

# --- Initialize Log File ---
exec > >(tee -i pan_genome_analysis_log.txt)
exec 2>&1

echo "Starting Pan-Genome and Phylogenetic Analysis Pipeline"

# --- 1. Genome Annotation with Prokka ---
echo "[Step 1/5] Running Prokka for genome annotation..."
PROKKA_DIR="$OUTPUT_DIR/1_prokka"
GFF_DIR="$PROKKA_DIR/gff"
mkdir -p "$PROKKA_DIR" "$GFF_DIR"

source /root/miniconda3/etc/profile.d/conda.sh # Adjust path to your conda installation
conda activate prokka

for GENOME in "$INPUT_DIR"/*.fna; do
  BASENAME=$(basename "$GENOME" .fna)
  echo "Annotating: $BASENAME"
  prokka --outdir "$PROKKA_DIR/$BASENAME" --cpus "$THREADS" --force --prefix "$BASENAME" "$GENOME"
  cp "$PROKKA_DIR/$BASENAME/$BASENAME.gff" "$GFF_DIR/"
done
conda deactivate
echo "[Step 1/5] Prokka annotation complete."


# --- 2. Pan-Genome Analysis with PEPPAN ---
echo "[Step 2/5] Running PEPPAN for pan-genome analysis..."
PEPPAN_DIR="$OUTPUT_DIR/2_peppan"
mkdir -p "$PEPPAN_DIR"
cd "$PEPPAN_DIR" # PEPPAN runs in the current directory

conda activate peppan
# Run PEPPAN on all GFF files
PEPPAN -t "$THREADS" "$GFF_DIR"/*.gff

# Parse PEPPAN results to get the gene presence/absence matrix and core genes
PEPPAN_parser -g "PEPPAN.PEPPAN.gff" -p 1st_parser
PEPPAN_parser -g "PEPPAN.PEPPAN.gff" -p 2nd_parser -m -t -c -a 100 # Extracts core genes found in 100% of strains
conda deactivate
echo "[Step 2/5] PEPPAN analysis complete."


# --- 3. Core Gene Alignment with MAFFT ---
echo "[Step 3/5] Aligning core genes with MAFFT..."
COREGENE_DIR="$PEPPAN_DIR/2nd_parser/core_gene_fastas" # Default PEPPAN output folder
ALIGN_DIR="$OUTPUT_DIR/4_aligned_genes"
mkdir -p "$ALIGN_DIR"
NUM_STRAINS=$(find "$INPUT_DIR" -name "*.fna" | wc -l)
echo "[INFO] Number of strains detected: $NUM_STRAINS"

conda activate mafft
for fasta_file in "$COREGENE_DIR"/*.fasta; do
    base=$(basename "$fasta_file" .fasta)
    output_file="$ALIGN_DIR/${base}_aligned.fasta"
    echo "Aligning: $base"
    mafft --auto --quiet --thread "$THREADS" "$fasta_file" > "$output_file"
done
conda deactivate
echo "[Step 3/5] MAFFT alignment complete."


# --- 4. Concatenate Alignments with FASconCAT-G ---
echo "[Step 4/5] Concatenating alignments with FASconCAT-G..."
CONCAT_DIR="$OUTPUT_DIR/5_fasconcat"
mkdir -p "$CONCAT_DIR"
cd "$ALIGN_DIR" # FASconCAT-G runs on files in the current directory

# Note: Ensure FASconCAT-G script is executable and in your PATH or provide the full path
FASconCAT-G_v1.06.1.pl -s -l -p -n
mv FcC_supermatrix.fas "$CONCAT_DIR/"
cd "$BASE_DIR"
echo "[Step 4/5] FASconCAT-G concatenation complete. Supermatrix created."


# --- 5. Phylogenetic Tree with RAxML-NG ---
echo "[Step 5/5] Building phylogenetic tree with RAxML-NG..."
RAXML_DIR="$OUTPUT_DIR/7_raxml/PEPPAN_core_gene"
mkdir -p "$RAXML_DIR"
SUPERMATRIX_PATH="$CONCAT_DIR/FcC_supermatrix.fas"

conda activate raxml-ng
raxml-ng --all \
         --msa "$SUPERMATRIX_PATH" \
         --model GTR+I+G4 \
         --prefix "$RAXML_DIR/core_gene_tree" \
         --bs-trees 1000 \
         --threads "$THREADS"
conda deactivate
echo "[Step 5/5] RAxML-NG tree construction complete."

echo "Pipeline finished successfully."
