#!/bin/bash
set -e

# === [Clean FASTA headers using Biopython] ===
echo "Cleaning FASTA headers to retain only genome IDs..."

python3 - <<EOF
from Bio import SeqIO
from pathlib import Path

fasta_dir = Path("/data/Microcystis/3_strict_single_copy_core_genes")

for fasta_file in fasta_dir.glob("*.fasta"):
    updated_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome_id = record.id.split(":")[0]
        record.id = genome_id
        record.name = ""
        record.description = ""
        updated_records.append(record)
    SeqIO.write(updated_records, fasta_file, "fasta")

print("FASTA headers cleaned successfully.")
EOF

# === [MAFFT alignment] ===
# Set path variables
BASE_DIR="/data/Microcystis"
RAW_DIR="$BASE_DIR/0_rawdata"
COREGENE_DIR="$BASE_DIR/3_strict_single_copy_core_genes"
ALIGN_DIR="$BASE_DIR/4_aligned_genes"
mkdir -p "$ALIGN_DIR"

# Count number of input strains based on raw .fna files
NUM_STRAINS=$(find "$RAW_DIR" -name "*.fna" | wc -l)
echo "[INFO] Number of strains detected: $NUM_STRAINS"

# Activate Conda environment for MAFFT
conda activate mafft

echo "Starting MAFFT alignments..."

for fasta_file in "$COREGENE_DIR"/*.fasta; do
    base=$(basename "$fasta_file" .fasta)
    output_file="$ALIGN_DIR/${base}_aligned.fasta"

    echo "[Aligning] $base"
    mafft --auto --quiet --thread 64 "$fasta_file" > "$output_file"
 
    aligned_count=$(grep -c "^>" "$output_file")
    if [ "$aligned_count" -ne "$NUM_STRAINS" ]; then
        echo "[Warning] $base: Number of aligned sequences is $aligned_count (expected: $NUM_STRAINS)"
    fi
done

conda deactivate

echo "All MAFFT alignments completed successfully."
