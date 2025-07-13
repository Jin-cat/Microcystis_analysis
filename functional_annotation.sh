#!/bin/bash
set -e

# --- Configuration ---
BASE_DIR="/data2/250401_Microcystis"
PROKKA_DIR="$BASE_DIR/1_prokka"
FUNC_DIR="$BASE_DIR/9_functional_annotation"
THREADS=64
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # Assumes scripts are in the same folder

# --- Initialize Log File ---
exec > >(tee -a functional_annotation_log.txt)
exec 2>&1

echo "Starting Functional Annotation Pipeline"
source /root/miniconda3/etc/profile.d/conda.sh # Adjust path to your conda installation

# --- 1. Functional Annotation with eggNOG-mapper ---
echo "[Step 1/3] Running eggNOG-mapper..."
EGGNOG_OUTPUT="$FUNC_DIR/eggnog"
FAA_DIR="$PROKKA_DIR" # Directory containing .faa files from Prokka
mkdir -p "$EGGNOG_OUTPUT"

conda activate func # Assumes eggnog-mapper is in 'func' env
for strain_dir in "$FAA_DIR"/*/; do
    strain_name=$(basename "$strain_dir")
    faa_file="${strain_dir}/${strain_name}.faa"
    output_dir_strain="$EGGNOG_OUTPUT/$strain_name"
    
    if [ -f "$output_dir_strain/${strain_name}.emapper.annotations" ]; then
        echo "eggNOG results for $strain_name already exist. Skipping."
        continue
    fi
    
    echo "Running eggNOG-mapper for $strain_name"
    mkdir -p "$output_dir_strain"
    emapper.py -i "$faa_file" \
               -o "$strain_name" \
               --output_dir "$output_dir_strain" \
               --tax_scope Bacteria \
               --cpu "$THREADS" \
               -m diamond \
               --data_dir /path/to/eggnog-db # IMPORTANT: Set path to eggNOG database
done
conda deactivate
echo "[Step 1/3] eggNOG-mapper annotation complete."


# --- 2. BGC Detection with antiSMASH ---
echo "[Step 2/3] Running antiSMASH for BGC detection..."
ANTISMASH_OUTPUT="$FUNC_DIR/antismash"
GBK_DIR="$PROKKA_DIR" # Directory containing .gbk files from Prokka
mkdir -p "$ANTISMASH_OUTPUT"

conda activate antismash
for strain_dir in "$GBK_DIR"/*/; do
    strain_name=$(basename "$strain_dir")
    gbk_file="${strain_dir}/${strain_name}.gbk"
    output_dir_strain="$ANTISMASH_OUTPUT/$strain_name"

    if [ -d "$output_dir_strain" ]; then
        echo "antiSMASH results for $strain_name already exist. Skipping."
        continue
    fi

    echo "Running antiSMASH for $strain_name"
    antismash --cpus "$THREADS" \
              --output-dir "$output_dir_strain" \
              --genefinding-tool none \
              "$gbk_file"
done
conda deactivate
echo "[Step 2/3] antiSMASH analysis complete."


# --- 3. Parse antiSMASH Results ---
echo "[Step 3/3] Parsing and summarizing antiSMASH results..."
SUMMARY_DIR="$FUNC_DIR/antismash_summary"
mkdir -p "$SUMMARY_DIR"

conda activate beautifulsoup4 # Assumes beautifulsoup4 is in this env
python3 "$SCRIPT_DIR/parse_antismash_results.py" \
    --input_dir "$ANTISMASH_OUTPUT" \
    --output_dir "$SUMMARY_DIR"

# Merge all individual summaries into one file
cat "$SUMMARY_DIR"/*.txt > "$SUMMARY_DIR/BGC_summary_all.txt"
echo "[Step 3/3] antiSMASH summary created."

echo "Functional annotation pipeline finished."
