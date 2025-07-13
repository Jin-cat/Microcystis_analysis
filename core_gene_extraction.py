import pandas as pd
from Bio import SeqIO
import re
from pathlib import Path
from collections import defaultdict
import os
import csv
import sys

# Path setup
# Receive BASE_DIR as the first command line argument
if len(sys.argv) > 1:
    BASE_DIR = Path(sys.argv[1])
else:
    # If no argument, set default based on script location (for manual execution)
    # scripts/02_core_gene_processing/core_genes_SEQ_parsing.py -> Microcystis_Genome_Analysis_Pipeline
    BASE_DIR = Path(__file__).resolve().parent.parent.parent

input_dir = BASE_DIR / "results" / "2_peppan"
output_dir = BASE_DIR / "results" / "3_core_genes"
output_dir.mkdir(parents=True, exist_ok=True)

# File paths
gene_content_path = input_dir / "PEPPAN.PEPPAN.gene_content.csv"
allele_fasta_path = input_dir / "PEPPAN.allele.fna"
gff_path = input_dir / "PEPPAN.PEPPAN.gff"
mapping_output_path = output_dir / "fasta_key_mapping.tsv"
fasta_key_file = output_dir / "fasta_key_mapping.cleaned.tsv" # Use the cleaned version

# Function to remove parentheses and their contents (keep existing code)
def clean_locus_tag(tag):
    # Remove parentheses and content within them
    tag = re.sub(r"\(.*?\)", "", tag)
    return tag

# --- Section related to core gene sequence extraction ---

# 1. Determine strict single-copy core genes (first part of core_genes_SEQ_parsing.py)
print("Identifying strict single-copy core genes...")
strict_core_dict = {}  # gene_name -> [peppan_id1, peppan_id2, ...]
with open(gene_content_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_name = row["Gene"]
        # Get values from columns other than the first column (Gene)
        peppan_ids = [row[col] for col in reader.fieldnames[1:]]

        # If all values are non-empty (i.e., not '-') and do not contain semicolons (;)
        if all(peppan_ids) and all(";" not in v for v in peppan_ids):
            strict_core_dict[gene_name] = peppan_ids

print(f"Found {len(strict_core_dict)} strict single-copy core genes.")

# 2. Read fasta_key_mapping.cleaned.tsv
print("Loading PEPPAN ID → FASTA key mapping...")
peppan_to_key = {}
if not fasta_key_file.exists():
    print(f"Error: {fasta_key_file} not found. Ensure 0_parsing.py or an equivalent script has been run to create this file.")
    sys.exit(1)

with open(fasta_key_file, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        peppan_to_key[row["PEPPAN_ID"]] = row["FASTA_Key"]

# 3. Load PEPPAN.allele.fna (all FASTA keys → SeqRecord)
print("Indexing allele fasta...")
if not allele_fasta_path.exists():
    print(f"Error: {allele_fasta_path} not found. Ensure PEPPAN has been run successfully.")
    sys.exit(1)

allele_index = SeqIO.index(str(allele_fasta_path), "fasta") # SeqIO.index prefers string paths

# 4. Generate FASTA files per gene_name
print("Writing fasta files for each strict core gene...")
for gene_name, peppan_list in strict_core_dict.items():
    records = []
    for peppan_id in peppan_list:
        fasta_key = peppan_to_key.get(peppan_id)
        if fasta_key and fasta_key in allele_index:
            record = allele_index[fasta_key]
            # FASTA header will be cleaned later by 0_re_header.py, so keep original or set simply for now
            # Example: record.id = f"{gene_name}|{record.id}"
            records.append(record)
        else:
            print(f"Warning: FASTA key for PEPPAN ID {peppan_id} not found in mapping or allele fasta. Skipping.")

    if records:
        output_fasta_file = output_dir / f"{gene_name}.fasta"
        SeqIO.write(records, output_fasta_file, "fasta")
    else:
        print(f"Warning: No sequences found for gene {gene_name}. Skipping writing file.")

print("✅ Strict single-copy core gene sequences extraction completed.")
