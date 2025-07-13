#!/usr/bin/env python

import pandas as pd
from Bio import SeqIO
import re
from pathlib import Path
from collections import defaultdict
import os
import csv

# === [User Editable Paths] ===
input_dir = Path("./results/2_peppan")                 # Directory containing PEPPAN outputs
output_dir = Path("./results/3_core_genes")            # Output directory for FASTA and mappings
output_dir.mkdir(parents=True, exist_ok=True)
# ==============================

# Define input files
gene_content_path = input_dir / "PEPPAN.PEPPAN.gene_content.csv"
allele_fasta_path = input_dir / "PEPPAN.allele.fna"
gff_path = input_dir / "PEPPAN.PEPPAN.gff"
mapping_output_path = output_dir / "fasta_key_mapping.tsv"

# Helper function to clean locus tags (remove parentheses)
def clean_locus_tag(tag):
    return re.sub(r"\(.*?\)", "", tag)

# Step 1: Identify strict single-copy core genes
gene_content = pd.read_csv(gene_content_path, index_col=0)
strict_core = gene_content[
    gene_content.apply(lambda row: row.str.contains(";").sum() == 0 and (row == "-").sum() == 0, axis=1)
]

# Step 2: Map PEPPAN_g_xxx → genome:locus_tag using GFF
peppan_to_fasta_key = defaultdict(dict)

with open(gff_path) as gff_file:
    for line in gff_file:
        if not line.startswith("GCA"):
            continue
        fields = line.strip().split("\t")
        genome = fields[0]
        attributes = fields[8]
        id_match = re.search(r"ID=(PEPPAN_g_\d+);", attributes)
        locus_match = re.search(r"locus_tag=([^;]+)", attributes)
        if id_match and locus_match:
            peppan_id = id_match.group(1)
            locus_tag = clean_locus_tag(locus_match.group(1))
            fasta_key = f"{genome}:{locus_tag}"
            peppan_to_fasta_key[peppan_id][genome] = fasta_key

# Step 3: Load allele FASTA records
allele_records = SeqIO.to_dict(SeqIO.parse(str(allele_fasta_path), "fasta"))

# Step 4: Write FASTA files per strict core gene
with open(mapping_output_path, "w", newline="") as out_map:
    writer = csv.writer(out_map, delimiter="\t")
    writer.writerow(["Genome", "Gene", "FASTA_Header", "Sequence_Length"])

    for gene_id, row in strict_core.iterrows():
        for genome, peppan_id in row.items():
            if peppan_id == "-":
                continue
            fasta_key = peppan_to_fasta_key.get(peppan_id, {}).get(genome)
            if fasta_key and fasta_key in allele_records:
                record = allele_records[fasta_key]
                record.id = f"{genome}_{gene_id}"
                record.description = ""
                out_path = output_dir / f"{record.id}.fasta"
                SeqIO.write(record, out_path, "fasta")
                writer.writerow([genome, gene_id, fasta_key, len(record.seq)])

print(f"[✓] FASTA files written to: {output_dir}")
print(f"[✓] Mapping table saved to: {mapping_output_path}")
