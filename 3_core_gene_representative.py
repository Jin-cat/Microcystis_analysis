import pandas as pd
from Bio import SeqIO
from pathlib import Path
from copy import deepcopy

# === [User Editable Paths] ===
metadata_path = "./results/3_core_genes/core_gene_metadata.tsv"       # Metadata with representative_gene info
allele_fasta_path = "./results/2_peppan/PEPPAN.allele.fna"            # PEPPAN allele FASTA
output_dir = Path("./results/3_core_genes/representatives")           # Output directory for filtered FASTA
output_dir.mkdir(parents=True, exist_ok=True)
# ==============================

# Step 1: Load metadata
meta_df = pd.read_csv(metadata_path, sep="\t", dtype=str)
meta_df.columns = [c.strip() for c in meta_df.columns]  # Remove potential whitespace

# Step 2: Load allele FASTA into a dictionary (ignore duplicate keys)
allele_dict = {}
for record in SeqIO.parse(allele_fasta_path, "fasta"):
    if record.id not in allele_dict:
        allele_dict[record.id] = record

# Step 3: Extract representative sequences
written_count = 0
missing_count = 0

for _, row in meta_df.iterrows():
    genome = row["Genome"]
    rep_gene = row["representative_gene"]
    fasta_key = f"{genome}:{rep_gene}"

    if fasta_key in allele_dict:
        record = deepcopy(allele_dict[fasta_key])
        record.id = f"{genome}_{rep_gene}"
        record.name = ""
        record.description = ""
        out_path = output_dir / f"{record.id}.fasta"
        SeqIO.write(record, out_path, "fasta")
        written_count += 1
    else:
        missing_count += 1

print(f"[✓] {written_count} representative gene FASTA files written to: {output_dir}")
if missing_count > 0:
    print(f"[!] {missing_count} records were missing in PEPPAN.allele.fna")
