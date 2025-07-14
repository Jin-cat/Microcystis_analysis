import pandas as pd
import re
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path
from copy import deepcopy

# === [User Editable Paths] ===
gene_content_path = "/data/Microcystis/2_peppan/strict_core_genes.csv"
gff_path = "/data/Microcystis/2_peppan/PEPPAN.PEPPAN.gff"
allele_fasta_path = "/data/Microcystis/2_peppan/PEPPAN.allele.fna"
output_dir = Path("./data/Microcystis/3_strict_single_copy_core_genes/")
output_dir.mkdir(parents=True, exist_ok=True)

metadata_path = output_dir / "peppan_strict_single_copy_core_genes_metadata.tsv"
# ==============================

### Step 1: Generate metadata for strict single-copy core genes
df = pd.read_csv(gene_content_path)
genomes = df.columns[1:]

# Strict core genes: Present in all genomes and only once in each
strict_df = df[df[genomes].notnull().all(axis=1) & 
               df[genomes].applymap(lambda x: ";" not in str(x)).all(axis=1)]
long_df = strict_df.melt(id_vars="Gene", var_name="Genome", value_name="PEPPAN_ID").dropna()

# Parse GFF: Map PEPPAN ID to gene information
peppan_to_ortholog_map = defaultdict(list)

with open(gff_path) as f:
    for line in f:
        if "ID=PEPPAN_g_" not in line or "inference=ortholog_group:" not in line:
            continue
        fields = line.strip().split("\t")
        contig_full = fields[0]
        attr = fields[8]

        peppan_match = re.search(r"ID=(PEPPAN_g_\d+);", attr)
        inference_match = re.search(r"inference=ortholog_group:([^;]+)", attr)

        if not peppan_match or not inference_match:
            continue

        peppan_id = peppan_match.group(1)
        genome_id = contig_full.split(":")[0]
        contig_id = contig_full.split(":")[1] if ":" in contig_full else contig_full
        original_genome = f"{genome_id}:{contig_id}"
        orthologs = inference_match.group(1).split(",")

        for orth in orthologs:
            parts = orth.split(":")
            if len(parts) < 3:
                continue
            ortho_genome_id = parts[0]
            gene_id = parts[1].split("(")[0]
            allele_part = parts[2]
            if "(" in allele_part:
                allele_part = re.sub(r"\([^)]*\)", "", allele_part)

            allele_match = re.match(r"(t?\d+)", allele_part)
            allele_number = allele_match.group(1) if allele_match else ""

            peppan_to_ortholog_map[peppan_id].append({
                "gene_id": gene_id,
                "allele_number": allele_number,
                "original_genome": original_genome,
                "genome_id": ortho_genome_id
            })

# Create metadata table
records = []

for _, row in long_df.iterrows():
    rep_gene = row["Gene"]
    genome = row["Genome"]
    peppan_id = row["PEPPAN_ID"]

    allele_number = "NA"
    original_genome = "NA"
    rep_genome = "NA"

    if peppan_id in peppan_to_ortholog_map:
        for ortho in peppan_to_ortholog_map[peppan_id]:
            if ortho["gene_id"] == rep_gene:
                rep_genome = ortho["genome_id"]
            if ortho["gene_id"] == rep_gene and ortho["original_genome"].startswith(genome):
                allele_number = ortho["allele_number"]
                original_genome = ortho["original_genome"]

    records.append({
        "representative_genome": rep_genome,
        "representative_gene": rep_gene,
        "allele_number": allele_number,
        "peppan_id": peppan_id,
        "original_genome": original_genome
    })

meta_df = pd.DataFrame(records)
meta_df = meta_df.sort_values(by=["peppan_id", "representative_genome", "representative_gene"])
meta_df.to_csv(metadata_path, sep="\t", index=False)
print(f"[✔] Metadata saved to: {metadata_path}")

### Step 2: Extract representative core gene FASTA sequences
# Load allele sequences from FASTA
allele_dict = {}
for record in SeqIO.parse(allele_fasta_path, "fasta"):
    if record.id not in allele_dict:
        allele_dict[record.id] = record

# Extract sequences for representative genes
gene_to_records = {}
seen_fasta_ids = set()

for _, row in meta_df.iterrows():
    rep_genome = row["representative_genome"]
    rep_gene = row["representative_gene"]
    allele_num = row["allele_number"]
    original_genome = row["original_genome"]

    fasta_key = f"{rep_genome}:{rep_gene}_{allele_num}"
    if fasta_key not in allele_dict:
        continue

    fasta_id = f"{original_genome}_{rep_gene}_{allele_num}"
    if fasta_id in seen_fasta_ids:
        continue
    seen_fasta_ids.add(fasta_id)

    seq = deepcopy(allele_dict[fasta_key])
    seq.id = fasta_id
    seq.name = ""
    seq.description = ""

    if rep_gene not in gene_to_records:
        gene_to_records[rep_gene] = []
    gene_to_records[rep_gene].append(seq)

# Write representative core gene sequences to separate FASTA files
for gene_id, seq_list in gene_to_records.items():
    SeqIO.write(seq_list, output_dir / f"{gene_id}.fasta", "fasta")

print("✅ All representative core gene FASTA files generated.")
print(f"Total representative genes: {len(gene_to_records)}")
