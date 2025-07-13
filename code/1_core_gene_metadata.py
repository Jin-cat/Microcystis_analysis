import pandas as pd
import re
from collections import defaultdict

# === [User Editable Paths] ===
gene_content_path = "./results/2_peppan/PEPPAN.PEPPAN.gene_content.csv"  # PEPPAN gene presence/absence matrix
gff_path = "./results/2_peppan/PEPPAN.PEPPAN.gff"                        # PEPPAN GFF with ortholog annotations
output_path = "./results/3_core_genes/core_gene_metadata.tsv"           # Output: core gene metadata
# ==============================

# Step 1: Load and filter strict single-copy core genes
df = pd.read_csv(gene_content_path)
genomes = df.columns[1:]

# Filter rows where all genomes have a gene, and each genome has only one gene (no ';')
strict_df = df[df[genomes].notnull().all(axis=1) & 
               df[genomes].applymap(lambda x: ";" not in str(x)).all(axis=1)]

# Convert to long format
long_df = strict_df.melt(id_vars="Gene", var_name="Genome", value_name="PEPPAN_ID").dropna()

# Step 2: Parse GFF to get ortholog mappings and representative genes
peppan_to_ortholog_map = defaultdict(list)
peppan_to_rep_genome = {}

with open(gff_path) as f:
    for line in f:
        if "ID=PEPPAN_g_" not in line or "inference=ortholog_group:" not in line:
            continue
        fields = line.strip().split("\t")
        genome = fields[0]
        attributes = fields[8]

        # Extract gene ID and ortholog ID
        peppan_match = re.search(r"ID=(PEPPAN_g_\d+);", attributes)
        inference_match = re.search(r"inference=ortholog_group:(OG\d+)", attributes)
        if peppan_match and inference_match:
            gene_id = peppan_match.group(1)
            ortholog_id = inference_match.group(1)
            peppan_to_ortholog_map[gene_id].append((genome, gene_id))
            if gene_id not in peppan_to_rep_genome:
                peppan_to_rep_genome[gene_id] = genome

# Step 3: Annotate metadata with representative gene
long_df["ortholog_genes"] = long_df["PEPPAN_ID"].map(lambda pid: peppan_to_ortholog_map.get(pid, []))
long_df["representative_genome"] = long_df["PEPPAN_ID"].map(lambda pid: peppan_to_rep_genome.get(pid, "NA"))
long_df["representative_gene"] = long_df.apply(
    lambda row: next((gene for genome, gene in row["ortholog_genes"]
                      if genome == row["representative_genome"]), "NA"), axis=1)

# Step 4: Save metadata
long_df[["Gene", "Genome", "PEPPAN_ID", "representative_genome", "representative_gene"]]\
    .to_csv(output_path, sep="\t", index=False)

print(f"[✓] Core gene metadata written to: {output_path}")
