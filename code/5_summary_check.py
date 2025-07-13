import pandas as pd

# === [User Editable Paths] ===
metadata_path = "./results/3_core_genes/core_gene_metadata.tsv"
strict_core_path = "./results/2_peppan/PEPPAN.PEPPAN.gene_content.csv"
# ==============================

# Step 1: Load metadata
meta = pd.read_csv(metadata_path, sep="\t", dtype=str)

# Step 2: Count unique representative gene IDs
uniq_genes = meta["representative_gene"].dropna().unique()
print(f"[✓] Number of unique representative genes (non-redundant): {len(uniq_genes)}")

# Step 3: Find PEPPAN_IDs associated with multiple representative genes
dupe_peppan = meta.groupby("PEPPAN_ID")["representative_gene"].nunique()
dupe_peppan = dupe_peppan[dupe_peppan > 1]

if not dupe_peppan.empty:
    print(f"\n[!] PEPPAN_IDs with multiple representative genes:")
    print(dupe_peppan.head())

    print("\n[Examples]")
    for pid in dupe_peppan.index[:3]:
        subset = meta[meta["PEPPAN_ID"] == pid]
        print(subset[["Genome", "PEPPAN_ID", "representative_gene"]].drop_duplicates())
        print("-" * 40)
else:
    print("[✓] No PEPPAN_IDs with multiple representative genes found.")

# Step 4: Cross-check with original strict_core gene list
try:
    strict_core_df = pd.read_csv(strict_core_path, index_col=0)
    strict_gene_ids = set(strict_core_df.index)
    meta_gene_ids = set(meta["Gene"].unique())

    print(f"\n[✓] Number of strict core genes in gene_content.csv: {len(strict_gene_ids)}")
    print(f"[✓] Number of genes in metadata: {len(meta_gene_ids)}")
    print(f"[✓] Number of matching Gene IDs: {len(strict_gene_ids & meta_gene_ids)}")
except Exception as e:
    print(f"[!] Skipped strict_core comparison due to error: {e}")
