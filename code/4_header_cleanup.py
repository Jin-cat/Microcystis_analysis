from Bio import SeqIO
from pathlib import Path

# === [User Editable Path] ===
fasta_dir = Path("./results/3_core_genes/representatives")  # Directory containing representative FASTA files
# ==============================

# Step 1: Clean headers for all .fasta files in the directory
for fasta_file in fasta_dir.glob("*.fasta"):
    updated_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Keep only the genome ID before the colon (e.g., "GCA_xxx" from "GCA_xxx:gene")
        genome_id = record.id.split(":")[0]
        record.id = genome_id
        record.name = ""
        record.description = ""
        updated_records.append(record)

    # Step 2: Overwrite original file with updated records
    SeqIO.write(updated_records, fasta_file, "fasta")

print("✅ FASTA headers have been cleaned and updated by genome ID.")
