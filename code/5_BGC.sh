#!/bin/bash
set -e

source ~/miniconda3/etc/profile.d/conda.sh

conda activate antismash

input_dir="/data/Microcystis/1_prokka/gbk"
output_dir="/data/Microcystis/9_antismash/output"
mkdir -p "$output_dir"

for gbk in "$input_dir"/*.gbk; do
  strain=$(basename "$gbk" .gbk)
  outdir="$output_dir/$strain"

  if [ -d "$outdir" ]; then
    echo "$strain already processed, skipping..."
    continue
  fi

  echo "Running antiSMASH for $strain"
  antismash -c 16 \
    --output-dir "$outdir" \
    --output-basename "$strain" \
    --hmmdetection-strictness strict --fullhmmer --clusterhmmer \
    --cb-general --cb-knownclusters --pfam2go --smcog-trees \
    --genefinding-tool none "$gbk"
done

conda deactivate


conda activate beautifulsoup4

input_root="/data/Microcystis/9_antismash/output"
summary_out="/data/Microcystis/9_antismash/html"
mkdir -p "$summary_out"

for each_BGC_output in "$input_root"/*; do
    html_file="$each_BGC_output/index.html"
    species_name=$(basename "$each_BGC_output")
    output_file="$summary_out/${species_name}_summary.txt"

    python3 - <<EOF
from bs4 import BeautifulSoup
import re

html_file = "$html_file"
output_file = "$output_file"
species = "$species_name"

with open(html_file, encoding="utf-8") as f:
    soup = BeautifulSoup(f, "html.parser")

species_id = soup.select_one("strong").text if soup.select_one("strong") else "Unknown"
trs = soup.select("tbody tr")

with open(output_file, "w") as out:
    out.write(f"# Species = {species} ({species_id})\n")
    out.write("Species\tRegion\tType\tFrom\tTo\tCluster\tCluster_add_info\tPercentage\n")
    for tr in trs:
        cells = tr.select("td")
        if len(cells) < 7: continue
        out.write(f"{species}\t{cells[0].text.strip()}\t{cells[1].text.strip()}\t{cells[2].text.strip()}\t{cells[3].text.strip()}\t{cells[4].text.strip()}\t{cells[-2].text.strip()}\t{cells[-1].text.strip()}\n")
EOF
done

cat "$summary_out"/*.txt | grep -v "^#" | grep -v "^Species\t" > /data/Microcystis/9_antismash/BGC_summary_cleaned.txt

conda deactivate
