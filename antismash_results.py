# -*- coding: utf-8 -*-
import os
import argparse
from bs4 import BeautifulSoup
from datetime import datetime

def parse_antismash_html(html_file, species_name):
    """Parses an antiSMASH index.html file to extract BGC information."""
    summary_data = []
    with open(html_file, 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f, "html.parser")

    # Find the table body containing the results
    results_table = soup.find("table", id="results-table")
    if not results_table:
        return summary_data
        
    tbody = results_table.find("tbody")
    if not tbody:
        return summary_data

    trs = tbody.find_all("tr")
    for tr in trs:
        tds = tr.find_all('td')
        if len(tds) < 5:
            continue
            
        region = tds[0].text.strip()
        bgc_type = tds[1].text.strip()
        start = tds[2].text.strip()
        end = tds[3].text.strip()
        
        # Handle knownclusterblast results which have more columns
        known_cluster = tds[4].text.strip() if len(tds) > 5 else "N/A"
        
        summary_data.append([
            species_name, region, bgc_type, start, end, known_cluster
        ])
        
    return summary_data

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Parse antiSMASH HTML results to create summary files.")
    parser.add_argument("-i", "--input_dir", required=True, help="Root directory of antiSMASH outputs.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save summary .txt files.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    date_today = datetime.now().strftime("%Y-%m-%d")

    # Header for the summary files
    header = "Species\tRegion\tType\tFrom\tTo\tMost_Similar_Known_Cluster\n"

    for species_folder in os.listdir(args.input_dir):
        species_path = os.path.join(args.input_dir, species_folder)
        if os.path.isdir(species_path):
            html_path = os.path.join(species_path, "index.html")
            if os.path.exists(html_path):
                print(f"Parsing results for: {species_folder}")
                summary_data = parse_antismash_html(html_path, species_folder)
                
                if summary_data:
                    output_file = os.path.join(args.output_dir, f"{species_folder}_summary.txt")
                    with open(output_file, "w", encoding="utf-8") as f_out:
                        f_out.write(f"# antiSMASH result {date_today}\n")
                        f_out.write(f"# Species = {species_folder}\n")
                        f_out.write(header)
                        for row in summary_data:
                            f_out.write("\t".join(row) + "\n")

    print(f"All summaries saved to: {args.output_dir}")

if __name__ == "__main__":
    main()
