# **Assembly and Comparative Genome Analysis of *Microcystis* Strains (FBCC-A68, FBCC-A270, and FBCC-A1114)**

This repository contains the complete bioinformatics workflow for the **comparative genomic analysis** of three *Microcystis* strains: *M. aeruginosa* FBCC-A68, *M. ichthyoblabe* FBCC-A1114, and *M. protocystis* FBCC-A270. The analysis pipeline starts from raw sequencing reads and proceeds through genome assembly, annotation, pan-genome analysis, phylogenomics, and secondary metabolite profiling.<br>

> **Note**: This repository contains the analysis scripts used in the study  
> _"The Complete Genomes of *Microcystis ichthyoblabe* Kützing and *Microcystis protocystis* (Crow) Komárek & Anagnostidis Reveal the Complexity and Plasticity of *Microcystis* Genomes"_  
> (Currently under review at *Microorganisms* journal)

---

## **Analysis Workflow**

The analysis is structured as a series of sequential scripts, each performing a distinct stage of the pipeline:

1.  `0_assembly.sh`: Raw read basecalling, hybrid genome assembly, and two-step polishing.
2.  `1_annotation_and_peppan.sh`: Genome annotation and pan-genome analysis.
3.  `2_strict_single_copy_core_genes.py`: Identification and extraction of strict single-copy core genes.
4.  `3_MSA.sh`: Multiple sequence alignment of each core gene.
5.  `4_phylogenetic_tree.sh`: Concatenation of alignments, phylogenetic inference, and ANI calculation.
6.  `5_BGC.sh` & `6_BGC_dot_plot.R`: BGC prediction, summary, and visualization.

---

## **Dependencies & Software**

This workflow relies on several bioinformatics tools and libraries. 

| Category            | Tool/Library      | Version                                |
| :------------------ | :---------------- | :------------------------------------- |
| Assembly & Polishing| Dorado            | `v0.7.1`                               |
|                     | Flye              | `v2.9.5`                               |
|                     | Medaka            | `v2.0.1`                               |
|                     | Pilon             | `v1.24`                                |
|                     | Bowtie2           |                                        |
| Annotation & Pan-Genome| Prokka           | `v1.14.5`                              |
|                     | PEPPAN            | `v1.0.5`                               |
| Phylogenetics       | MAFFT             | `v7.526`                               |
|                     | FASconCAT-G       | `v1.06.1`                              |
|                     | RAxML-NG          | `v1.2.2`                               |
|                     | pyani             | `v0.2.10`                              |
| BGC & Visualization | antiSMASH         | `v7.1.0`                               |
|                     | R (with ggplot2)  | `v3.5.2`                               |
| Genome Quality      | BUSCO             | `v5.8.2` (with `chroococcales_odb10`) |
|                     | CheckM2           | `v1.1.0`                               |
| Python Libraries    | pandas, biopython, beautifulsoup4 |                        |
| Genome Comparison   | Mauve             | `v2.4.0`                               |

---

## **Usage**

### **Step 1: Genome Assembly and Polishing**

This script automates the assembly process from raw reads.

* **Script**: `0_assembly.sh`
* **Description**:
    * Basecalls Nanopore reads using **Dorado** with the `super accuracy model (sup@v5)` and adapter trimming.
    * Assembles genomes using **Flye** with the `--nano-hq` parameter. For *M. aeruginosa* A68, a meta-assembly approach (`--meta`) was used to handle symbiotic microorganisms.
    * Conducts a two-step polishing process: first with **Medaka** using Nanopore reads, then with **Pilon** using Illumina reads mapped by **Bowtie2**.<br><br>

### **Step 2: Annotation and Pan-Genome Analysis**

This script annotates the assembled genomes and initiates the pan-genome comparison.

* **Script**: `1_annotation_and_peppan.sh`
* **Description**:
    * Predicts genes and annotates the polished assemblies using **Prokka**.
    * Runs **PEPPAN** on the Prokka-generated GFF files to perform the pan-genome analysis.<br><br>

### **Step 3: Strict Single-Copy Core Gene Extraction**

This script identifies and extracts the core genes shared by across all strains.

* **Script**: `2_strict_single_copy_core_genes.py`
* **Description**:
    * Parses the `gene_presence_absence.csv` from PEPPAN to identify genes present exactly once in all genomes.
    * Extracts the corresponding nucleotide sequences for these strict single-copy core genes, creating a separate FASTA file for each.<br><br>

### **Step 4: Multiple Sequence Alignment (MSA)**

This script aligns the sequences for each core gene.

* **Script**: `3_MSA.sh`
* **Description**:
    * Iterates through each core gene FASTA file generated in the previous step.
    * Performs a multiple sequence alignment for each gene using **MAFFT**.<br><br>

### **Step 5: Phylogenetic Inference**

This script builds the phylogenetic tree and calculates genomic similarity metrics.

* **Script**: `4_phylogenetic_tree.sh`
* **Description**:
    * Concatenates the aligned core gene sequences into a supermatrix using **FASconCAT-G**.
    * Infers a maximum likelihood phylogeny from the supermatrix using **RAxML-NG** with the `GTR+G` model and 1,000 bootstrap replicates.
    * Calculates pairwise Average Nucleotide Identity (**ANI**) values among genomes using **pyani**.<br><br>

### **Step 6: BGC Analysis and Visualization**

These scripts predict and visualize the distribution of secondary metabolite gene clusters.

* **Scripts**: `5_BGC.sh` and `6_BGC_dot_plot.R`
* **Description**:
    * `5_BGC.sh` runs **antiSMASH** in strict mode to predict BGCs from the annotated GenBank files. It then parses the HTML output to create a unified summary table.
    * `6_BGC_dot_plot.R` reads the summary table and generates a dot plot with **ggplot2** to visualize the presence, type, and similarity confidence of BGCs across all strains.

---

## **Data Availability**

The raw sequencing data and complete genome assemblies of the following *Microcystis* strains have been deposited in the NCBI database under **BioProject accession number [PRJNA1262101](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1262101)**:

| Assembly Accession   | Organism Name                         | Assembly Level     |
|----------------------|----------------------------------------|--------------------|
| GCA_050500255.1      | *Microcystis aeruginosa* FBCC-A68     | Complete Genome    |
| GCA_050500245.1      | *Microcystis ichthyoblabe* FBCC-A1114 | Complete Genome    |
| GCA_050500235.1      | *Microcystis protocystis* FBCC-A270   | Complete Genome    |

---
