# Single-Cell RNA-Seq Analysis Guide: From FASTQ to DEGs

## To Do's for this Guide

 - I overcomplicated this whole thing and started out with doing a lot of the clustering and analysis in seurat. Changed it to just use ASAP once count matrix and metadata is ready.
 - I want to add a 'Further Investigation' section about what else can be done with the data or other analysis that could be done.
 - The whole processing and visualization of up/down-regulation needs work.
 - I had significant help from Claude Sonnet 3.7 for the R portions of this document because it's not my forté. I want to properly review and clean it later.

## Introduction

We'll learn how to process FASTQ files, cluster cell types, visualize the data, and identify differentially expressed genes (DEGs) by analyzing scRNA-seq data from [GSE183590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183590), a dataset from the paper [Single-cell RNA sequencing for the identification of early-stage lung cancer biomarkers from circulating blood](https://www.nature.com/articles/s41525-021-00248-y). 

(I'm a newbie to bioinformatics and molecular biology so please pardon my mistakes. I'm updating this document as I go.)



We'll use Docker containers to manage the software environment, making it easier to run tools without installing them directly on our system.

### Prerequisites

- Basic knowledge of command line interface
- Docker should be installed (Dockerfiles are included in repo)

### Overview of Experiment Workflow

This diagram represents the experiments and observations completed in the paper as I understand them.

![Experiment Flow Chart](./experiment-flow.svg)

### What's Covered in this Exercise

This exercise aims to cover Steps 1 and 2 in the Experiment Workflow above.
Within those steps, we'll handle the following tasks:

1. Setup and Preparation
2. Raw Data Processing
3. Exploratory Data Analysis
4. Cell Clustering and Annotation 
5. Differential Expression Analysis
6. Saving and Exporting Results

Let's begin!

## 1. Setup and Preparation

### Project Directory Structure

```
/deg-practice/                     # Main project directory
├── data/                          # Contains our FASTQ files
│   ├── aligned/                   # Contains SAM and BAM files
│   │   └── ...                    
│   ├── SRR15740035.fastq.gz
│   ├── SRR15740036.fastq.gz
│   └── ...
├── metadata/                      # GEO metadata files
│   └── GSE183590_family.soft.gz   # SOFT format metadata
├── output/                        # Contains output of various scripts
│   └── ...
├── reference/                     # Contains reference genome files
│   ├── genome.fa                  # Reference genome in gzipped FASTA format
│   ├── annotation.gtf             # Annotation file for featureCount
│   └── ...
├── scripts/                       # Scripts to process data and metadata
│   ├── extract_meta.sh            # Extracts metadata we want from SOFT file
│   ├── map_srx_srr.sh             # Maps SRX to SRR metadata
│   ├── process_samples.awk        # Creates tabular file with SRX and cell line
│   └── sra-to-fastq.sh            # Converts all SRA files to FASTQ
├── fastqc_results/                # Quality control results (optional)
│   ├── SRR15740035_fastqc.html
│   ├── SRR15740035_fastqc.zip
│   └── ...
├── seurat_analysis.rds            # Saved Seurat object
├── cell_clusters.csv              # Cluster assignments
├── umap_coordinates.csv           # UMAP coordinates for visualization
├── top_markers.csv                # Top marker genes by cluster
└── all_degs.csv                   # All differentially expressed genes
```

This structure separates our data by function and keeps processed results separate from raw data. We'll create these directories as needed throughout the analysis.

### Setting Up a Docker Environment for scRNA-seq Analysis

Instead of installing various tools, let's create a reproducible Docker environment with all the necessary software for scRNA-seq analysis:

Check that we have the Dockerfiles, `Dockerfile` and `Dockerfile.eutils`.

Build the primary Docker image:
```
docker build -t scrnaseq-analysis:1.0 .

# or if using MacOS with mX chips:
docker build --platform linux/arm64 -t scrnaseq-analysis:1.0 .
```

Verify the image was built successfully
```
docker images | grep scrnaseq-analysis
```

Take a look at the docker-compose.yml file understand the analysis environment a little better.
To run our analysis environment using docker-compose:

```bash
# Start the container and enter an interactive R session
docker-compose run --rm analysis R

# Or for an interactive bash shell
docker-compose run --rm analysis bash
```

To exit or stop our Docker container when using docker-compose, we have a few options:

If we're in an interactive R session:
1. Type `q()` and press Enter
2. When prompted "Save workspace image?", type `n` (or our preference) and press Enter

If we're in an interactive bash shell:
1. Type exit and press Enter, or Press Ctrl+D

### Downloading and Using GEO Metadata

The SOFT and MINiML files from GEO provide crucial metadata about our samples. Let's download these and use them to understand the experimental design:

```bash
# Navigate to our metadata directory
cd metadata

# Download the SOFT and MINiML files from GEO
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183590/soft/GSE183590_family.soft.gz

# Uncompress them for easier viewing
gunzip GSE183590_family.soft.gz
```

Take a look at the SOFT file to understand sample information

```bash
grep -A 20 "!Sample_title" GSE183590_family.soft | head -n 40
```

If we're using zsh, use single quotes instead of double in the grep command.

This file helps in identifying which SRA accession (and thus FASTQ file) corresponds to which cell line and experimental condition, crucial for downstream analysis.

### Download and Organize Our FASTQ Files

```bash
# Use SRA Toolkit to all runs associated with the project from GEO
cd data
prefetch SRP336022

# Convert all downloaded SRA files to FASTQ
for sra_file in $(find . -name "*.sra"); do
  fasterq-dump $sra_file
done

# Or, do it with specified threads
fasterq-dump SRP336022 --threads 8

# Compress the resulting FASTQ files
gzip *.fastq
```

### Extracting Cell Line Information from GEO Metadata

To integrate cell line information with our expression data, we need to create a mapping between SRA run accessions and their corresponding cell lines using the GEO metadata.

Extract sample information from the SOFT file: 

```bash
# Make sure we're in the scripts folder
cd scripts

# Run script to extract metadata
./extract_meta.sh
```

Create a mapping file between SRA accessions and cell lines, using the process_samples.awk script:

```bash
awk -f process_samples.awk output/sample_info.txt > output/cell_line_mapping.txt
```

This script:
- Sets the output field separator to a tab `(OFS="\t")`
- Looks for lines containing `!Sample_title` and saves the content to a variable
- Looks for lines containing `!Sample_geo_accession` and saves the content
- Looks for lines with `!Sample_relation` and then prints the collected data

View the mapping file. This mapping file connects each GSM accession to its corresponding cell line and SRX. However, we're still missing the SRR, since the metadata files did not include them. We will need to query and map each SRX to SRR using the NCBI e-utils tool.

We'll do this with a one-time use Docker container. See Dockerfile.eutils in this repo.

```bash
docker build -t ncbi-eutils-py -f Dockerfile.eutils .
```

```bash
docker run -it --rm \
  -v $(pwd)/metadata:/data \
  ncbi-eutils-py \
  python3 /data/map_srx_srr.py
```

If you look at the newly made cell_line_metadata.txt, the cell line names contain a read number. We should remove those with this quick awk script:

```awk
awk 'BEGIN {OFS="\t"}
     NR==1 {print $0}  # Keep the header as-is
     NR>1 {
         # Print sample_id, then just the cell line name without the number
         split($2, parts, " ")
         print $1, parts[1]  # Print sample_id and just the first part of cell_line
     }' output/cell_line_metadata.txt > output/cell_line_metadata_asap.txt
```

The result is a clean, tabular metadata file that maps each SRA accession number to its corresponding cell line, which will be used in the analysis to annotate cells. 

To summarize: we started with the SOFT file, `used extract_meta.sh` to generate a `sample_info.txt` file. Then we used `process_samples.awk` to transform it to `cell_line_mapping.txt`. Then, we pulled the SSR for each SRX in it via NCBI using the `map_srx_srr.py` script to output `cell_line_metadata.txt`. We also generated `srx_to_srr_mapping.txt` which is mostly just for our own reference, and won't be used anywhere. We finally cleaned the file again and produced `cell_line_metadata_asap.txt` which we'll later use to import as metadata in ASAP.

## 2. Trimmed Data Processing

### Processing scRNA-seq Data with HISAT2 and featureCount

For data derived from the Fluidigm C1 platform, we will use HISAT2 for alignment and featureCount for generating the count matrix. (As I understand it, I might require other tools for aligning and counts if it were data from 10x.)

```bash
# Create necessary directories if they don't exist
mkdir -p data/aligned
mkdir -p reference

# Download the GRCh38.p13 reference genome with gget
cd reference
gget ref -w dna homo_sapiens -d
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > genome.fa

# Download the GRCh38.p13 gene annotations
gget ref -w dna homo_sapiens -d
gunzip -c Homo_sapiens.GRCh38.113.gtf.gz > annotation.gtf

# Create a HISAT2 index
time hisat2-build genome.fa hisat2_index

# Return to project root
cd ..

# Align each trimmed FASTQ file using HISAT2
for file in data/*.fastq.gz; do
  base=$(basename "$file" .fastq.gz)
  hisat2 -p 4 -x reference/hisat2_index -U $file -S data/aligned/${base}.sam
done

# Convert SAM to BAM and sort using samtools
for file in data/aligned/*.sam; do
  base=$(basename "$file" .sam)
  samtools view -bS $file | samtools sort -o data/aligned/${base}.sorted.bam
done

# Run featureCount to generate count matrix
featureCounts -T 4 -a reference/annotation.gtf -o output/counts.txt data/aligned/*.sorted.bam
```

The result is `counts.txt` and `counts.txt.summary`. To prepare the counts file for ASAP, we'll need to clean some of the file content. You'll notice that the first row is meant to be just a comment, and the second row uses the path name instead of sample name from column 7 and on. We can clean it using the following:

```awk
# Step 1: Skip the first line (metadata/command) and process from the second line
# Step 2: For the header row, keep "Geneid" and clean sample names (columns 7+)
# Step 3: Keep all data rows unchanged

awk 'BEGIN {OFS="\t"} 
     NR==2 {
         # Process the header row
         printf "Geneid";
         for(i=7; i<=NF; i++) {
             # Clean the sample name - extract just SRR ID
             sample = $i;
             gsub(".*/", "", sample);        # Remove path
             gsub("\\.sorted\\.bam", "", sample);  # Remove extension
             printf "\t%s", sample;
         }
         printf "\n";
     }
     NR>2 {
         # Process data rows - keep gene ID and count values
         printf "%s", $1;                    # Gene ID
         for(i=7; i<=NF; i++) {              # Count values
             printf "\t%s", $i;
         }
         printf "\n";
     }' output/counts.txt > output/counts_for_asap.txt
```

The processed file for import into ASAP is now ready.

## 3. Exploratory Data Analysis in ASAP

Upload our resulting `output/counts_for_asap.txt` file into the web-based ASAP platform.

Also upload the metadata that maps cell lines to same sample with `output/cell_line_metadata_asap.txt`

## 4. Cell Clustering and Annotation

TODO

### Hierarchical Clustering with SC3

Let's perform unsupervised single-cell consensus clustering using SC3, which implements complete-linkage hierarchical clustering.

### Visualizing SC3 Clustering Results

TODO

### Cell Line Annotation Using Metadata

Use GEO metadata from the SOFT file to further annotate clusters with their corresponding cell lines.

## 5. Differential Expression Analysis

TODO

### Identifying Differentially Expressed Genes (DEGs)

Following the methodology from the original paper, we'll set Cluster 1 as our reference/control cluster and identify DEGs by comparing it with other clusters. We'll apply the specific fold-change threshold of ≥\|2\| and FDR-corrected P-value < 0.05 in R:

## Conclusion

This guide has walked us through the essential steps for analyzing scRNA-seq data—from processing FASTQ files to identifying differentially expressed genes. Congratulations! 

This guide may be modified later to walk through some of the latter steps in the experiment's workflow. I will update it as I can.

## Further Investigation

- We could imputate with MAGIC, SAVER, or scImpute and look for any interesting changes in clustering or DEF results. 
- Try running analysis with UMAP instead of t-SNE with deep learning-based autoencoders (SCvi?) to find other cell subpopulations missed by SC3
- Try using LASSO or random forest importance scores to identify genes with the strongest discriminatory power between clusters or conditions

## Additional Resources

- [Seurat vignettes](https://satijalab.org/seurat/vignettes.html)
- [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/)
- [HISAT2 documentation](https://daehwankimlab.github.io/hisat2/)
- [Subread/featureCounts documentation](http://subread.sourceforge.net/)

## Troubleshooting Tips

- If we encounter memory issues, consider adjusting the Docker container's resource allocation.
- For large datasets, subsample the data for initial analyses.
- If clustering does not clearly separate cell types, try modifying the resolution parameter or number of principal components.
- Save intermediate steps to avoid losing progress.
