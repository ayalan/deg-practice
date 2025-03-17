# Single-Cell RNA-Seq Analysis Guide: From FASTQ to DEGs

## 1. Introduction and Background

This guide will walk you through analyzing single-cell RNA sequencing (scRNA-seq) data from [GSE183590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183590), a dataset from the paper [Single-cell RNA sequencing for the identification of early-stage lung cancer biomarkers from circulating blood](https://www.nature.com/articles/s41525-021-00248-y).

You'll learn how to process FASTQ files, cluster cell types, visualize the data, and identify differentially expressed genes (DEGs).

We'll use Docker containers to manage the software environment, making it easier to run tools without installing them directly on your system.

### Prerequisites

- Basic knowledge of command line interface
- Docker installed on your macOS system
- The FASTQ files you've already downloaded and processed from SRA
- The matrix.txt metadata file

### Overview of the Workflow

1. Setup and Preparation
2. Raw Data Processing
3. Exploratory Data Analysis
4. Cell Clustering and Annotation 
5. Differential Expression Analysis
6. Saving and Exporting Results

Let's begin!

## 2. Setup and Preparation

### Project Directory Structure

Before we begin, let's establish a clear directory structure for your project. This will help keep your analysis organized and make it easier to track your progress.

```
/deg-practice/                     # Main project directory
├── data/                          # Contains all your FASTQ files and the reference genome
│   ├── SRR15740035.fastq.gz
│   ├── SRR15740036.fastq.gz
│   ├── ...
│   └── matrix.txt                 # Metadata file from the study
├── metadata/                      # GEO metadata files
│   ├── GSE183590_family.soft.gz   # SOFT format metadata
│   └── GSE183590_family.xml.gz    # MINiML format metadata
├── reference/                     # Reference genome files
│   ├── hg38.fa.gz                 # Reference genome in gzipped FASTA format
│   ├── annotation.gtf             # Annotation file for featureCount
│   └── (other reference files as needed)
├── fastqc_results/                # Quality control results
│   ├── SRR15740035_fastqc.html
│   ├── SRR15740035_fastqc.zip
│   └── ...
├── seurat_analysis.rds            # Saved Seurat object
├── cell_clusters.csv              # Cluster assignments
├── umap_coordinates.csv           # UMAP coordinates for visualization
├── top_markers.csv                # Top marker genes by cluster
└── all_degs.csv                   # All differentially expressed genes
```

This structure separates your data by function and keeps processed results separate from raw data. We'll create these directories as needed throughout the analysis.

### Setting Up a Docker Environment for scRNA-seq Analysis

Instead of installing various tools, let's create a reproducible Docker environment with all the necessary software for scRNA-seq analysis:

Make a Dockerfile with `touch Dockerfile` and insert the following:

```bash
# Use the latest version of Bioconductor
FROM bioconductor/bioconductor_docker:latest

# Set maintainer information
LABEL maintainer="Your Name <your.email@example.com>"
LABEL description="A practice environment for scRNA-seq analysis"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies in a single layer to reduce image size
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    hisat2 \
    subread \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libfftw3-dev \
    libgsl-dev \
    parallel \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages with specific versions for reproducibility
RUN R -e "options(repos = c(CRAN = 'https://cran.r-project.org')); \
    BiocManager::install(version = '3.15', ask = FALSE); \
    BiocManager::install(c( \
      'Seurat', \
      'SingleCellExperiment', \
      'scater', \
      'scran', \
      'SC3', \
      'MAST', \
      'clusterProfiler', \
      'limma', \
      'edgeR', \
      'org.Hs.eg.db', \
      'patchwork' \
    ), ask = FALSE)"

# Create directories for data and results
RUN mkdir -p /data /results /reference

# Set working directory
WORKDIR /data

# Command to run when the container starts
CMD ["R"]
```

Build the Docker image
```
docker build -t scrnaseq-analysis:1.0 .
```

Verify the image was built successfully
```
docker images | grep scrnaseq-analysis
```

Now, let's create a docker-compose.yml file to make running the analysis environment easier:

```bash
cat > docker-compose.yml << EOF
version: '3'
services:
  analysis:
    image: scrnaseq-analysis:1.0
    volumes:
      - ./data:/data
      - ./results:/results
      - ./reference:/reference
    working_dir: /data
    # Uncomment the next line if you want to use the container interactively
    # command: /bin/bash
EOF
```

To run your analysis environment using docker-compose:

```bash
# Start the container and enter an interactive R session
docker-compose run --rm analysis R

# Or for an interactive bash shell
docker-compose run --rm analysis bash
```

If you prefer direct docker commands:

```bash
# Run an interactive R session
docker run -it --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  -v $(pwd)/reference:/reference \
  scrnaseq-analysis:1.0 R

# Or for an interactive bash shell
docker run -it --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  -v $(pwd)/reference:/reference \
  scrnaseq-analysis:1.0 bash
```

### Downloading and Using GEO Metadata

The SOFT and MINiML files from GEO provide crucial metadata about your samples. Let's download these and use them to understand the experimental design:

```bash
# Navigate to your metadata directory
cd metadata

# Download the SOFT and MINiML files from GEO
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183590/soft/GSE183590_family.soft.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183590/miniml/GSE183590_family.xml.gz

# Uncompress them for easier viewing
gunzip GSE183590_family.soft.gz
gunzip GSE183590_family.xml.gz

# Take a look at the SOFT file to understand sample information
grep -A 20 "!Sample_title" GSE183590_family.soft | head -n 40
```

This helps in identifying which SRA accession (and thus FASTQ file) corresponds to which cell line and experimental condition, crucial for downstream analysis.

### Download and Organize Your FASTQ Files

```bash
# Download all runs associated with the project
prefetch --option-file SRXXXXXX

# Convert all downloaded SRA files to FASTQ
for sra_file in $(find . -name "*.sra"); do
  fasterq-dump $sra_file
done

# Compress the resulting FASTQ files
gzip *.fastq
```

### Extracting Cell Line Information from GEO Metadata

To integrate cell line information with your expression data, you need to create a mapping between SRA run accessions and their corresponding cell lines using the GEO metadata:

```bash
# Extract sample information from the SOFT file
cd metadata
grep -A 5 "!Sample_title" GSE183590_family.soft > sample_info.txt
grep "!Sample_geo_accession" GSE183590_family.soft >> sample_info.txt
grep "!Sample_relation" GSE183590_family.soft >> sample_info.txt

# Create a mapping file between SRA accessions and cell lines
awk 'BEGIN {OFS="\t"} 
     /!Sample_title/ {gsub("!Sample_title = ", ""); title=$0} 
     /!Sample_geo_accession/ {gsub("!Sample_geo_accession = ", ""); gsm=$0} 
     /!Sample_relation/ {gsub("!Sample_relation = SRA: ", ""); print gsm, title, $0}' sample_info.txt > cell_line_mapping.txt
```
This script:
- Sets the output field separator to a tab `(OFS="\t")`
- Looks for lines containing `!Sample_title` and saves the content to a variable
- Looks for lines containing `!Sample_geo_accession` and saves the content
- Looks for lines with `!Sample_relation` and then prints the collected data

View the mapping file. This mapping file connects each GSM accession to its corresponding cell line and SRA run. 
```
cat cell_line_mapping.txt
```

Now create a metadata file that you can use to annotate your Seurat object:

```bash
# Create a metadata file for Seurat integration
cd ..
awk 'BEGIN {print "sample_id\tcell_line"} 
     {split($3, sra, ";"); 
      for (i in sra) {
        gsub(" ", "", sra[i]); 
        if (sra[i] != "") print sra[i], $2
      }
     }' metadata/cell_line_mapping.txt > data/cell_line_metadata.txt
```

This script:
- Creates a header line with column names
- Takes the SRA accession numbers from column 3 (which might contain multiple accessions separated by semicolons)
- Splits them into individual accession numbers
- For each accession number, pairs it with the cell line name from column 2
- Saves these pairs to a new file that can be used by Seurat

The result is a clean, tabular metadata file that maps each SRA accession number to its corresponding cell line, which will be used in the analysis to annotate cells.Later, after creating your Seurat object, you'll integrate this cell line information.

### Understanding Cell Lines and Experimental Design

The dataset GSE183590 contains scRNA-seq data from 4 cell lines. To understand the experimental setup, examine the metadata and the matrix file:

```bash
# Run a container to examine the matrix.txt file (located in your data directory)
docker run -it --rm -v $(pwd)/data:/data scrnaseq-analysis:1.0 bash
head -n 20 /data/matrix.txt
```

This file typically contains information about cell barcodes, gene IDs, and counts.

## 3. Raw Data Processing

### Understanding FASTQ Files

FASTQ files contain sequence reads from high-throughput sequencing experiments, along with quality scores for each base. Each read in a FASTQ file is represented by 4 lines:

1. A header line starting with '@' that contains the sequence identifier
2. The actual nucleotide sequence
3. A line starting with '+' (sometimes with additional information)
4. Quality scores encoded in ASCII characters

Let's first examine one of your FASTQ files to understand its structure:

```bash
# Navigate to your project directory
cd /path/to/deg-practice

# Look at the first few lines of a FASTQ file
zcat data/SRR15740035.fastq.gz | head -n 12
```

### Quality Control with FastQC

Before proceeding with analysis, it's essential to check the quality of your sequencing data:

```bash
# Create a directory for FastQC results
mkdir -p fastqc_results

# Run FastQC on all FASTQ files
fastqc data/*.fastq.gz -o fastqc_results

# Optionally, generate a summary report with MultiQC
# Install MultiQC if not already installed
# pip install multiqc
multiqc fastqc_results -o fastqc_results
```

This will generate HTML reports for each FASTQ file and a summary report if you use MultiQC. Look for:
- Per base sequence quality
- Overrepresented sequences
- GC content
- Sequence duplication levels

### Quality Trimming with BBDuk

Based on the quality assessment, we'll perform adapter removal and quality trimming using BBDuk. This step improves downstream analysis by removing low-quality bases and adapter sequences.

```bash
# Create a directory for trimmed reads
mkdir -p data/trimmed

# Download BBMap (which includes BBDuk)
wget https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz
tar -xvzf BBMap_38.90.tar.gz

# Create a Docker container with BBDuk
docker run -it --rm \
  -v $(pwd)/data:/input \
  -v $(pwd)/data/trimmed:/output \
  -v $(pwd)/bbmap:/bbmap \
  ubuntu:20.04 bash

# Inside the container, run BBDuk on each file
for file in /input/*.fastq.gz; do
  basename=$(basename "$file" .fastq.gz)
  /bbmap/bbduk.sh \
    in=$file \
    out=/output/${basename}_trimmed.fastq.gz \
    ref=/bbmap/resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 \
    qtrim=r trimq=30 minlen=30 \
    threads=4
done
```

This command does the following:
- `ktrim=r`: Trims adapters from the right side (3' end)
- `k=23 mink=11 hdist=1`: Parameters for adapter matching
- `qtrim=r trimq=30`: Trims low-quality bases (Phred score Q<30) from the right side, equivalent to the 10^-3 error rate mentioned in the paper
- `minlen=30`: Discards reads shorter than 30 bases after trimming, matching the paper's specification
- `threads=4`: Uses 4 threads for processing

After trimming, you can run FastQC again on the trimmed reads to verify quality improvement:

```bash
fastqc data/trimmed/*_trimmed.fastq.gz -o fastqc_results/post_trimming
```

In subsequent steps, we'll use these trimmed FASTQ files instead of the original ones.

### Processing scRNA-seq Data with HISAT2 and featureCount

For data derived from the Fluidigm C1 platform, we will use HISAT2 for alignment and featureCount for generating the count matrix:

```bash
# Create necessary directories if they don't exist
mkdir -p data/aligned
mkdir -p reference
mkdir -p fastqc_results
mkdir -p metadata

# Prepare the reference genome
gunzip -c reference/hg38.fa.gz > reference/hg38.fa

# Download GTF annotation file for human genome (hg38/GRCh38) from Ensembl
wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
gunzip Homo_sapiens.GRCh38.104.gtf.gz
mv Homo_sapiens.GRCh38.104.gtf reference/annotation.gtf

# Create a HISAT2 index
hisat2-build reference/hg38.fa reference/hisat2_index

# Align each trimmed FASTQ file using HISAT2
for file in data/trimmed/*_trimmed.fastq.gz; do
  base=$(basename "$file" _trimmed.fastq.gz)
  hisat2 -p 4 -x reference/hisat2_index -U $file -S data/aligned/${base}.sam
done

# Convert SAM to BAM and sort using samtools
for file in data/aligned/*.sam; do
  base=$(basename "$file" .sam)
  samtools view -bS $file | samtools sort -o data/aligned/${base}.sorted.bam
done

# Run featureCount to generate count matrix
featureCounts -T 4 -a reference/annotation.gtf -o data/counts.txt data/aligned/*.sorted.bam
```

## 4. Exploratory Data Analysis

### Count Matrix Generation and Quality Control in R

Now, let's analyze the count matrix in R using the Seurat package:

```bash
# Run R in the container with all project directories mounted
docker run -it --rm \
    -v $(pwd)/data:/data \
    -v $(pwd):/project \
    scrnaseq-analysis:1.0 R
```

Within R:

```R
# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(scater)
library(scran)
library(SingleCellExperiment)

# Read the count matrix generated by featureCount
counts <- read.table("/data/counts.txt", header = TRUE, row.names = 1)

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = "GSE183590", min.cells = 3, min.features = 200)

# Add cell line information from the metadata file
cell_line_meta <- read.table("/data/cell_line_metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Map SRA accessions to cell lines
# This assumes your cell/sample names in the Seurat object contain the SRA run accession
# Adjust the matching logic if your naming convention is different
sample_ids <- colnames(seurat_obj)
cell_lines <- rep(NA, length(sample_ids))

for (i in 1:length(sample_ids)) {
  # Extract the SRA accession from the sample name
  sra_match <- regexpr("SRR[0-9]+", sample_ids[i])
  if (sra_match != -1) {
    sra_id <- substr(sample_ids[i], sra_match, sra_match + attr(sra_match, "match.length") - 1)
    # Find the cell line for this SRA accession
    match_idx <- which(cell_line_meta$sample_id == sra_id)
    if (length(match_idx) > 0) {
      cell_lines[i] <- cell_line_meta$cell_line[match_idx[1]]
    }
  }
}

# Add the cell line information to the Seurat object
seurat_obj$cell_line <- cell_lines

# Verify cell line assignment success
table(seurat_obj$cell_line, useNA = "ifany")

# Calculate mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC metrics
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
```

### Normalization and Dimensionality Reduction

After filtering, normalize the data and perform dimensionality reduction:

```R
# Normalize data
seurat_obj <- NormalizeData(seurat_obj)

# Identify highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

# Perform PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Visualize PCA results
ElbowPlot(seurat_obj)

# Perform t-SNE dimensionality reduction (as used in the original paper)
seurat_obj <- RunTSNE(
  seurat_obj,
  dims = 1:15,         # Use the first 15 PCs, adjust based on ElbowPlot
  perplexity = 30,     # Default perplexity parameter
  reduction.name = "tsne"
)

# Visualize t-SNE results
p1 <- DimPlot(seurat_obj, reduction = "tsne", pt.size = 1, label = TRUE) + 
  ggtitle("t-SNE Visualization of Clusters")

# Save the t-SNE plot
ggsave(filename = "/project/tsne_visualization.pdf", plot = p1, width = 10, height = 8)

# Color t-SNE by cell line instead of cluster (similar to Fig 1b in the paper)
# Using the cell_line metadata column we added earlier
p2 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "cell_line", pt.size = 1) + 
  ggtitle("t-SNE Visualization by Cell Line")

# Save the cell line t-SNE plot
ggsave(filename = "/project/tsne_by_cellline.pdf", plot = p2, width = 10, height = 8)
```

## 5. Cell Clustering and Annotation

### Clustering and Visualization

Cluster the cells and visualize them using both UMAP and t-SNE:

```R
# Determine the number of PCs to use based on ElbowPlot
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)

# Visualize clusters using UMAP
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP Visualization of Clusters")
print(p_umap)

# Save the UMAP plot
ggsave(filename = "/project/umap_visualization.pdf", plot = p_umap, width = 10, height = 8)

# Visualize t-SNE results (already generated in the previous step)
p_tsne <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE) + 
  ggtitle("t-SNE Visualization of Clusters")
print(p_tsne)

# Create a side-by-side comparison of both visualizations
combined_plot <- p_umap + p_tsne + plot_layout(ncol = 2)
print(combined_plot)
ggsave(filename = "/project/umap_tsne_comparison.pdf", plot = combined_plot, width = 16, height = 8)

# Try different resolutions to find optimal clustering
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)

# Show both visualizations with the new resolution
p_umap_03 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP (Resolution 0.3)")
p_tsne_03 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE) + 
  ggtitle("t-SNE (Resolution 0.3)")
combined_03 <- p_umap_03 + p_tsne_03 + plot_layout(ncol = 2)
print(combined_03)

seurat_obj <- FindClusters(seurat_obj, resolution = 0.7)

# Show both visualizations with the new resolution
p_umap_07 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP (Resolution 0.7)")
p_tsne_07 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE) + 
  ggtitle("t-SNE (Resolution 0.7)")
combined_07 <- p_umap_07 + p_tsne_07 + plot_layout(ncol = 2)
print(combined_07)

# Add a note to load patchwork if not already installed
# This can be run beforehand if needed: install.packages("patchwork")
library(patchwork)
```

### Cell Line Annotation Using Metadata

To identify what cell types your clusters represent, examine marker genes using R:

```R
# Identify marker genes using limma
library(limma)
library(edgeR)

# Convert count data to DGEList object
dge <- DGEList(counts = counts)

# Normalize data
dge <- calcNormFactors(dge)

# Design matrix for differential expression
design <- model.matrix(~0 + seurat_obj$seurat_clusters)
colnames(design) <- levels(seurat_obj$seurat_clusters)

# Voom transformation
v <- voom(dge, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design)

# Contrast matrix for comparisons
contrast.matrix <- makeContrasts(Cluster1vsCluster2 = Cluster1 - Cluster2, levels = design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DEGs
degs <- topTable(fit2, adjust = "BH", number = Inf)

# View top DEGs
head(degs, n = 20)

# Save DEGs to a file
write.csv(degs, "/project/all_degs.csv")

# Find marker genes for each cluster using Seurat's implementation
# This is different from the limma-based approach above and provides markers specific to each cluster
top_markers <- FindAllMarkers(seurat_obj, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25)

top_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top_markers

# View the top markers
print(head(top_markers, n = 20))

# Save top markers to a file (this means we don't need the write.csv in Section 7)
write.csv(top_markers, "/project/top_markers.csv", row.names = FALSE)
```

Use GEO metadata from the MINiML or SOFT file to further annotate clusters with their corresponding cell lines.

## 6. Differential Expression Analysis

### Identifying Differentially Expressed Genes (DEGs)

Following the methodology from the original paper, we'll set Cluster 1 as our reference/control cluster and identify DEGs by comparing it with other clusters. We'll apply the specific fold-change threshold of ≥|2| and FDR-corrected P-value < 0.05 in R:

```R
# Load required libraries for DEG analysis
library(limma)
library(edgeR)
library(biomaRt)
library(pheatmap)

# Assuming Cluster 1 corresponds to cluster label 0 in Seurat
# (adjust accordingly if your cluster numbers are different)
reference_cluster <- 1 # Used to be 0

# Create a list to store results from all comparisons
deg_results <- list()

# Set fold-change and p-value thresholds as specified in the paper
logfc_threshold <- 1  # log2(2) = 1
pval_threshold <- 0.05

# Perform pairwise comparisons using Cluster 1 as reference
for (cluster_id in unique(Idents(seurat_obj))) {
  if (cluster_id != reference_cluster) {
    # Calculate DEGs between reference cluster and current cluster
    comparison_name <- paste0("Cluster", reference_cluster, "_vs_Cluster", cluster_id)
    
    deg_results[[comparison_name]] <- FindMarkers(
      seurat_obj,
      ident.1 = reference_cluster,   # Reference cluster (Cluster 1)
      ident.2 = cluster_id,          # Target cluster
      logfc.threshold = logfc_threshold,
      min.pct = 0.1,                 # Detect genes that are expressed in at least 10% of cells
      test.use = "wilcox",           # Non-parametric Wilcoxon rank sum test
      only.pos = FALSE               # Get both up and down regulated genes
    )
    
    # Apply FDR correction and filtering
    deg_results[[comparison_name]]$FDR <- p.adjust(
      deg_results[[comparison_name]]$p_val, 
      method = "fdr"
    )
    
    # Filter DEGs by FDR and fold-change
    filtered_degs <- deg_results[[comparison_name]] %>%
      as.data.frame() %>%
      filter(FDR < pval_threshold & abs(avg_log2FC) >= logfc_threshold)
    
    # Add to the results list
    deg_results[[comparison_name]] <- filtered_degs
    
    # Print summary
    cat(paste0("DEGs in ", comparison_name, ": ", nrow(filtered_degs), 
               " (", sum(filtered_degs$avg_log2FC > 0), " up-regulated, ", 
               sum(filtered_degs$avg_log2FC < 0), " down-regulated)\n"))
  }
}

# Combine all DEGs into a single data frame
all_degs <- do.call(rbind, lapply(names(deg_results), function(comp_name) {
  if (nrow(deg_results[[comp_name]]) > 0) {
    df <- deg_results[[comp_name]]
    df$comparison <- comp_name
    df$gene <- rownames(df)
    return(df)
  } else {
    return(NULL)
  }
}))

# Save the combined DEGs to a file
write.csv(all_degs, "/project/all_degs.csv", row.names = FALSE)

# Cross-check gene IDs between Ensembl and NCBI gene databases
# Connect to Ensembl database
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Get unique gene symbols from DEGs
gene_symbols <- unique(rownames(all_degs))

# Map gene symbols to Ensembl and NCBI IDs
gene_id_mapping <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = ensembl
)

# Save gene ID mapping to a file
write.csv(gene_id_mapping, "/project/gene_id_mapping.csv", row.names = FALSE)

# Categorize DEGs into cluster-specific and cell-line-specific groups
# First, create a matrix of expression values for all DEGs across cells
deg_expr_matrix <- GetAssayData(seurat_obj, slot = "data")[rownames(all_degs), ]

# Add metadata for cell line information
cell_metadata <- data.frame(
  Cluster = Idents(seurat_obj),
  CellLine = seurat_obj$cell_line  # Using our mapped cell line information
)

# Generate heatmap for visual categorization
pheatmap(
  deg_expr_matrix,
  annotation_col = cell_metadata,
  show_rownames = FALSE,
  show_colnames = FALSE,
  scale = "row",
  clustering_method = "complete",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  filename = "/project/deg_heatmap.pdf",
  width = 10,
  height = 15
)

# Separate DEGs by up/down regulation per cluster comparison
up_c1_vs_c2 <- rownames(deg_results[["Cluster1_vs_Cluster2"]][deg_results[["Cluster1_vs_Cluster2"]]$avg_log2FC > 0, ])
down_c1_vs_c2 <- rownames(deg_results[["Cluster1_vs_Cluster2"]][deg_results[["Cluster1_vs_Cluster2"]]$avg_log2FC < 0, ])
up_c1_vs_c3 <- rownames(deg_results[["Cluster1_vs_Cluster3"]][deg_results[["Cluster1_vs_Cluster3"]]$avg_log2FC > 0, ])
down_c1_vs_c3 <- rownames(deg_results[["Cluster1_vs_Cluster3"]][deg_results[["Cluster1_vs_Cluster3"]]$avg_log2FC < 0, ])
up_c1_vs_c4 <- rownames(deg_results[["Cluster1_vs_Cluster4"]][deg_results[["Cluster1_vs_Cluster4"]]$avg_log2FC > 0, ])
down_c1_vs_c4 <- rownames(deg_results[["Cluster1_vs_Cluster4"]][deg_results[["Cluster1_vs_Cluster4"]]$avg_log2FC < 0, ])

# Find unique DEGs for each category (similar to Fig 2a in the paper)
unique_up_c1_vs_c2 <- setdiff(up_c1_vs_c2, c(up_c1_vs_c3, down_c1_vs_c3, up_c1_vs_c4, down_c1_vs_c4))
unique_down_c1_vs_c2 <- setdiff(down_c1_vs_c2, c(up_c1_vs_c3, down_c1_vs_c3, up_c1_vs_c4, down_c1_vs_c4))
unique_up_c1_vs_c3 <- setdiff(up_c1_vs_c3, c(up_c1_vs_c2, down_c1_vs_c2, up_c1_vs_c4, down_c1_vs_c4))
unique_down_c1_vs_c3 <- setdiff(down_c1_vs_c3, c(up_c1_vs_c2, down_c1_vs_c2, up_c1_vs_c4, down_c1_vs_c4))
unique_up_c1_vs_c4 <- setdiff(up_c1_vs_c4, c(up_c1_vs_c2, down_c1_vs_c2, up_c1_vs_c3, down_c1_vs_c3))
unique_down_c1_vs_c4 <- setdiff(down_c1_vs_c4, c(up_c1_vs_c2, down_c1_vs_c2, up_c1_vs_c3, down_c1_vs_c3))

# Create a table summarizing unique DEGs per category
unique_degs_summary <- data.frame(
  Comparison = c("Cluster1 vs Cluster2 (Up)", "Cluster1 vs Cluster2 (Down)",
                "Cluster1 vs Cluster3 (Up)", "Cluster1 vs Cluster3 (Down)",
                "Cluster1 vs Cluster4 (Up)", "Cluster1 vs Cluster4 (Down)"),
  Total_DEGs = c(length(up_c1_vs_c2), length(down_c1_vs_c2),
                length(up_c1_vs_c3), length(down_c1_vs_c3),
                length(up_c1_vs_c4), length(down_c1_vs_c4)),
  Unique_DEGs = c(length(unique_up_c1_vs_c2), length(unique_down_c1_vs_c2),
                 length(unique_up_c1_vs_c3), length(unique_down_c1_vs_c3),
                 length(unique_up_c1_vs_c4), length(unique_down_c1_vs_c4)),
  Percentage = c(length(unique_up_c1_vs_c2)/length(up_c1_vs_c2)*100,
                length(unique_down_c1_vs_c2)/length(down_c1_vs_c2)*100,
                length(unique_up_c1_vs_c3)/length(up_c1_vs_c3)*100,
                length(unique_down_c1_vs_c3)/length(down_c1_vs_c3)*100,
                length(unique_up_c1_vs_c4)/length(up_c1_vs_c4)*100,
                length(unique_down_c1_vs_c4)/length(down_c1_vs_c4)*100)
)

# Save the unique DEGs summary
write.csv(unique_degs_summary, "/project/unique_degs_summary.csv", row.names = FALSE)

# Generate volcano plots for each comparison
for (comparison in names(deg_results)[sapply(deg_results, nrow) > 0]) {
  # Create volcano plot data
  volcano_data <- deg_results[[comparison]]
  volcano_data$gene <- rownames(volcano_data)
  
  # Add significance and fold-change categories
  volcano_data$significance <- ifelse(volcano_data$FDR < pval_threshold, 
                                     ifelse(volcano_data$avg_log2FC > logfc_threshold, "Up-regulated",
                                           ifelse(volcano_data$avg_log2FC < -logfc_threshold, "Down-regulated", "Not Significant")),
                                     "Not Significant")
  
  # Plot
  p <- ggplot(volcano_data, aes(x = avg_log2FC, y = -log10(FDR), color = significance)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "gray")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", comparison),
         x = "Log2 Fold Change",
         y = "-Log10 FDR",
         color = "Regulation") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed")
  
  # Save plot
  ggsave(paste0("/project/volcano_", gsub(" ", "_", comparison), ".pdf"), plot = p, width = 10, height = 8)
}
```

This approach closely follows the DEG analysis strategy from the paper:
1. Uses Cluster 1 as the reference cluster for all comparisons
2. Applies the exact fold-change threshold (≥|2|) and FDR-corrected P-value (<0.05) from the paper
3. Cross-checks gene IDs between Ensembl and NCBI databases
4. Categorizes DEGs into cluster-specific and cell-line-specific groups
5. Identifies unique DEGs for each cluster comparison, similar to Figure 2a in the paper
6. Creates volcano plots to visualize the statistical significance of DEGs

The code will also generate several output files that help visualize and interpret the results:
- `all_degs.csv`: Combined list of all DEGs from all comparisons
- `gene_id_mapping.csv`: Mapping between gene symbols, Ensembl IDs, and NCBI IDs
- `deg_heatmap.pdf`: Heatmap for visual categorization of DEGs
- `unique_degs_summary.csv`: Summary of unique DEGs per category
- Volcano plots for each comparison

## 7. Saving and Exporting Results

Finally, save your Seurat object and export the results:

```R
# Save the Seurat object to the project directory
saveRDS(seurat_obj, file = "/project/seurat_analysis.rds")

# Export cluster assignments
clusters <- data.frame(Cell = names(Idents(seurat_obj)), Cluster = as.numeric(Idents(seurat_obj)))
write.csv(clusters, file = "/project/cell_clusters.csv", row.names = FALSE)

# Export UMAP coordinates
umap_coords <- data.frame(Cell = rownames(seurat_obj@reductions$umap@cell.embeddings),
                          UMAP1 = seurat_obj@reductions$umap@cell.embeddings[,1],
                          UMAP2 = seurat_obj@reductions$umap@cell.embeddings[,2])
write.csv(umap_coords, file = "/project/umap_coordinates.csv", row.names = FALSE)

# Note: top_markers have already been exported in Section 5
```

## 8. Conclusion

This guide has walked you through the essential steps for analyzing scRNA-seq data—from processing FASTQ files to identifying differentially expressed genes. With these modifications, the workflow now assumes the data is derived from the Fluidigm C1 platform and uses HISAT2 for alignment and featureCount for count matrix generation. Adjust parameters as needed based on your specific dataset.

## Additional Resources

- [Seurat vignettes](https://satijalab.org/seurat/vignettes.html)
- [Orchestrating Single-Cell Analysis with Bioconductor](https://osca.bioconductor.org/)
- [HISAT2 documentation](https://daehwankimlab.github.io/hisat2/)
- [Subread/featureCounts documentation](http://subread.sourceforge.net/)

## Troubleshooting Tips

- If you encounter memory issues, consider adjusting the Docker container's resource allocation.
- For large datasets, subsample the data for initial analyses.
- If clustering does not clearly separate cell types, try modifying the resolution parameter or number of principal components.
- Save intermediate steps to avoid losing progress.
