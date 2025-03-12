# RNA-seq Analysis: Differential Gene Expression and Cell Clustering Using Docker

Follow this guide to analyze RNA-seq data from a RNA-seq GSE dataset to identify differentially expressed genes (DEGs) and cluster cells based on their expression profiles. We'll use Docker containers to keep your system clean while leveraging powerful bioinformatics tools.

## Table of Contents
- [Prerequisites](#prerequisites)
- [Project Setup](#project-setup)
- [Understanding Your Dataset](#understanding-your-dataset)
- [Quality Control and Preprocessing](#quality-control-and-preprocessing)
- [Quantifying Gene Expression](#quantifying-gene-expression)
- [Differential Gene Expression Analysis](#differential-gene-expression-analysis)
- [Cell Clustering](#cell-clustering)
- [Visualizing Results](#visualizing-results)
- [Complete Workflow Script](#complete-workflow-script)
- [Troubleshooting](#troubleshooting)

## Prerequisites

Before starting, make sure you have:
- Docker installed and running
- Your GSE dataset downloaded (FASTQ files)
- Basic familiarity with the command line
- Approximately 50GB of free disk space

## Project Setup

### 1. Create Project Structure

```bash
# Create main project directory
mkdir -p ~/deg-practice/{data,results,scripts,metadata}
cd ~/deg-practice

# Create subdirectories for results
mkdir -p results/{qc,counts,degs,clusters,visualization}
```

### 2. Organize Your FASTQ Files

```bash
# If your FASTQ files are spread across multiple directories, consolidate them
find /path/to/downloaded/gse/data -name "*.fastq.gz" -exec ln -s {} data/ \;
```

### 3. Create a Sample Metadata File

Create a file `metadata/samples.csv` with sample information:

```csv
sample_id,condition,replicate
SRR15739641,A549,1
SRR15739642,A549,2
# ... add all samples
```

## Understanding the Dataset

GSE183590 is a single-cell RNA-seq dataset from non-small cell lung cancer. Before proceeding, it's important to understand:

```bash
# Count your samples
ls data/*.fastq.gz | wc -l

# View file naming pattern
ls data/*.fastq.gz | head -5

# Check file contents (beginning of a FASTQ file)
zcat data/$(ls data/*.fastq.gz | head -1) | head -12
```

## Quality Control and Preprocessing

### 1. Run FastQC on All Samples

```bash
# Pull FastQC Docker image
docker pull biocontainers/fastqc:v0.11.9_cv8

# Run FastQC on all samples
docker run --rm -v $(pwd):/data -w /data biocontainers/fastqc:v0.11.9_cv8 \
  fastqc -t 4 -o results/qc data/*.fastq.gz
```

### 2. Generate QC Summary with MultiQC

```bash
# Pull MultiQC Docker image
docker pull ewels/multiqc:latest

# Run MultiQC to summarize FastQC results
docker run --rm -v $(pwd):/data -w /data ewels/multiqc:latest \
  multiqc results/qc -o results/qc
```

### 3. Trim Reads if Necessary

If FastQC reveals adapter contamination or poor quality bases:

```bash
# Pull Trimmomatic Docker image
docker pull quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2

# Run Trimmomatic on each sample
for fq in data/*.fastq.gz; do
  sample=$(basename $fq .fastq.gz)
  docker run --rm -v $(pwd):/data -w /data quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2 \
    trimmomatic SE \
    /data/$fq \
    /data/data/${sample}.trimmed.fastq.gz \
    ILLUMINACLIP:/data/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

## Quantifying Gene Expression

### 1. Pull Required Container

```bash
# Pull Salmon container for transcript quantification
docker pull combinelab/salmon:latest
```

### 2. Download Reference Transcriptome

```bash
# Create reference directory
mkdir -p reference

# Download human reference transcriptome
docker run --rm -v $(pwd):/data -w /data ubuntu:20.04 bash -c \
  "apt-get update && apt-get install -y wget && \
   wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz -O /data/reference/gencode.v38.transcripts.fa.gz"
```

### 3. Create Salmon Index

```bash
# Build Salmon index
docker run --rm -v $(pwd):/data -w /data combinelab/salmon:latest \
  salmon index -t /data/reference/gencode.v38.transcripts.fa.gz \
               -i /data/reference/salmon_index \
               -p 4
```

### 4. Quantify Transcript Expression

```bash
# Create directory for Salmon results
mkdir -p results/counts/salmon

# Run Salmon quantification for each sample
for fq in data/*.fastq.gz; do
  sample=$(basename $fq .fastq.gz)
  docker run --rm -v $(pwd):/data -w /data combinelab/salmon:latest \
    salmon quant -i /data/reference/salmon_index \
                -l A \
                -r /data/$fq \
                -p 4 \
                --validateMappings \
                -o /data/results/counts/salmon/$sample
done
```

## Differential Gene Expression Analysis

For this step, we'll use R with DESeq2 inside a Docker container.

### 1. Create an R Script for DEG Analysis

Create a file `scripts/deg_analysis.R` with:

```R
# Load libraries
library(tximport)
library(DESeq2)
library(dplyr)
library(ggplot2)

# Read sample metadata
samples <- read.csv("/data/metadata/samples.csv")

# Get Salmon quant files
files <- file.path("/data/results/counts/salmon", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

# Import transcript-level estimates
tx2gene <- read.csv("/data/reference/tx2gene.csv")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]  # Filter low-expression genes

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast=c("condition", "A549", "control"))
res <- res[order(res$padj), ]

# Export results
write.csv(as.data.frame(res), "/data/results/degs/deseq2_results.csv")

# Create MA plot
png("/data/results/degs/ma_plot.png", width=8, height=6, units="in", res=300)
plotMA(res, ylim=c(-5,5))
dev.off()

# Create PCA plot
vsd <- vst(dds, blind=FALSE)
png("/data/results/degs/pca_plot.png", width=8, height=6, units="in", res=300)
plotPCA(vsd, intgroup="condition")
dev.off()

# Create heatmap of top DEGs
library(pheatmap)
top_genes <- rownames(res)[1:50]
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
png("/data/results/degs/heatmap.png", width=10, height=8, units="in", res=300)
pheatmap(mat, annotation_col=samples[,c("condition"),drop=FALSE])
dev.off()

# Create volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")

png("/data/results/degs/volcano_plot.png", width=8, height=6, units="in", res=300)
ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("No"="gray", "Yes"="red")) +
  theme_bw() +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted P-value")
dev.off()

# Summary statistics
sig_genes <- sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm=TRUE)
up_genes <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm=TRUE)
down_genes <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm=TRUE)

write(paste("Significant DEGs:", sig_genes), "/data/results/degs/summary.txt")
write(paste("Upregulated:", up_genes), "/data/results/degs/summary.txt", append=TRUE)
write(paste("Downregulated:", down_genes), "/data/results/degs/summary.txt", append=TRUE)
```

### 2. Prepare tx2gene File

Create a transcript-to-gene mapping file:

```bash
# Create tx2gene file
docker run --rm -v $(pwd):/data -w /data --entrypoint=bash rocker/tidyverse:4.1.0 -c \
  "Rscript -e \"library(biomaRt); \
  mart <- useMart('ensembl', dataset='hsapiens_gene_ensembl'); \
  tx <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name'), mart=mart); \
  write.csv(tx, '/data/reference/tx2gene.csv', row.names=FALSE)\""
```

### 3. Run DEG Analysis

```bash
# Pull R container with Bioconductor
docker pull rocker/tidyverse:4.1.0

# Install required R packages
docker run --rm -v $(pwd):/data -w /data --entrypoint=bash rocker/tidyverse:4.1.0 -c \
  "R -e \"install.packages('BiocManager'); \
  BiocManager::install(c('tximport', 'DESeq2', 'pheatmap'))\""

# Run DEG analysis
docker run --rm -v $(pwd):/data -w /data rocker/tidyverse:4.1.0 \
  Rscript /data/scripts/deg_analysis.R
```

## Cell Clustering

For cell clustering based on gene expression patterns, we'll use Seurat in R:

### 1. Create an R Script for Clustering

Create a file `scripts/cell_clustering.R` with:

```R
# Load libraries
library(Seurat)
library(tidyverse)
library(ggplot2)

# Read expression data
# First, create a matrix from Salmon counts
samples <- read.csv("/data/metadata/samples.csv")
files <- file.path("/data/results/counts/salmon", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

# Function to convert Salmon counts to a matrix
get_counts <- function(files) {
  all_counts <- list()
  for (i in seq_along(files)) {
    quant <- read.delim(files[i], sep="\t")
    all_counts[[names(files)[i]]] <- quant$NumReads
    rownames(all_counts[[names(files)[i]]]) <- quant$Name
  }
  counts_matrix <- do.call(cbind, all_counts)
  return(counts_matrix)
}

counts <- get_counts(files)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = "GSE183590", min.cells = 3, min.features = 200)

# Add metadata
seurat_obj$condition <- samples$condition[match(colnames(seurat_obj), samples$sample_id)]

# Normalize data
seurat_obj <- NormalizeData(seurat_obj)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Determine dimensionality
png("/data/results/clusters/elbow_plot.png", width=8, height=6, units="in", res=300)
ElbowPlot(seurat_obj, ndims = 50)
dev.off()

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Run UMAP and t-SNE
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20)

# Save the Seurat object
saveRDS(seurat_obj, "/data/results/clusters/seurat_object.rds")

# Generate cluster plots
# UMAP colored by cluster
png("/data/results/clusters/umap_clusters.png", width=10, height=8, units="in", res=300)
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")
dev.off()

# UMAP colored by condition
png("/data/results/clusters/umap_condition.png", width=10, height=8, units="in", res=300)
DimPlot(seurat_obj, reduction = "umap", group.by = "condition")
dev.off()

# t-SNE colored by cluster
png("/data/results/clusters/tsne_clusters.png", width=10, height=8, units="in", res=300)
DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters")
dev.off()

# Find markers for each cluster
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, "/data/results/clusters/cluster_markers.csv", row.names = FALSE)

# Generate heatmap of top markers
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

png("/data/results/clusters/markers_heatmap.png", width=12, height=10, units="in", res=300)
DoHeatmap(seurat_obj, features = top_markers$gene) + NoLegend()
dev.off()

# Generate feature plots for top markers
top_genes <- top_markers %>%
  top_n(n = 9, wt = avg_log2FC) %>%
  pull(gene)

png("/data/results/clusters/feature_plot.png", width=12, height=10, units="in", res=300)
FeaturePlot(seurat_obj, features = top_genes, ncol = 3)
dev.off()

# Generate violin plots for top markers
png("/data/results/clusters/violin_plot.png", width=12, height=10, units="in", res=300)
VlnPlot(seurat_obj, features = top_genes[1:6], ncol = 3)
dev.off()
```

### 2. Run Cell Clustering

```bash
# Install required R packages for clustering
docker run --rm -v $(pwd):/data -w /data --entrypoint=bash rocker/tidyverse:4.1.0 -c \
  "R -e \"install.packages('BiocManager'); \
  BiocManager::install('Seurat')\""

# Run clustering analysis
docker run --rm -v $(pwd):/data -w /data rocker/tidyverse:4.1.0 \
  Rscript /data/scripts/cell_clustering.R
```

## Visualizing Results

After running the analyses, you'll have several visualizations:

### Differential Expression Results
- MA Plot: `results/degs/ma_plot.png`
- PCA Plot: `results/degs/pca_plot.png`
- Heatmap: `results/degs/heatmap.png`
- Volcano Plot: `results/degs/volcano_plot.png`

### Clustering Results
- UMAP Plots: `results/clusters/umap_clusters.png` and `results/clusters/umap_condition.png`
- t-SNE Plot: `results/clusters/tsne_clusters.png`
- Marker Heatmap: `results/clusters/markers_heatmap.png`
- Feature Plots: `results/clusters/feature_plot.png`
- Violin Plots: `results/clusters/violin_plot.png`

## Complete Workflow Script

Create `scripts/run_analysis.sh` to automate the workflow:

```bash
#!/bin/bash
# RNA-seq Analysis Workflow

set -e  # Exit on error

# Create necessary directories
mkdir -p ~/deg-practice/{data,results,scripts,metadata,reference}
cd ~/deg-practice
mkdir -p results/{qc,counts,degs,clusters,visualization}

# Step 1: Quality Control
echo "=== Running FastQC ==="
docker pull biocontainers/fastqc:v0.11.9_cv8
docker run --rm -v $(pwd):/data -w /data biocontainers/fastqc:v0.11.9_cv8 \
  fastqc -t 4 -o results/qc data/*.fastq.gz

echo "=== Generating QC Summary ==="
docker pull ewels/multiqc:latest
docker run --rm -v $(pwd):/data -w /data ewels/multiqc:latest \
  multiqc results/qc -o results/qc

# Step 2: Prepare Reference
echo "=== Downloading Reference Transcriptome ==="
docker run --rm -v $(pwd):/data -w /data ubuntu:20.04 bash -c \
  "apt-get update && apt-get install -y wget && \
   wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz -O /data/reference/gencode.v38.transcripts.fa.gz"

echo "=== Creating Salmon Index ==="
docker pull combinelab/salmon:latest
docker run --rm -v $(pwd):/data -w /data combinelab/salmon:latest \
  salmon index -t /data/reference/gencode.v38.transcripts.fa.gz \
               -i /data/reference/salmon_index \
               -p 4

# Step 3: Quantify Gene Expression
echo "=== Quantifying Gene Expression ==="
mkdir -p results/counts/salmon
for fq in data/*.fastq.gz; do
  sample=$(basename $fq .fastq.gz)
  echo "Processing $sample"
  docker run --rm -v $(pwd):/data -w /data combinelab/salmon:latest \
    salmon quant -i /data/reference/salmon_index \
                -l A \
                -r /data/$fq \
                -p 4 \
                --validateMappings \
                -o /data/results/counts/salmon/$sample
done

# Step 4: Create tx2gene mapping
echo "=== Creating Transcript-to-Gene Mapping ==="
docker pull rocker/tidyverse:4.1.0
docker run --rm -v $(pwd):/data -w /data --entrypoint=bash rocker/tidyverse:4.1.0 -c \
  "Rscript -e \"if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
  if (!requireNamespace('biomaRt', quietly = TRUE)) BiocManager::install('biomaRt'); \
  library(biomaRt); \
  mart <- useMart('ensembl', dataset='hsapiens_gene_ensembl'); \
  tx <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name'), mart=mart); \
  write.csv(tx, '/data/reference/tx2gene.csv', row.names=FALSE)\""

# Step 5: Install required R packages
echo "=== Installing R Packages ==="
docker run --rm -v $(pwd):/data -w /data --entrypoint=bash rocker/tidyverse:4.1.0 -c \
  "R -e \"if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); \
  BiocManager::install(c('tximport', 'DESeq2', 'pheatmap', 'Seurat'))\""

# Step 6: Differential Gene Expression Analysis
echo "=== Running Differential Expression Analysis ==="
docker run --rm -v $(pwd):/data -w /data rocker/tidyverse:4.1.0 \
  Rscript /data/scripts/deg_analysis.R

# Step 7: Cell Clustering
echo "=== Running Cell Clustering ==="
docker run --rm -v $(pwd):/data -w /data rocker/tidyverse:4.1.0 \
  Rscript /data/scripts/cell_clustering.R

echo "=== Analysis Complete! ==="
echo "Results are available in the 'results' directory."
```

Make the script executable and run it:

```bash
chmod +x scripts/run_analysis.sh
./scripts/run_analysis.sh
```

## Troubleshooting

### Common Issues and Solutions:

1. **Docker permission issues**: 
   If you get permission errors, try running the Docker commands with `sudo` or add your user to the Docker group.

2. **Memory limitations**:
   - For large datasets, increase Docker's memory allocation in Docker Desktop settings
   - For Salmon and R analyses, add memory parameters (e.g., `--max-memory=16G` for Salmon)

3. **Missing packages in R**:
   If the workflow stops due to missing R packages, you can install them separately:
   ```bash
   docker run --rm -v $(pwd):/data -w /data --entrypoint=bash rocker/tidyverse:4.1.0 -c \
     "R -e \"BiocManager::install('packageName')\""
   ```

4. **Unmapped reads in Salmon**:
   If Salmon reports a high percentage of unmapped reads, check:
   - That you're using the correct reference transcriptome for your species
   - That your library type is correctly specified (-l parameter)
   - For contamination in your samples

5. **Error in DESeq2/Seurat**:
   - Check that your metadata file correctly matches your sample IDs
   - Ensure the columns and format match what the scripts expect

### Getting Help:

If you encounter issues not covered here, these resources can help:
- Docker documentation: https://docs.docker.com/
- Salmon documentation: https://salmon.readthedocs.io/
- DESeq2 vignette: http://bioconductor.org/packages/DESeq2
- Seurat vignette: https://satijalab.org/seurat/articles/pbmc3k_tutorial

---

This guide provides a comprehensive workflow to analyze RNA-seq data, identify differentially expressed genes, and perform cell clustering using Docker. The workflow is scalable to handle your GSE dataset with many FASTQ files while keeping your system clean by containing all tools within Docker.