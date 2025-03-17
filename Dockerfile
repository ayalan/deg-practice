# Use a specific version of Bioconductor for reproducibility
FROM bioconductor/bioconductor_docker:RELEASE_3_20

# Set maintainer information
LABEL maintainer="Aya Walraven <ayalan@gmail.com>"
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

# Install R packages with compatible version
RUN R -e "options(repos = c(CRAN = 'https://cran.r-project.org')); \
    BiocManager::install(ask = FALSE); \
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

