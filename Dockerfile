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
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for genomic data retrieval in a virtual environment
RUN python3 -m venv /opt/venv && \
    . /opt/venv/bin/activate && \
    pip3 install --no-cache-dir gget ffq && \
    # Make sure venv is activated when container starts
    echo '. /opt/venv/bin/activate' >> ~/.bashrc

# Add the virtual environment to PATH
ENV PATH="/opt/venv/bin:$PATH"

# Install SRA toolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    cp -r sratoolkit.*/bin/* /usr/local/bin/ && \
    rm -rf sratoolkit.* && \
    mkdir -p /root/.ncbi && \
    printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg && \
    printf '/repository/user/main/public/root = "%s"\n' "/data/sra-cache" >> /root/.ncbi/user-settings.mkfg && \
    printf '/repository/user/default-path = "%s"\n' "/data/sra-cache" >> /root/.ncbi/user-settings.mkfg

# Install R packages with compatible version
RUN R -e "options(repos = c(CRAN = 'https://cran.r-project.org')); \
    BiocManager::install(ask = FALSE); \
    BiocManager::install(c( \
      'Seurat', \
      'SingleCellExperiment', \
      'scater', \
      'SC3', \
      'clusterProfiler', \
      'limma' \
    ), ask = FALSE)"

# Create directories for data and results
RUN mkdir -p /data /results /reference /metadata /output /scripts

# Set working directory
WORKDIR /data

# Command to run when the container starts
CMD ["R"]
