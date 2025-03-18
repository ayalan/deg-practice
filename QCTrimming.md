## (Supplemental) 3. Raw Data Processing

If the data is raw and hasn't gone through QC or trimming yet, follow these instructions.

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