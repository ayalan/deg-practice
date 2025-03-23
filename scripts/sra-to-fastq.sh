# Convert all downloaded SRA files to FASTQ
for sra_file in $(find data/ -name "*.sra"); do
  fasterq-dump $sra_file
done

# Compress the resulting FASTQ files
gzip data/*.fastq