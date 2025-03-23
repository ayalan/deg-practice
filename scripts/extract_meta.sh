#!/bin/bash

# This script will extract only the lines we want from the SOFT file. It should
# output a file with content that has an entry of each sample, like :
#  
# !Sample_title = A549 2
# !Sample_geo_accession = GSM5561685
# !Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX12032383
#
# We actually want a simpler output for "SRA: blah blah" but we can handle it in
# the AWK script later.

# Start with a clean file
> output/sample_info.txt

# Process each sample
for gsm in $(grep '!Sample_geo_accession' GSE183590_family.soft | cut -d ' ' -f 3); do
  # Find the section for this specific GSM (using a simpler approach)
  section=$(grep -A 50 "!Sample_geo_accession = $gsm" GSE183590_family.soft | grep -m 50 "^!")
  
  # Extract title and GSM
  title=$(echo "$section" | grep '!Sample_title' | head -n 1)
  echo "$title" >> output/sample_info.txt
  echo "!Sample_geo_accession = $gsm" >> sample_info.txt
  
  # Extract SRA relation
  sra_line=$(echo "$section" | grep '!Sample_relation = SRA:' | head -n 1)
  if [[ -n "$sra_line" ]]; then
    echo "$sra_line" >> sample_info.txt
  else
    echo "!Sample_relation = SRA: Not_Found" >> sample_info.txt
  fi
  
  echo "" >> output/sample_info.txt  # Add an empty line between samples
done