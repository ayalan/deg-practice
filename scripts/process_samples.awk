# This script will process the content of sample_info.txt and make a cleaner
# file (cell_line_mapping.txt) that, for each sample, has a tab-separated entry 
# with accession number, title, and SRX like:
# 
# GSM5561685	A549 2	SRX12032383
# GSM5561686	A549 3	SRX12032384
# GSM5561687	A549 4	SRX12032385
# 
# Ideally we'd also have SSR in here too, but it was not included with
# metadata provided in SOFT files.

# Syntax: Takes an input and output folder. 
# Use command like awk -f process_samples.awk sample_info.txt > cell_line_mapping.txt

# Start by defining the output format - fields will be separated by tabs
BEGIN {OFS="\t"} 

# When a line contains "!Sample_title", process it
/!Sample_title/ {
    # Remove the prefix "!Sample_title = " from the line
    gsub("!Sample_title = ", "");
    # Store the cell line name in the variable "title"
    title=$0;
} 

# When a line contains "!Sample_geo_accession", process it
/!Sample_geo_accession/ {
    # Remove the prefix "!Sample_geo_accession = " from the line
    gsub("!Sample_geo_accession = ", "");
    # Store the GSM accession number in the variable "gsm"
    gsm=$0;
} 

# When a line contains "!Sample_relation =" followed by "SRA:", process it
/!Sample_relation =.*SRA:/ {
    # Remove the entire URL prefix, leaving only the SRX identifier
    # The \\? escapes the ? character results is special in regex
    gsub("!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=", "");
    
    # Store the SRX accession in the variable "srx"
    srx=$0;
    
    # Output a line with the GSM accession, cell line name, and SRX accession
    # These will be separated by tabs due to the OFS setting above
    print gsm, title, srx;
}
