#!/usr/bin/env python3
import os
import time
import subprocess
import sys

# Set your NCBI API key as environment variable
os.environ["NCBI_API_KEY"] = "3a7625b342ed39e26b00676476f453ce4e08"
print("Starting SRX to SRR mapping with API key set in environment")

# Create output files
with open("srx_to_srr_mapping.txt", "w") as f:
    f.write("srx_id\tsrr_id\n")

with open("cell_line_metadata.txt", "w") as f:
    f.write("sample_id\tcell_line\n")

# Extract SRX IDs
srx_ids = []
with open("cell_line_mapping.txt", "r") as f:
    for line in f:
        if line.startswith("sample_id"):
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            srx_ids.append(parts[2])

print(f"Found {len(srx_ids)} SRX IDs to process")

# Read cell line mapping
srx_to_cell = {}
with open("cell_line_mapping.txt", "r") as f:
    for line in f:
        if line.startswith("sample_id"):
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            srx_to_cell[parts[2]] = parts[1]

# Process each SRX ID
count = 0
for srx in srx_ids:
    count += 1
    print(f"[{count}/{len(srx_ids)}] Processing {srx}")
    
    try:
        # Basic query without API key flag (using env var)
        cmd = f'esearch -db sra -query "{srx}"'
        print(f"Running: {cmd}")
        search_result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        print(f"Search result exit code: {search_result.returncode}")
        if search_result.returncode != 0:
            print(f"Error in search: {search_result.stderr}")
            continue
            
        # Extract count to confirm we found something
        if "<Count>0</Count>" in search_result.stdout:
            print(f"No results found for {srx}")
            continue
            
        # Now fetch the run info
        fetch_cmd = f'echo "{search_result.stdout}" | efetch -format runinfo'
        print(f"Running: {fetch_cmd}")
        fetch_result = subprocess.run(fetch_cmd, shell=True, capture_output=True, text=True)
        
        if fetch_result.returncode != 0:
            print(f"Error in fetch: {fetch_result.stderr}")
            continue
            
        # Process the result to extract SRR
        lines = fetch_result.stdout.strip().split("\n")
        if len(lines) > 1:
            srr = lines[1].split(",")[0]
            if srr.startswith("SRR"):
                print(f"{srx} -> {srr}")
                
                # Save to mapping file
                with open("srx_to_srr_mapping.txt", "a") as f:
                    f.write(f"{srx}\t{srr}\n")
                
                # Save to metadata file
                if srx in srx_to_cell:
                    with open("cell_line_metadata.txt", "a") as f:
                        f.write(f"{srr}\t{srx_to_cell[srx]}\n")
            else:
                print(f"{srx} -> No valid SRR found in result")
        else:
            print(f"{srx} -> No data lines in result")
            
    except Exception as e:
        print(f"Error processing {srx}: {e}")
    
    # Be nice to NCBI servers
    time.sleep(1.0)

print(f"Processed {count}/{len(srx_ids)} SRX IDs.")
print("Results saved to srx_to_srr_mapping.txt and cell_line_metadata.txt")