#!/usr/bin/env python3
"""

1. Unzips files matching the pattern source_<digits>.zip.
2. Reads the TSV file (matches_out_all.csv) from the unzipped directory.
3. Filters the rows based on the 18th column ('status (0.5)')—only retaining rows where
   the status is 'new_sh_in' or 'new_singleton_in'.
4. Extracts the sequence accession numbers (seq_accno, column 2) from these rows.
5. Opens the corresponding FASTA file from the indata directory (named identically to the source directory).
6. Extracts sequences whose headers match the filtered seq_accnos.
7. Appends the string ";sample=barcodeXX" (where XX is the numeric barcode extracted from the source file name) to each FASTA header.
8. Writes the filtered sequences into a new file (named source_<digits>_filtered) in a temporary directory.
9. Concatenates all filtered files into one file "concatenated_source" written in the current working directory.
10. Removes the intermediate extraction and filtered directories unless the flag --keep-temp is provided.

Requirements:
- Python 3
- Biopython (for FASTA parsing via Bio.SeqIO)
- argparse (standard with Python 3)
"""

import os
import glob
import zipfile
import csv
import re
import shutil
import argparse
from Bio import SeqIO

# Set up argument parser
parser = argparse.ArgumentParser(
    description="Filter sequences from Promethion runs and optionally clean up intermediate files."
)
parser.add_argument(
    "--keep-temp", 
    action="store_true", 
    help="Keep the intermediate extraction and filtered directories (default: remove them)."
)
args = parser.parse_args()

# Define base directories
OUTDATA_DIR = "path/to/outdata"
INDATA_DIR = "path/to/indata"
FILTERED_DIR = "filtered_sequences"

# Create output directory for filtered sequences if it doesn't exist
if not os.path.exists(FILTERED_DIR):
    os.makedirs(FILTERED_DIR)

# Keep track of extraction directories for later cleanup
extraction_dirs = []

# Process each zip file matching the pattern source_<digits>.zip
zip_files = glob.glob(os.path.join(OUTDATA_DIR, "source_[0-9]*.zip"))
for zip_file in zip_files:
    # Get the base file name and derive the file identifier (e.g., source_22 or source_123)
    base_name = os.path.basename(zip_file)      # e.g., source_22.zip
    file_id = os.path.splitext(base_name)[0]      # e.g., source_22

    # Extract the numeric barcode from the file_id using regex (matches any number of digits)
    match = re.search(r"source_(\d+)", file_id)
    if match:
        barcode = match.group(1)
    else:
        print(f"Barcode not found in filename {file_id}. Skipping.")
        continue

    # Create a temporary extraction directory for this zip file
    extract_dir = os.path.join(OUTDATA_DIR, file_id + "_extracted")
    extraction_dirs.append(extract_dir)
    if not os.path.exists(extract_dir):
        os.makedirs(extract_dir)

    # Unzip the file into the extraction directory
    with zipfile.ZipFile(zip_file, 'r') as zf:
        zf.extractall(path=extract_dir)

    # Build the correct path to the CSV file inside the extracted directory:
    # Expected path: <extract_dir>/matches/matches_out_all.csv
    csv_path = os.path.join(extract_dir, "matches", "matches_out_all.csv")
    if not os.path.exists(csv_path):
        print(f"CSV file not found: {csv_path}. Skipping this file.")
        continue

    # Read the TSV (tab-separated) file and filter for the desired statuses
    seq_accnos = set()
    with open(csv_path, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)  # Skip header row
        # Assuming columns: index 1 -> seq_accno and index 17 -> status (0.5)
        for row in reader:
            if len(row) < 18:
                continue
            status = row[17].strip()
            if status in ["new_sh_in", "new_singleton_in"]:
                seq_accno = row[1].strip()
                seq_accnos.add(seq_accno)

    # Open the corresponding FASTA file from the indata directory
    fasta_file = os.path.join(INDATA_DIR, file_id)
    if not os.path.exists(fasta_file):
        print(f"FASTA file not found: {fasta_file}. Skipping file {file_id}.")
        continue

    # Parse the FASTA file and filter sequences whose header (record.id) matches the seq_accnos
    filtered_sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in seq_accnos:
            # Append the barcode sample string to the header
            record.id = record.id + f";sample=barcode{barcode}"
            record.description = record.id  # Update description to match the new header
            filtered_sequences.append(record)

    # Write the filtered sequences to a new file named source_<digits>_filtered in FILTERED_DIR
    output_fasta = os.path.join(FILTERED_DIR, file_id + "_filtered")
    with open(output_fasta, "w") as out_f:
        SeqIO.write(filtered_sequences, out_f, "fasta")
    print(f"Processed {file_id}: {len(filtered_sequences)} sequences written to {output_fasta}")

# Concatenate all filtered FASTA files into a single file "concatenated_source" in the working directory
final_concatenated_file = os.path.join(os.getcwd(), "concatenated_source")
with open(final_concatenated_file, "w") as outfile:
    for filtered_fasta in glob.glob(os.path.join(FILTERED_DIR, "source_[0-9]*_filtered")):
        with open(filtered_fasta, "r") as infile:
            for line in infile:
                if line.startswith(">"):
                    # Replace only the first semicolon with an underscore
                    line = ">" + line[1:].replace(";", "_", 1)
                outfile.write(line)
print(f"All filtered sequences concatenated into {final_concatenated_file}")

# Clean up intermediate directories if the user has not set --keep-temp
if not args.keep_temp:
    # Remove extraction directories
    for ed in extraction_dirs:
        if os.path.exists(ed):
            shutil.rmtree(ed)
            print(f"Removed extraction directory: {ed}")
    # Remove the filtered sequences directory
    if os.path.exists(FILTERED_DIR):
        shutil.rmtree(FILTERED_DIR)
        print(f"Removed filtered sequences directory: {FILTERED_DIR}")
else:
    print("Intermediate files and directories retained (--keep-temp flag set).")
