#!/bin/bash
set -euo pipefail

# -------------------------------------------------------------------
# Configuration (manually set variables)
# -------------------------------------------------------------------
THREADS=63                   # Number of threads to use per task
INPUT_FILE="path/meta_SH/barcode_230P1.fasta"   # Input FASTA file
MAX_SEQUENCES=50000          # Maximum number of sequences per split file
MODE="sequential"            # Execution mode: "parallel" or "sequential"
NODES=1                      # Number of nodes (only used in parallel mode)
SH_MATCHING_DIR="path/meta_SH/sh_matching_pub"  # Directory for SH-matching files and tools
OUTPUT_DIR="path/meta_SH/first_tryout"  # Directory for output files

# New parameter: choose whether to re-run unmatched sequences ("yes" or "no")
RERUN_UNMATCHED="yes"

# -------------------------------------------------------------------
# Check and Clear OUTPUT_DIR if it exists
# -------------------------------------------------------------------
if [[ -d "$OUTPUT_DIR" ]]; then
  echo "WARNING: Output directory $OUTPUT_DIR already exists. Clearing its contents..."
  rm -rf "$OUTPUT_DIR"/*
fi
mkdir -p "$OUTPUT_DIR"

# -------------------------------------------------------------------
# Logging: Create a timestamp-based log file in the OUTPUT_DIR
# -------------------------------------------------------------------
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
INPUT_BASENAME=$(basename "$INPUT_FILE")
LOG_FILE="$OUTPUT_DIR/job_${TIMESTAMP}_${INPUT_BASENAME}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date)] Starting processing."
echo "Input file: $INPUT_FILE"
echo "Max sequences per file: $MAX_SEQUENCES"
echo "Mode: $MODE"
echo "Threads: $THREADS"
if [[ "$MODE" == "parallel" ]]; then
  echo "Nodes: $NODES"
fi
echo "SH_MATCHING_DIR: $SH_MATCHING_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "RERUN_UNMATCHED: $RERUN_UNMATCHED"

# -------------------------------------------------------------------
# (1) Directory Checks and Setup
# -------------------------------------------------------------------
if [[ ! -d "$SH_MATCHING_DIR" ]]; then
  echo "Fatal error: SH_MATCHING_DIR ($SH_MATCHING_DIR) does not exist."
  exit 1
fi

required_subdirs=(indata outdata userdir sh_matching_analysis data_udb)
for sub in "${required_subdirs[@]}"; do
  if [[ ! -d "$SH_MATCHING_DIR/$sub" ]]; then
    echo "Fatal error: Required subdirectory '$sub' not found in $SH_MATCHING_DIR."
    exit 1
  fi
done

# -------------------------------------------------------------------
# (2) Sequence Counting and FASTA File Subdivision
# -------------------------------------------------------------------
TOTAL_SEQS=$(grep -c '^>' "$INPUT_FILE")
echo "Total sequences in FASTA: $TOTAL_SEQS"

NUM_FILES=$(( (TOTAL_SEQS + MAX_SEQUENCES - 1) / MAX_SEQUENCES ))
echo "Will split into $NUM_FILES file(s)."

if (( NUM_FILES > 999 )); then
  echo "Fatal error: Attempting to split into more than 999 files ($NUM_FILES)."
  exit 1
elif (( NUM_FILES > 50 )); then
  echo "Warning: Splitting into $NUM_FILES files. Are you certain about your max sequence number?"
fi

# Use a file-wide search for "sample=" and print the first five headers for debugging.
SAMPLE_FLAG=0
if grep -m 1 "sample=" "$INPUT_FILE" > /dev/null; then
  SAMPLE_FLAG=1
  echo "Detected 'sample=' in the input file; sample info will be mapped."
else
  echo "No 'sample=' string detected in the input file."
fi

TMP_SPLIT_DIR=$(mktemp -d)
MAPPING_FILE="$OUTPUT_DIR/mapping_${TIMESTAMP}.txt"
# New header includes MATCH_STATUS_005 column.
echo -e "INPUT_FILE\tSOURCE_FILE\tHEADER\tSAMPLE_INFO\tSOURCE_NAME\tRUN_TYPE\tMATCH_STATUS_005" > "$MAPPING_FILE"

awk -v max_seq="$MAX_SEQUENCES" \
    -v outdir="$TMP_SPLIT_DIR" \
    -v prefix="source_" \
    -v mapping_file="$MAPPING_FILE" \
    -v input_file="$INPUT_FILE" \
    -v sample_flag="$SAMPLE_FLAG" '
BEGIN {
    seqCount = 0;
    fileCount = 0;
}
# Process header lines
/^>/ {
    seqCount++;
    if ((seqCount - 1) % max_seq == 0) {
        fileCount++;
        id = sprintf("%03d", fileCount);
        current_file = outdir "/" prefix id;
    }
    print $0 > current_file;
    sample = "";
    if (sample_flag == 1) {
        n = split($0, fields, ";")
        for (i = 1; i <= n; i++) {
            if (fields[i] ~ /^sample=/) {
                sub(/^sample=/, "", fields[i])
                sample = fields[i]
                break
            }
        }
    }
    # Derive source name from current_file (basename)
    srcname = current_file; sub(".*/", "", srcname);
    # Print mapping: the header, the extracted sample info, etc.
    print input_file "\t" current_file "\t" $0 "\t" sample "\t" srcname "\tinitial_run\tunmatched" >> mapping_file;
    next;
}
{
    print $0 >> current_file;
}
END {
    print "Created " fileCount " split files." > "/dev/stderr";
}' "$INPUT_FILE"

echo "Moving split files to $SH_MATCHING_DIR/indata..."
echo "from $TMP_SPLIT_DIR"
mv "$TMP_SPLIT_DIR"/source_* "$SH_MATCHING_DIR/indata/"
rm -r "$TMP_SPLIT_DIR"  # Remove the temporary split directory

# -------------------------------------------------------------------
# Prepare List of Unique Source Files (from mapping) to Process
# -------------------------------------------------------------------
# Store unique IDs in a Bash array using AWK's associative array (avoids sort)
mapfile -t source_ids < <(awk 'NR>1 {
    id = $2;
    sub(".*/", "", id);
    if (!seen[id]++) print id;
}' "$MAPPING_FILE")
echo "Unique source IDs to process: ${source_ids[*]}"

# -------------------------------------------------------------------
# (5) Running SH-Matching Pipeline (Initial Run)
# -------------------------------------------------------------------
if [[ "$MODE" == "sequential" ]]; then
  echo "Running in sequential mode..."
  pushd "$SH_MATCHING_DIR" > /dev/null
  for id in "${source_ids[@]}"; do
    num_id="${id#source_}"
    echo "Processing $id sequentially with run id $num_id..."
    ./sh_matching_echo.sif /sh_matching/run_pipeline.sh "$num_id" itsfull no yes no no
  done
  popd > /dev/null
else
  echo "Parallel mode not handled in this version."
fi

# -------------------------------------------------------------------
# (6) Output Aggregation: Unzip and Concatenate as TSV (Initial Run)
# -------------------------------------------------------------------
echo "Concatenating output files..."
TMP_OUT=$(mktemp -d)

for zipfile in "$SH_MATCHING_DIR"/outdata/source_*.zip; do
  if [[ -f "$zipfile" ]]; then
    base_zip=$(basename "$zipfile" .zip)
    target_dir="$TMP_OUT/$base_zip"
    echo "Unzipping $zipfile into $target_dir..."
    mkdir -p "$target_dir"
    unzip -o "$zipfile" -d "$target_dir"
  fi
done

base=$(basename "$INPUT_FILE")
output_file="$OUTPUT_DIR/${base%.*}_all_matches_out.tsv"
echo "Creating concatenated output file: $output_file"
> "$output_file"

for extracted_dir in "$TMP_OUT"/*; do
  echo "Checking directory: $(basename "$extracted_dir")"
  if [[ -d "$extracted_dir/matches" ]]; then
    matches_file="$extracted_dir/matches/matches_out_all.csv"
    if [[ -f "$matches_file" ]]; then
      echo "Found $matches_file"
      if [[ ! -s "$output_file" ]]; then
        echo "Appending full file from $matches_file"
        cat "$matches_file" >> "$output_file"
      else
        echo "Appending file from $matches_file (skipping header)"
        tail -n +2 "$matches_file" >> "$output_file"
      fi
    else
      echo "No matches_out_all.csv found in $(basename "$extracted_dir")"
    fi
  else
    echo "Directory $(basename "$extracted_dir") does not contain a matches folder, skipping."
  fi
done

echo "Concatenation complete. Final TSV output in $output_file"
rm -rf "$TMP_OUT"  # Remove the temporary output directory

# -------------------------------------------------------------------
# (6.5) Re-run Pipeline for Unmatched Sequences (if enabled)
# -------------------------------------------------------------------
if [[ "$RERUN_UNMATCHED" == "yes" ]]; then
  echo "Filtering unmatched sequences from original TSV..."
  unmatched_tsv="${OUTPUT_DIR}/${base%.*}_unmatched.tsv"
  matched_tsv="${OUTPUT_DIR}/${base%.*}_matched.tsv"
  new_total_tsv="${OUTPUT_DIR}/${base%.*}_new_total.tsv"
  
  # Filter unmatched rows (assuming column 18 holds the status)
  awk -F"\t" 'NR==1 || ($18=="new_sh_in" || $18=="new_singleton_in")' "$output_file" > "$unmatched_tsv"
  awk -F"\t" 'NR==1 || ($18!="new_sh_in" && $18!="new_singleton_in")' "$output_file" > "$matched_tsv"
  
  echo "Extracting sequence accession numbers from unmatched TSV..."
  unmatched_ids_file=$(mktemp)
  tail -n +2 "$unmatched_tsv" | awk -F"\t" '{print $2}' | sort | uniq > "$unmatched_ids_file"
  
  echo "Extracting unmatched sequences from the input FASTA..."
  rerun_fasta="${OUTPUT_DIR}/${base%.*}_rerun_input.fasta"
  python3 - <<EOF
import sys
from Bio import SeqIO
rerun_fasta = "$rerun_fasta"
unmatched_ids = set(line.strip() for line in open("$unmatched_ids_file"))
count = 0
with open("$INPUT_FILE") as infile, open(rerun_fasta, "w") as outfile:
    for rec in SeqIO.parse(infile, "fasta"):
        if rec.id in unmatched_ids:
            rec.id = rec.id  # Optionally modify header if needed.
            rec.description = rec.id
            count += 1
            SeqIO.write(rec, outfile, "fasta")
print(f"Extracted {count} unmatched sequences to {rerun_fasta}")
EOF
  
  # Determine new run id: one more than the number of initial run files
  initial_count=${#source_ids[@]}
  rerun_id=$(printf "%03d" $((initial_count + 1)))
  echo "Using run id $rerun_id for unmatched re-run."
  
  echo "Re-running SH-matching pipeline on unmatched sequences..."
  pushd "$SH_MATCHING_DIR" > /dev/null
  # Copy the new FASTA file into indata as source_<rerun_id>
  cp "$rerun_fasta" "indata/source_${rerun_id}"
  echo "Processing re-run pipeline with run id $rerun_id..."
  ./sh_matching_echo.sif /sh_matching/run_pipeline.sh "$rerun_id" itsfull no yes no no
  popd > /dev/null
  
  # Append new mapping line for re-run to the mapping file.
  # Extract the first header from the rerun FASTA file to record in the mapping file.
  first_header=$(grep '^>' "$rerun_fasta" | head -1)
  if [ -z "$first_header" ]; then
      first_header="-"
  fi
  # For SAMPLE_INFO, set a value such as "rerun_input".
  if [[ -s "$OUTPUT_DIR/${base%.*}_rerun.tsv" ]]; then
      match_status="rerun_match"
  else
      match_status="unmatched"
  fi
  echo -e "$INPUT_FILE\tindata/source_${rerun_id}\t$first_header\trerun_input\tsource_${rerun_id}\trerun_unmatched\t$match_status" >> "$MAPPING_FILE"
  
  echo "Aggregating output from re-run pipeline..."
  TMP_OUT_RERUN=$(mktemp -d)
  rerun_zip="$SH_MATCHING_DIR/outdata/source_${rerun_id}.zip"
  rerun_output="${OUTPUT_DIR}/${base%.*}_rerun.tsv"
  if [[ -f "$rerun_zip" ]]; then
      base_zip="source_${rerun_id}"
      target_dir="$TMP_OUT_RERUN/$base_zip"
      mkdir -p "$target_dir"
      unzip -o "$rerun_zip" -d "$target_dir"
      echo "Creating re-run TSV file: $rerun_output"
      > "$rerun_output"
      if [[ -f "$target_dir/matches/matches_out_all.csv" ]]; then
           cat "$target_dir/matches/matches_out_all.csv" >> "$rerun_output"
      else
           echo "No matches_out_all.csv found in re-run extraction."
      fi
      rm -rf "$TMP_OUT_RERUN"
  else
      echo "Re-run zip file not found: $rerun_zip"
      rerun_output=""
      rm -rf "$TMP_OUT_RERUN"  # Ensure removal even if zip file is missing
  fi
  
  echo "Creating new total TSV file..."
  cp "$matched_tsv" "$new_total_tsv"
  if [[ -n "$rerun_output" && -f "$rerun_output" ]]; then
      tail -n +2 "$rerun_output" >> "$new_total_tsv"
  fi
  
  echo "Re-run pipeline complete. Files generated:"
  echo "Original TSV: $output_file"
  echo "Unmatched TSV: $unmatched_tsv"
  echo "Matched TSV: $matched_tsv"
  echo "Re-run TSV: $rerun_output"
  echo "New Total TSV: $new_total_tsv"
  
  rm -f "$unmatched_ids_file"
else
  echo "Rerun of unmatched sequences is disabled. Skipping re-run pipeline."
fi

# -------------------------------------------------------------------
# (7) Cleanup: Remove Intermediate Files
# -------------------------------------------------------------------
echo "Cleaning up intermediate files..."

if [[ -f run_IDS_1.txt ]]; then
  rm -f run_IDS_1.txt
  echo "Removed run_IDS_1.txt."
fi

for id in "${source_ids[@]}"; do
  target="$SH_MATCHING_DIR/indata/$id"
  if [[ -f "$target" ]]; then
    rm -f "$target"
    echo "Removed intermediate file: $target"
  fi
done

# Also remove the rerun file from indata and its corresponding zip file from outdata (if rerun was enabled)
if [[ "$RERUN_UNMATCHED" == "yes" ]]; then
  rerun_indata="$SH_MATCHING_DIR/indata/source_${rerun_id}"
  if [[ -f "$rerun_indata" ]]; then
    rm -f "$rerun_indata"
    echo "Removed rerun intermediate file: $rerun_indata"
  fi
  
  rerun_outdata="$SH_MATCHING_DIR/outdata/source_${rerun_id}.zip"
  if [[ -f "$rerun_outdata" ]]; then
    rm -f "$rerun_outdata"
    echo "Removed rerun output zip file: $rerun_outdata"
  fi
fi

# Optionally, remove the mapping file if not needed.
# rm -f "$MAPPING_FILE"
# echo "Removed mapping file: $MAPPING_FILE."

echo "Processing complete."
