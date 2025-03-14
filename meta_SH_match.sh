#!/bin/bash
set -euo pipefail
# -------------------------------------------------------------------
# author: Glen Dierickx  March 2025
#
# Description
# script that takes a fasta file and does SH-matching on it, designed to process large (metabarcoding) input files,
# the script automatially splits the input fasta into default 30k seqs and rearranges the output in a single tsv table in a specified output dir.
# unmatched seqs at 0.5% can be rerun togheter to detect large non-existing SH and in that way, serve as a proxy to unidentified OTUs in analyses.
# To speed up there is a possibility to run the SH-matching on split files in parallel on different NODES (!not different threads!) using SLURM job manager.
# This script works on the Ghent University HPC and was not tested outside of that environment and was created because of both inode file issues
# and memory issues resulting from usearch hierarchical clustering which needs to allocate a full distance matrix into working mem (results in error filesize too big)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Usage function
# -------------------------------------------------------------------
usage() {
  echo "Usage: $0 -i INPUT_FILE -s SH_MATCHING_DIR -o OUTPUT_DIR [-t THREADS] [-m MAX_SEQUENCES] [-M MODE] [-n NODES] [-r RERUN_UNMATCHED]"
  echo ""
  echo "  -i   Path to the input FASTA file"
  echo "  -s   Directory for SH-matching files and tools"
  echo "  -o   Output directory for results"
  echo "  -t   Number of threads per task (default: all available cores via nproc)"
  echo "  -m   Maximum sequences per split file (default: 30000)"
  echo "  -M   Execution mode: \"sequential\" or \"parallel\" (default: sequential)"
  echo "  -n   Number of nodes (used only in parallel mode, default: 1)"
  echo "  -r   Re-run unmatched sequences: \"yes\" or \"no\" (default: yes)"
  exit 1
}

# -------------------------------------------------------------------
# Parse command-line arguments with getopts
# -------------------------------------------------------------------
THREADS=$(nproc)
MAX_SEQUENCES=30000
MODE="sequential"
NODES=1
RERUN_UNMATCHED="yes"

while getopts "i:t:m:M:n:s:o:r:" opt; do
  case "$opt" in
    i) INPUT_FILE="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MAX_SEQUENCES="$OPTARG" ;;
    M) MODE="$OPTARG" ;;
    n) NODES="$OPTARG" ;;
    s) SH_MATCHING_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    r) RERUN_UNMATCHED="$OPTARG" ;;
    *) usage ;;
  esac
done

if [[ -z "${INPUT_FILE:-}" || -z "${SH_MATCHING_DIR:-}" || -z "${OUTPUT_DIR:-}" ]]; then
  usage
fi

# -------------------------------------------------------------------
# Logging and initial configuration
# -------------------------------------------------------------------
if [[ -d "$OUTPUT_DIR" ]]; then
  echo "WARNING: Output directory $OUTPUT_DIR already exists. Clearing its contents..."
  rm -rf "$OUTPUT_DIR"/*
fi
mkdir -p "$OUTPUT_DIR"

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
# (A) Determine Offset for New Source IDs in indata
# -------------------------------------------------------------------
existing_files=("$SH_MATCHING_DIR/indata/source_"*)
max_existing=0
for file in "${existing_files[@]}"; do
  if [[ -f "$file" ]]; then
    fname=$(basename "$file")
    id_num=$((10#${fname#source_}))
    if (( id_num > max_existing )); then
      max_existing=$id_num
    fi
  fi
done
offset=$((max_existing + 10))
echo "Offset for new source ids is $offset"

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

SAMPLE_FLAG=0
if grep -m 1 "sample=" "$INPUT_FILE" > /dev/null; then
  SAMPLE_FLAG=1
  echo "Detected 'sample=' in the input file; sample info will be mapped."
else
  echo "No 'sample=' string detected in the input file."
fi

TMP_SPLIT_DIR=$(mktemp -d)
MAPPING_FILE="$OUTPUT_DIR/mapping_${TIMESTAMP}.txt"
echo -e "INPUT_FILE\tSOURCE_FILE\tHEADER\tSAMPLE_INFO\tSOURCE_NAME\tRUN_TYPE\tMATCH_STATUS_005" > "$MAPPING_FILE"

awk -v max_seq="$MAX_SEQUENCES" \
    -v outdir="$TMP_SPLIT_DIR" \
    -v prefix="source_" \
    -v mapping_file="$MAPPING_FILE" \
    -v input_file="$INPUT_FILE" \
    -v sample_flag="$SAMPLE_FLAG" \
    -v offset="$offset" '
BEGIN { seqCount = 0; fileCount = 0; }
/^>/ {
    seqCount++;
    if ((seqCount - 1) % max_seq == 0) {
        fileCount++;
        id = sprintf("%03d", fileCount + offset);
        current_file = outdir "/" prefix id;
    }
    print $0 > current_file;
    sample = "";
    if (sample_flag == 1) {
        n = split($0, fields, ";")
        for (i = 1; i <= n; i++) {
            if (fields[i] ~ /^sample=/) {
                sub(/^sample=/, "", fields[i])
                sample = fields[i];
                break;
            }
        }
    }
    srcname = current_file; sub(".*/", "", srcname);
    print input_file "\t" current_file "\t" $0 "\t" sample "\t" srcname "\tinitial_run\tunmatched" >> mapping_file;
    next;
}
{ print $0 >> current_file; }
END { print "Created " fileCount " split files." > "/dev/stderr"; }
' "$INPUT_FILE"

echo "Moving split files to $SH_MATCHING_DIR/indata..."
mv "$TMP_SPLIT_DIR"/source_* "$SH_MATCHING_DIR/indata/"
rm -r "$TMP_SPLIT_DIR"

# -------------------------------------------------------------------
# Prepare List of Unique Source Files (from mapping) to Process
# -------------------------------------------------------------------
mapfile -t source_ids < <(awk 'NR>1 { id = $2; sub(".*/", "", id); if (!seen[id]++) print id }' "$MAPPING_FILE")
echo "Unique source IDs to process: ${source_ids[*]}"

# Export SOURCE_IDS and OMP_NUM_THREADS for use by parallel tasks.
export SOURCE_IDS="${source_ids[*]}"
export OMP_NUM_THREADS="$THREADS"

# -------------------------------------------------------------------
# (5) Running SH-Matching Pipeline (Initial Run)
# -------------------------------------------------------------------
if [[ "$MODE" == "sequential" ]]; then
  echo "Running in sequential mode..."
  pushd "$SH_MATCHING_DIR" > /dev/null
  for id in "${source_ids[@]}"; do
    num_id="${id#source_}"
    echo "Processing $id sequentially with run id $num_id..."
    /scratch/gent/vo/001/gvo00142/sh_matching_pub_tweak/sh_matching_echo_rawumi.sif /sh_matching/run_pipeline.sh "$num_id" itsfull no yes yes no
  done
  popd > /dev/null
elif [[ "$MODE" == "parallel" ]]; then
  if (( NODES <= 1 )); then
    echo "Error: Parallel mode requires more than 1 node. Skipping pipeline execution."
  else
    echo "Running in parallel mode on $NODES nodes..."
    pushd "$SH_MATCHING_DIR" > /dev/null

    # --- Revised srun block for debugging ---
    PARALLEL_COMMAND=$(cat <<'EOF'
set -x
# Set defaults for SLURM variables if not set.
: ${SLURM_NTASKS:=1}
: ${SLURM_PROCID:=0}
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
echo "DEBUG: Current time is $TIMESTAMP"
# Reconstruct the source IDs array from the environment variable.
read -r -a ids <<< "$SOURCE_IDS"
echo "DEBUG: Source IDs received: ${ids[@]}"
TOTAL_FILES=${#ids[@]}
echo "DEBUG: TOTAL_FILES = $TOTAL_FILES"
BASE_FILES=$(( TOTAL_FILES / SLURM_NTASKS ))
REMAINDER=$(( TOTAL_FILES % SLURM_NTASKS ))
echo "DEBUG: BASE_FILES = $BASE_FILES, REMAINDER = $REMAINDER"
if [ "$SLURM_PROCID" -lt "$REMAINDER" ]; then
    start_index=$(( SLURM_PROCID * (BASE_FILES + 1) ))
    nfiles=$(( BASE_FILES + 1 ))
else
    start_index=$(( REMAINDER * (BASE_FILES + 1) + (SLURM_PROCID - REMAINDER) * BASE_FILES ))
    nfiles=$(( BASE_FILES ))
fi
end_index=$(( start_index + nfiles ))
echo "DEBUG: Task \$SLURM_PROCID on \$(hostname) assigned indices: \$start_index to \$end_index"
for (( i = start_index; i < end_index; i++ )); do
    echo "DEBUG: In for loop, i=\$i"
    RUNID="${ids[i]}"
    if [ -z "$RUNID" ]; then
        echo "DEBUG: No RUNID found at index \$i. Terminating loop."
        break
    fi
    num_id="${RUNID#source_}"
    echo "DEBUG: Task \$SLURM_PROCID on \$(hostname) processing RUNID index \$i: \$num_id"
    /scratch/gent/vo/001/gvo00142/sh_matching_pub_tweak/sh_matching_echo_rawumi.sif /sh_matching/run_pipeline.sh "$num_id" itsfull no yes yes no
done
echo "DEBUG: Done with parallel processing"
EOF
)

echo "DEBUG: The following parallel command will be executed by srun:"
echo "$PARALLEL_COMMAND"

srun --nodes="$NODES" --ntasks="$NODES" --export=ALL --label /bin/bash -c "$PARALLEL_COMMAND" || { echo "srun command failed, but continuing..."; }



    popd > /dev/null
  fi
else
  echo "Invalid mode: $MODE. Exiting."
  exit 1
fi

# -------------------------------------------------------------------
# (6) Output Aggregation: Unzip and Concatenate as TSV (Initial Run)
# -------------------------------------------------------------------
echo "Concatenating output files from current run..."
TMP_OUT=$(mktemp -d)
for id in "${source_ids[@]}"; do
    zipfile="$SH_MATCHING_DIR/outdata/${id}.zip"
    if [[ -f "$zipfile" ]]; then
        target_dir="$TMP_OUT/${id}"
        echo "Unzipping $zipfile into $target_dir..."
        mkdir -p "$target_dir"
        unzip -o "$zipfile" -d "$target_dir"
    else
        echo "Zip file for ${id} not found; skipping."
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
rm -rf "$TMP_OUT"

# -------------------------------------------------------------------
# (6.5) Re-run Pipeline for Unmatched Sequences (if enabled)
# -------------------------------------------------------------------
if [[ "$RERUN_UNMATCHED" == "yes" ]]; then
  echo "Filtering unmatched sequences from original TSV..."
  unmatched_tsv="${OUTPUT_DIR}/${base%.*}_unmatched.tsv"
  matched_tsv="${OUTPUT_DIR}/${base%.*}_matched.tsv"
  new_total_tsv="${OUTPUT_DIR}/${base%.*}_new_total.tsv"
  
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
            rec.id = rec.id
            rec.description = rec.id
            count += 1
            SeqIO.write(rec, outfile, "fasta")
print(f"Extracted {count} unmatched sequences to {rerun_fasta}")
EOF
  
  existing_files=("$SH_MATCHING_DIR/indata/source_"*)
  max_id=0
  for file in "${existing_files[@]}"; do
    if [[ -f "$file" ]]; then
      fname=$(basename "$file")
      id_num=$((10#${fname#source_}))
      if (( id_num > max_id )); then
        max_id=$id_num
      fi
    fi
  done
  rerun_id=$(printf "%03d" $((max_id + 1)))
  echo "Using run id $rerun_id for unmatched re-run."
  
  echo "Re-running SH-matching pipeline on unmatched sequences..."
  pushd "$SH_MATCHING_DIR" > /dev/null
  cp "$rerun_fasta" "indata/source_${rerun_id}"
  
  echo "Processing re-run pipeline with run id $rerun_id..."
  (
    /scratch/gent/vo/001/gvo00142/sh_matching_pub_tweak/sh_matching_echo_rawumi.sif /sh_matching/run_pipeline.sh "$rerun_id" itsfull no yes yes yes
  )
  popd > /dev/null
  
  first_header=$(grep '^>' "$rerun_fasta" | head -1)
  if [ -z "$first_header" ]; then
      first_header="-"
  fi
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
      rm -rf "$TMP_OUT_RERUN"
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
echo "Removing files from indata..."
for id in "${source_ids[@]}"; do
  target="$SH_MATCHING_DIR/indata/$id"
  if [[ -f "$target" ]]; then
    rm -f "$target"
    echo "Removed indata file: $target"
  else
    echo "Indata file $target not found."
  fi
done

echo "Removing files from outdata..."
for id in "${source_ids[@]}"; do
  target="$SH_MATCHING_DIR/outdata/${id}.zip"
  if [[ -f "$target" ]]; then
    rm -f "$target"
    echo "Removed outdata zip file: $target"
  else
    echo "Outdata zip file $target not found."
  fi
done

if [[ "$RERUN_UNMATCHED" == "yes" ]]; then
  rerun_indata="$SH_MATCHING_DIR/indata/source_${rerun_id}"
  if [[ -f "$rerun_indata" ]]; then
    rm -f "$rerun_indata"
    echo "Removed rerun indata file: $rerun_indata"
  fi
  
  rerun_outdata="$SH_MATCHING_DIR/outdata/source_${rerun_id}.zip"
  if [[ -f "$rerun_outdata" ]]; then
    rm -f "$rerun_outdata"
    echo "Removed rerun outdata zip file: $rerun_outdata"
  fi
fi

echo "Processing complete."
