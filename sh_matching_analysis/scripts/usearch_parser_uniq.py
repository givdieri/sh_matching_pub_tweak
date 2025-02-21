import argparse
import csv
import logging
import os
from pathlib import Path

from Bio import SeqIO

parser = argparse.ArgumentParser(description="Script to create mapping file for duplicates")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")

# infiles
cov100_uniq_file = user_dir / f"source_{run_id}_fastanames"

# outfiles
duplic_seqs_file = user_dir / "duplic_seqs.txt"

log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

seq_counter = 0  # sequence count
seq_counter_all = 0  # sequence count with duplicates included

# include duplicate sequences (tmp_files/sequences.names)
print("Collecting duplicate sequences (cov100) ...")
cov100_uniq_dict = {}
with open(cov100_uniq_file, "r") as f, open(duplic_seqs_file, "w") as dupl:
    for line in f:
        seq_counter += 1
        # Remove trailing newline characters
        line = line.rstrip("\n")
        # Split on the first tab only (guarantees two columns)
        parts = line.split("\t", 1)
        if len(parts) < 2:
            continue  # skip lines that don't have two columns
        key, value = parts
        # Process only rows where the key and value differ
        if key != value:
            # Split the second column by commas
            cov100_list = value.split(",")
            try:
                cov100_list.remove(key)
            except ValueError:
                pass  # key might not be in the list; ignore if it's not
            for seq in cov100_list:
                seq_counter_all += 1
                dupl.write(f"{seq}\t{key}\t\n")
                
logging.info(f"USEARCH_PARSER_UNIQ\tNumber of sequences at the end: {seq_counter} and {seq_counter_all}")
