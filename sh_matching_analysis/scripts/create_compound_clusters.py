import argparse
import csv
import gzip
import logging
import os
import subprocess
from pathlib import Path

from Bio import SeqIO

# for debugging add echo of subscript name
if '__file__' in globals():
    script_name = os.path.basename(__file__)
else:
    script_name = 'Interactive session or unknown'

echo_message = f"Running script: {script_name}"
print(echo_message)
# echo done 

parser = argparse.ArgumentParser(description="Script to create compound clusters for hits")
parser.add_argument("run_id", help="Need run id in numeric format!")
args = parser.parse_args()

# read in args
run_id = args.run_id
if not run_id.isdigit():
    raise ValueError("Run id is not numeric", run_id)

user_dir = Path(f"{os.getcwd()}/userdir/{run_id}")
data_dir = Path("/sh_matching/data")
infile = user_dir / "iupac_out_full.fasta"
tmp_infile = user_dir / "closedref.80-best-hits.map.uc"

sanger_refs_file = data_dir / "sanger_refs_sh_full.fasta"
true_compound2seq_file = data_dir / "compound2seq_mapping.txt"
compound2seq_file = data_dir / "sh030_2seq_mapping.txt"
sh2compound_file = data_dir / "sh005_to_sh030_mappings.txt"

# Logging conf
log_file = user_dir / f"err_{run_id}.log"
logging.basicConfig(
    filename=log_file, filemode="a", format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level="INFO",
)

uniques = {}
sh_ucl_dict = {}
seq_ucl_dict = {}
true_seq_ucl_dict = {}
original = {}
ucl_seq_dict = {}
seq_ucl_map_dict = {}
true_seq_ucl_map_dict = {}
seq_ucl_length_dict = {}
must_refs_dict = {}
true_sh2compound_dict = {}

with open(sanger_refs_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        name = record.id
        fields = name.split("_")
        new_name = fields[0]
        uniques[new_name] = str(record.seq)
        must_refs_dict[new_name] = 1

with open(infile, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        uniques[record.id] = str(record.seq)

# get compound and SH mappings
with open(sh2compound_file) as sh2c:
    dataReader_sh2c = csv.DictReader(sh2c, delimiter="\t", fieldnames=["field_1", "field_2", "field_3"])
    for row in dataReader_sh2c:
        sh_ucl_dict[row["field_1"]] = row["field_2"]

with open(true_compound2seq_file) as true_c2seq:
    # UCL10_000001     1465615 CCCAT...     1327    SH1392229.10FU
    true_dataReader_c2seq = csv.reader(true_c2seq, delimiter="\t")
    for row in true_dataReader_c2seq:
        new_name = f"i{row[1]}i"
        true_sh2compound_dict[row[4]] = row[0]
        if row[4] in sh_ucl_dict:
            if row[0] in true_seq_ucl_map_dict:
                if sh_ucl_dict[row[4]] in true_seq_ucl_map_dict[row[0]]:
                    if not row[0] in true_seq_ucl_dict:
                        true_seq_ucl_dict[row[0]] = true_seq_ucl_map_dict[row[0]][sh_ucl_dict[row[4]]]
                    true_seq_ucl_map_dict[row[0]][sh_ucl_dict[row[4]]] = new_name
                else:
                    true_seq_ucl_map_dict[row[0]][sh_ucl_dict[row[4]]] = new_name
            else:
                true_seq_ucl_map_dict[row[0]] = {}
                true_seq_ucl_map_dict[row[0]][sh_ucl_dict[row[4]]] = new_name
            seq_string = row[2].replace("-", "").replace(" ", "").replace("\r\n", "").upper()
            original[new_name] = seq_string

with open(compound2seq_file) as c2seq:
    # SH0000001.10FU  6816002 ACACT...        480     SH1434604.10FU  UCL10_000001
    dataReader_c2seq = csv.reader(c2seq, delimiter="\t")
    for row in dataReader_c2seq:
        new_name = f"i{row[1]}i"

        if row[0] in seq_ucl_map_dict:
            if row[4] in seq_ucl_map_dict[row[0]]:
                if row[0] in seq_ucl_dict:
                    seq_ucl_dict[row[0]] = f"{seq_ucl_dict[row[0]]},{seq_ucl_map_dict[row[0]][row[4]]}"
                else:
                    seq_ucl_dict[row[0]] = seq_ucl_map_dict[row[0]][row[4]]
                seq_ucl_map_dict[row[0]][row[4]] = new_name
            else:
                seq_ucl_map_dict[row[0]][row[4]] = new_name
        else:
            seq_ucl_map_dict[row[0]] = {}
            seq_ucl_map_dict[row[0]][row[4]] = new_name
        seq_string = row[2].replace("-", "").replace(" ", "").replace("\r\n", "").upper()
        original[new_name] = seq_string

for ucl in seq_ucl_map_dict:
    if ucl in seq_ucl_dict:
        ucl_present_arr = seq_ucl_dict[ucl].split(",")
    else:
        ucl_present_arr = []
    for sh in seq_ucl_map_dict[ucl]:
        if ucl in seq_ucl_dict:
            seq_ucl_dict[ucl] = f"{seq_ucl_dict[ucl]},{seq_ucl_map_dict[ucl][sh]}"
        else:
            seq_ucl_dict[ucl] = seq_ucl_map_dict[ucl][sh]
        ucl_present_arr.append(seq_ucl_map_dict[ucl][sh])

        if true_sh2compound_dict[sh] in true_seq_ucl_map_dict:
            for sh_050 in true_seq_ucl_map_dict[true_sh2compound_dict[sh]]:
                if not true_seq_ucl_map_dict[true_sh2compound_dict[sh]][sh_050] in ucl_present_arr:
                    seq_ucl_dict[ucl] = f"{seq_ucl_dict[ucl]},{true_seq_ucl_map_dict[true_sh2compound_dict[sh]][sh_050]}"
                    ucl_present_arr.append(true_seq_ucl_map_dict[true_sh2compound_dict[sh]][sh_050])

with open(tmp_infile) as tmp_inf:
    dataReader_tmp_inf = csv.reader(tmp_inf, delimiter="\t")
    for row in dataReader_tmp_inf:
        if row[0] == "H":
            fields2 = row[9].split("_")
            if fields2[1] in sh_ucl_dict:
                if sh_ucl_dict[fields2[1]] in ucl_seq_dict:
                    ucl_seq_dict[sh_ucl_dict[fields2[1]]] = ucl_seq_dict[sh_ucl_dict[fields2[1]]] + "," + row[8]
                else:
                    ucl_seq_dict[sh_ucl_dict[fields2[1]]] = row[8]
            else:
                logging.info(f"COMP\tUCL for {fields2[0]} not found.")

# print out mapping table
for key, value in ucl_seq_dict.items():
    compound_file = user_dir / "compounds" / f"{key}.fas"
    compound_file_gz = user_dir / "compounds" / f"{key}.fas.gz"
    with open(compound_file, "w") as c:
        seqs = value.split(",")
        for i, val in enumerate(seqs):
            c.write(f">{val}\n")
            c.write(f"{uniques[val]}\n")

        seqs2 = seq_ucl_dict[key].split(",")
        for i, val in enumerate(seqs2):
            c.write(f">{val}\n")
            c.write(f"{original[val]}\n")

    # gzip the file
    with open(compound_file, 'rb') as c:
        with gzip.open(compound_file_gz, 'wb') as z:
            z.writelines(c)

    # rm original file
    rm_cmd_1 = subprocess.run(["rm", str(compound_file)], stdout=subprocess.DEVNULL)
