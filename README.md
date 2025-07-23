**!! NB This fork is in development for metabarcoding  !!**  
  
**Please use and [cite the main branch of SH-matching](https://github.com/TU-NHM/sh_matching_pub) if you intend to use this tool!!.**  

Adjustments include:  
-using all available threads.  
-usearch11 64 bit version (for large clusters).  
-replaced some csv dictionaries as they run into field limits.  
-added script for re-running unmatched seqs togheter (detect larger 'new_SH' across input files).  
-increased verbosity through echo statements.  
-script for metabarcoding to run a fasta file with automatic splitting, possibility to run split files in parallel on multiple compute nodes.

# SH MATCHING analysis tool

[![run with singularity](https://img.shields.io/badge/run%20with-singularity-blue?style=flat&logo=singularity)](https://sylabs.io/docs/)
[![Github_Status_Badge](https://img.shields.io/badge/GitHub-2.0.0-blue.svg)](https://github.com/TU-NHM/sh_matching_pub)
[![GitHub license](https://img.shields.io/github/license/TU-NHM/sh_matching_pub)](https://github.com/TU-NHM/sh_matching_pub/blob/master/LICENSE.md)

**NB! master branch is used as development branch. Please check out [Releases](https://github.com/TU-NHM/sh_matching_pub/releases) to download a specific version of the SH matching tool.**

Developed as part of [EOSC-Nordic](https://www.eosc-nordic.eu/) project (task 5.2.1: Cross-border data processing workflows), UNITE SH matching analysis is a digital service for the global species discovery from eDNA (environmental DNA). SH matching service is based on the [UNITE](https://unite.ut.ee) datasets hosted in [PlutoF](https://plutof.ut.ee). Its output includes information about what species are present in eDNA samples, are they potentially undescribed new species, where are they found in other studies, are they alien or threatened species, etc. The output will provide DOI (Digital Object Identifier) based stable identifiers for the communicating species found in eDNA. DOIs are connected to the taxonomic backbone of [PlutoF](https://plutof.ut.ee) and [GBIF](https://www.gbif.org). In this way every DOI is accompanied by a taxon name which is still widely used for the communication of species. In the case of undescribed species, DOIs will soon be issued by the [PlutoF](https://plutof.ut.ee) system (only if SH matching service integrated with the [PlutoF](https://plutof.ut.ee) platform is used for the analysis). SH matching service covers all Eukaryota by using rDNA ITS marker sequences accompanied by sample metadata.

The script expects input files in FASTA format. Outdata files are described in [sh_matching_analysis/readme.txt](https://github.com/TU-NHM/sh_matching_pub/blob/master/sh_matching_analysis/readme.txt).

## Third-party software used by this tool

* [USEARCH](https://www.drive5.com/usearch/)
* [VSEARCH](https://github.com/torognes/vsearch)
* [ITSx](https://microbiology.se/software/itsx/)
* [KronaTools](https://github.com/marbl/Krona/wiki/KronaTools)

## Setup

### Pre-requisites

* [Singularity](https://sylabs.io/singularity/) - install Singularity (tested with version 3.5) and obtain API key for remote build

### Setup steps

1. Create Singularity Image File (SIF)
```console
    git clone https://github.com/MycoMatics/sh_matching_pub_tweak
    export APPTAINER_TMPDIR=/tmp
    mkdir -p $APPTAINER_TMPDIR
    chmod 777 /tmp
    cd sh_matching_pub_tweak/
    apptainer build sh_matching_tweak.sif sh_matching.def 2>&1 | tee build.log
```
After cloning you may wish to track the original repository to stay up to date:
```bash
git remote add upstream https://github.com/TU-NHM/sh_matching_pub.git
git fetch upstream
git merge upstream/master     # or use rebase if preferred
```
If a merge reports conflicts, inspect the files listed by `git status`, edit them to resolve the differences and then complete the merge with `git add` and `git commit`.

2. Create input, output and working data directories
    ```console
    mkdir userdir
    mkdir indata
    mkdir outdata
    ```

3. Download FASTA dbs (https://app.plutof.ut.ee/filerepository/view/7502694) and create UDB formatted dbs
    ```console
    wget https://s3.hpc.ut.ee/plutof-public/original/5eba59a4-8c8d-4364-8da6-58676df02474.zip
    mv 5eba59a4-8c8d-4364-8da6-58676df02474.zip sh_matching_data_udb_0_5.zip
    unzip sh_matching_data_udb_0_5.zip
    rm sh_matching_data_udb_0_5.zip
    cd data_udb/
    ml load VSEARCH
    vsearch --makeudb_usearch sanger_refs_sh.fasta --output sanger_refs_sh.udb
    rm sanger_refs_sh.fasta
    vsearch --makeudb_usearch sanger_refs_sh_full.fasta --output sanger_refs_sh_full.udb
    rm sanger_refs_sh_full.fasta
    ml purge
    ```

## Running the analysis

**NB! The script expects input files in FASTA format, named as source_[run_id] and placed in indata/ directory. Outdata files are described in [sh_matching_analysis/readme.txt](https://github.com/TU-NHM/sh_matching_pub/blob/master/sh_matching_analysis/readme.txt).**

4. Run the pipeline using SIF (example data with -

* run_id=11
* region=itsfull[default]|its2
* itsx_step=yes[default]|no - flag indicating whether to include the ITSx step in the analysis (default, "yes")
* remove_userdir=yes[default]|no - flag indicating whether to delete the user directory upon pipeline completion (default, "yes")
* include_vsearch_step=yes|no[default] - flag indicating whether to include the vsearch substring dereplication step (default, "no")
* conduct_usearch_05_step=yes|no[default] - flag indicating whether to conduct the usearch complete-linkage clustering at 0.5% dissimilarity (default, "no")

    ```console
    ./sh_matching_tweak.sif /sh_matching/run_pipeline.sh 11 itsfull yes yes no no  2>&1 | tee sh_matching.log
    ```
5. For a metabarcoding flavour, run the (meta_SH_match.sh)[https://raw.githubusercontent.com/MycoMatics/sh_matching_pub_tweak/refs/heads/master/meta_SH_match.sh] script which splits large input fasta files into smaller chunks and runs SH-matching on each. The results are combined. It is possible to re-run unmatched sequences at the 0.5% toghether (to detect larger 'new_sh')
    ```console
       chmod +x meta_SH_match.sh
    
       bash meta_SH_match.sh
       Usage: meta_SH_match.sh -i INPUT_FILE -s SH_MATCHING_DIR -o OUTPUT_DIR [-t THREADS] [-m MAX_SEQUENCES] [-M MODE] [-n NODES] [-r RERUN_UNMATCHED]
      -i   Path to the input FASTA file
      -s   Directory for SH-matching files and tools
      -o   Output directory for results
      -t   Number of threads per task (default: all available cores via nproc)
      -m   Maximum sequences per split file (default: 30000)
      -M   Execution mode: "sequential" or "parallel" (default: sequential)
      -n   Number of nodes (used only in parallel mode, default: 1)
      -r   Re-run unmatched sequences: "yes" or "no" (default: yes)
      ```
   For example, to loop over input fasta files in a directory:
    ```console

    for FILE in /path/to/*.fasta; do
        base=$(basename "$FILE" .fasta)
        OUTPUT_DIR="prefix_${base}"
    
        # Check if the output directory already exists
        if [ -d "$OUTPUT_DIR" ]; then
            echo "Skipping $FILE as it has already been processed."
        else
            echo "Processing $FILE with output in $OUTPUT_DIR"
            chmod +x meta_SH_match.sh
            bash meta_SH_match.sh -i "$FILE" \
            -s "/path/to/sh_matching_pub_tweak" \
            -M "sequential"  -o "$OUTPUT_DIR"
        fi
    done
      ```
    Or if you have multiple computing nodes available, it is possible to use the SLURM job manager to run in parallel (i.e. split files of maximum input sequences, default=30k seqs).
    For example using 4 nodes:
    ```console

    for FILE in /path/to/*.fasta; do
        base=$(basename "$FILE" .fasta)
        OUTPUT_DIR="prefix_${base}"
    
        # Check if the output directory already exists
        if [ -d "$OUTPUT_DIR" ]; then
            echo "Skipping $FILE as it has already been processed."
        else
            echo "Processing $FILE with output in $OUTPUT_DIR"
            chmod +x meta_SH_match.sh
            bash meta_SH_match.sh -i "$FILE" \
            -s "/path/to/sh_matching_pub_tweak" \
            -M "parallel" -n 4 -o "$OUTPUT_DIR"
        fi
    done
      ```
BONUS:
Similar to meta_SH_match.sh script, there is a seperate script that can iterate across output files from SH-matching.
From the output, sequences can be merged across barcodes/input files based on SH-code. That is not the case for new/non-existing clusters. Those we need to gather from output files. Running [filter_newsh.py]() iterates output files, collects sequences that are either new_singletons or new_sh at the 0.5% level, write them to a single file and appends RUNID information to the fasta headers to seperate them out later.
You can use the file as a new input file for SH-matching and use the preferred cluster level resulting from this SH-matching to finalize your OTU-table.

## Citing

When using this resource, please cite as:

Abarenkov K, Kõljalg U, Nilsson RH (2022) UNITE Species Hypotheses Matching Analysis. Biodiversity Information Science and Standards 6: e93856. [https://doi.org/10.3897/biss.6.93856](https://doi.org/10.3897/biss.6.93856)

## Funding

The work is supported by [EOSC-Nordic](https://eosc-nordic.eu/) and the Estonian Research Council grant (PRG1170).
