**!! This fork is in development for trying out some semi-streamlined metabarcoding workflow !! Please use and [cite the main branch of SH-matching](https://github.com/TU-NHM/sh_matching_pub) if you intend to use this tool!!.**

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
    apptainer build sh_matching.sif sh_matching.def
    ```

2. Create input, output and working data directories
    ```console
    mkdir userdir
    mkdir indata
    mkdir outdata
    ```

3. Download FASTA dbs (https://app.plutof.ut.ee/filerepository/view/6884701) and create UDB formatted dbs
    ```console
    wget https://s3.hpc.ut.ee/plutof-public/original/d3d8b3de-83af-4fb5-b82b-359f7b730f84.zip
    mv d3d8b3de-83af-4fb5-b82b-359f7b730f84.zip sh_matching_data_udb_0_5.zip
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
    ./sh_matching.sif /sh_matching/run_pipeline.sh 11 itsfull yes yes no no
    ```
5. For a metabarcoding flavour, rename your sourcefiles according to your pre-processed barcode files (e.g. barcode 28 could be source_28). Create a list of RUNID in a txt file and to run in parallel on mutliple nodes using SLURM, do:
    ```console
    srun --nodes=node_amount --ntasks=node_amount --label /bin/bash -c '
      # Generate a timestamp
      TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
      # Duplicate all output (stdout and stderr) to a log file and the screen:
      exec > >(tee -a job_${SLURM_PROCID}_$(hostname)_${TIMESTAMP}.log) 2>&1
      export OMP_NUM_THREADS=max_threads_per_node
      # Each node processes X input files sequentially:
      RUNID1=$(sed -n "$((SLURM_PROCID+1))p" source_numbers.txt)
      RUNID2=$(sed -n "$((SLURM_PROCID*2+2))p" source_numbers.txt)
      echo "Task $SLURM_PROCID on $(hostname) processing RUNID1: $RUNID1"
      ./sh_matching.sif /sh_matching/run_pipeline.sh $RUNID1 itsfull no yes no no
      echo "Task $SLURM_PROCID on $(hostname) processing RUNID1: $RUNID2"
      ./sh_matching.sif /sh_matching/run_pipeline.sh $RUNID2 itsfull no yes no no
    '
    ```
    The number of RUNID depends on how many input files you want to process per node. For example using 7 nodes, to process 14 source files you would need 14/7=2 RUNIDs

6. From the output, sequences can be merged across barcodes/input files based on SH-code. That is not the case for new/non-existing clusters. Those we need to gather from output files. Running [filter_newsh.py]() iterates output files, collects sequences that are either new_singletons or new_sh at the 0.5% level, write them to a single file and appends RUNID information to the fasta headers to seperate them out later.
You can use the file as a new input file for SH-matching and use the preferred cluster level resulting from this SH-matching to finalize your OTU-table.

## Citing

When using this resource, please cite as:

Abarenkov K, Kõljalg U, Nilsson RH (2022) UNITE Species Hypotheses Matching Analysis. Biodiversity Information Science and Standards 6: e93856. [https://doi.org/10.3897/biss.6.93856](https://doi.org/10.3897/biss.6.93856)

## Funding

The work is supported by [EOSC-Nordic](https://eosc-nordic.eu/) and the Estonian Research Council grant (PRG1170).
