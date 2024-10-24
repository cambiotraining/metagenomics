#!/bin/bash

set -e 

# Environments
mamba create -n alignment fastqc=0.12.1 cutadapt=4.9 trimmomatic=0.39 bowtie2=2.5.4 samtools=1.21 metaphlan=4.1.1 mash=2.3 multiqc==1.25.1

mamba create -n assembly fastqc=0.12.1 cutadapt=4.9 trimmomatic=0.39 bowtie2=2.5.4 samtools=1.21 spades=4.0.0 bbmap=39.10 flash=1.2.11 multiqc==1.25.1

mamba create -n mags maxbin2=2.2.7 prokka=1.14.6 gtdbtk=2.4.0 abricate=1.0.1 checkm-genome=1.2.3


# CheckM database

checkm_db="$HOME/Course_Materials/databases/checkmdb"
mkdir -p $checkm_db

wget -O checkm_db.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xzvf checkm_db.tar.gz -C $checkm_db
rm checkm_db.tar.gz

conda env config vars set CHECKM_DATA_PATH="$checkm_db" -n mags


# GTDB-tk database

gtdbtk_db="$HOME/Course_Materials/databases/gtdbtk"
mkdir -p $gtdbtk_db

wget -O gtdbtk_db.tar.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar -xzvf gtdbtk_db.tar.gz -C $gtdbtk_db
rm checkm_db.tar.gz

conda env config vars set GTDBTK_DATA_PATH="$gtdbtk_db" -n mags


# MetaPhlAn database

mamba activate alignment

metaphlan_db="$HOME/databases/metaphlan"
mkdir -p $metaphlan_db

metaphlan --install --bowtie2db $metaphlan_db

conda env config vars set DEFAULT_DB_FOLDER="$metaphlan_db" -n alignment