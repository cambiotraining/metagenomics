---
title: "Data & Setup"
number-sections: false
---

<!-- 
Note for Training Developers:
We provide instructions for commonly-used software as commented sections below.
Uncomment the sections relevant for your materials, and add additional instructions where needed (e.g. specific packages used).
Note that we use tabsets to provide instructions for all three major operating systems.
-->

::: {.callout-tip level=2}
## Workshop Attendees

If you are attending one of our workshops, we will provide a training environment with all of the required software and data.  
If you want to setup your own computer to run the analysis demonstrated on this course, you can follow the instructions below.
:::

## Data

You can download the data used in the workshop from the following link:

<a href="https://www.dropbox.com/scl/fo/fmzzgerrv98plnklzthh0/h?rlkey=x0zvsweshldl4m056q48jqe5d&st=qmrtfhso&dl=0">
  <button class="btn"><i class="fa fa-download"></i> Download</button>
</a>


## Software

### General setup

::: {.panel-tabset group="os"}
#### Windows

- Setup the **Windows Subsystem for Linux (WSL)** following [these instructions](https://cambiotraining.github.io/software-installation/materials/wsl.html).
- From a WSL terminal (previous step), install **Mamba** following the **Linux instructions** on [this page](https://cambiotraining.github.io/software-installation/materials/mamba.html).
- Install **R** and **RStudio** following [these instructions](https://cambiotraining.github.io/software-installation/materials/r-base.html).


#### Mac OS

- **Setup your macOS** by following [these instructions](https://cambiotraining.github.io/software-installation/materials/macos.html).
- Install **Mamba** using [these instructions](https://cambiotraining.github.io/software-installation/materials/mamba.html).
- Install **R** and **RStudio** following [these instructions](https://cambiotraining.github.io/software-installation/materials/r-base.html).

#### Linux

- Install **Mamba** using [these instructions](https://cambiotraining.github.io/software-installation/materials/mamba.html).
- Install **R** and **RStudio** following [these instructions](https://cambiotraining.github.io/software-installation/materials/r-base.html).

:::

### R Packages

Open RStudio. 
In the R console, run the following commands to install all the necessary packages: 

```r
install.packages("BiocManager")
BiocManager::install(c("dada2", 
                       "phyloseq", 
                       "Biostrings", 
                       "ggplot2", 
                       "reshape2", 
                       "readxl", 
                       "tidyverse"))
```


### Bioinformatics software

We can install the software used in the course using `mamba`. 
Due to the large number of programs, we recommend installing them in separate environments to avoid package version conflicts. 
The following commands install the latest version of each software at the time of writing. 
You may want to search [anaconda.org](https://anaconda.org/) for the latest versions available. 

```bash
mamba create -n alignment fastqc=0.12.1 cutadapt=4.9 trimmomatic=0.39 bowtie2=2.5.4 samtools=1.21 metaphlan=4.1.1 mash=2.3 multiqc==1.25.1

mamba create -n assembly fastqc=0.12.1 cutadapt=4.9 trimmomatic=0.39 bowtie2=2.5.4 samtools=1.21 spades=4.0.0 bbmap=39.10 flash=1.2.11 multiqc==1.25.1

mamba create -n mags maxbin2=2.2.7 prokka=1.14.6 gtdbtk=2.4.0 abricate=1.0.1 checkm-genome=1.2.3
```

From now on, you can use these packages, by activating the respective software environment using `mamba activate alignment`, `mamba activate assembly` or `mamba activate mags`.


### Databases

Some of the programs used require us to download public databases in addition to their installation. 
These files can be quite large, so we recommend that you use a shared storage if you're working in a team. 
We also recommend that you keep track of the database versions used (e.g. saving them in explicit folder names), in case new updates are released in the future and you want to reproduce an analysis. 

#### CheckM (1.5 GiB)

First activate the environment: 

```bash
mamba activate mags
```

The [CheckM documentation](https://github.com/Ecogenomics/CheckM/wiki/Installation#required-reference-data) gives the link to its database file. 

We will download this databases to a directory in our home called `~/databases/checkmdb_20150116`, but you can change this if you prefer to save it elsewhere. 
We use the date of the latest version of the database in the directory name for reference.

```bash
# create variable with output directory name for our database
# change this to be a directory of your choice
checkm_db="$HOME/databases/checkmdb_20150116"
mkdir -p $checkm_db
```

Download and decompress the file:

```bash
wget -O checkm_db.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xzvf checkm_db.tar.gz -C $checkm_db
rm checkm_db.tar.gz
```

After downloading, you need to run the following command to configure CheckM:

```bash
checkm data setRoot $checkm_db
```

Alternatively, you can set an environment variable specifically in your Conda/Mamba environment: 

```bash
conda env config vars set CHECKM_DATA_PATH="$checkm_db" -n mags
```


#### GTDB-Tk (40GB)

The [GTDB-tk documentation](https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data) gives the link to its database files. 

We will download this databases to a directory in our home called `~/databases/gtdbtk_r220`, but you can change this if you prefer to save it elsewhere. 
We use the name of the latest version of the database in the directory name for our reference.

```bash
# create variable with output directory name for our database
# change this to be a directory of your choice
gtdbtk_db="$HOME/databases/gtdbtk_r220"
mkdir -p $gtdbtk_db
```

Download and decompress the file:

```bash
wget -O gtdbtk_db.tar.gz https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar -xzvf gtdbtk_db.tar.gz -C $gtdbtk_db
rm checkm_db.tar.gz
```

Finally, we need to configure an environment variable to tell GTDB-tk where to find the database. 
We define this for our Conda/Mamba environment called `mags`:

```bash
conda env config vars set GTDBTK_DATA_PATH="$gtdbtk_db" -n mags
```

#### MetaPhlAn (24 GiB)

MetaPhlAn provides a command to download the latest database from its server ([instructions](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4#pre-requisites)). 
First we activate our environment and create the directory for the database: 

```bash
mamba activate alignment

# create variable with output directory name for our database
# change this to be a directory of your choice
metaphlan_db="$HOME/databases/metaphlan"
mkdir -p $metaphlan_db
```

Then we can run the download command (this can take a long time to finish):

```bash
metaphlan --install --bowtie2db $metaphlan_db
```

Finally, we need to configure an environment variable to tell MetaPhlAn where to find the database. 
We define this for our Conda/Mamba environment called `alignment`:

```bash
conda env config vars set DEFAULT_DB_FOLDER="$metaphlan_db" -n alignment
```