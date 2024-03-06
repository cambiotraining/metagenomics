---
title: Practical
---

## Discover an unknown genome in low complexity mixed community

In this practical we are simulating a virus infection that is caused by a completely unknown virus. In this simulated dataset a novel genome is hidden in the raw data with a huge background of human genome sequence. This simulates a real scenario when a pathogen is in the blood or other (otherwise sterile) body fluid. To make the bioinformatics step faster we generated the human "background" from chromosome 22, so the database will be relatively small.

### QC and Pre-processing

The raw date quality control and pre-processing is going the same way as we did with the mixed community data, for the details on these steps, please refer to [Day 2 practical material](22-pract.html#standard-quality-control-and-pre-processing-of-shotgun-metagenomics-raw-data).

```bash
# Deactivate the metagenomics environment if you are in that
conda deactivate

cd sg_raw_data/
fastqc unknown_pathogen_R1.fastq unknown_pathogen_R2.fastq

conda activate metagenomics

cutadapt -a CTGTCTCTTATACACATCT -A ATGTGTATAAGAGACA \
-o unknown_pathogen_noadapt_R1.fastq -p unknown_pathogen_noadapt_R2.fastq \
unknown_pathogen_R1.fastq unknown_pathogen_R2.fastq

trimmomatic PE -phred33 \
unknown_pathogen_noadapt_R1.fastq unknown_pathogen_noadapt_R2.fastq \
unknown_forward_paired.fq.gz unknown_forward_unpaired.fq.gz \
unknown_reverse_paired.fq.gz unknown_reverse_unpaired.fq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### Aligning reads to human chromosome 22

First let's try to align our raw data to the known genome content (in this case the human chromosome 22 sequence). This is again going similarly to our Day 2 practical, first we create a Bowtie2 database, then perform the alignment.

```bash
cd sg_reference/
bowtie2-build -q Homo_sapiens.GRCh38.dna.chromosome.22.fa human_chr22

cd ..

bowtie2 --qc-filter -p 8 --local -x sg_reference/human_chr22 \
-1 sg_raw_data/unknown_forward_paired.fq.gz \
-2 sg_raw_data/unknown_reverse_paired.fq.gz \
-S results/unknown_pathogen.sam

# Inspect the result .sam file
less results/unknown_pathogen.sam
```

If we look at the alignment result, we can see the aligning short reads to the human chromosome 22 reference sequence. Reads that are not aligning to the reference genome are also in the SAM file, but in this way not really in a usable format.

::: {.callout-exercise}
#### Extract reads that are not aligning to the human chromosome 22 for further processing
{{< level 3 >}}
Refer to the command line help (`bowtie2 -h`) to read about the numerous options the application provides. Re-run the bowtie2 alignment step with the following modifications:
- use all output fasta files from the `trimmomatic` output as unpaired reads
- Find and use the option that prints out short reads that are **not aligning** to the reference genome
- As we don't need the alignment SAM file, find a way to "throw out" the SAM file during the alignment
- name the output FASTQ file as `unknown_enriched.fastq`

::: {.callout-answer}
```bash
bowtie2 --qc-filter -p 8 --local \
-x sg_reference/human_chr22 \
-U sg_raw_data/unknown_forward_paired.fq.gz,sg_raw_data/unknown_reverse_paired.fq.gz,sg_raw_data/unknown_forward_unpaired.fq.gz,sg_raw_data/unknown_reverse_unpaired.fq.gz \
--un sg_raw_data/unknown_enriched.fastq > /dev/null
```

Alternatively, you can concatenate all input FASTQ files into one input file to avoid the long comma separated list.

```bash
cat sg_raw_data/unknown_forward_paired.fq.gz \
sg_raw_data/unknown_reverse_paired.fq.gz \
sg_raw_data/unknown_forward_unpaired.fq.gz \
sg_raw_data/unknown_reverse_unpaired.fq.gz > unknown_unpaired_library.fastq

bowtie2 --qc-filter -p 8 --local \
-x sg_reference/human_chr22 \
-U unknown_unpaired_library.fastq \
--un sg_raw_data/unknown_enriched.fastq > /dev/null
```
:::
:::

### Assemble the unknown genome

After eliminating the significant amount of human genome background from our raw data, we still only have the short next-generation sequencing reads for our unknown genome. In the next step we will use one of the most popular genome assembler [SPAdes](https://github.com/ablab/spades), the application that can be used for both single genome and multiple genomes (metagenome) assembly. After looking at the command line help, we will execute the `spades.py` script with very basic settings (input reads and output folder).

```bash
spades.py

# Run the assembly using all the default settings, only giving the input raw data file and the output folder
spades.py -s sg_raw_data/unknown_enriched.fastq -o results/unknown_genome

# The assembly results will be in the results/unknown_genome/ folder, check the output files, logs, warnings
cd results/unknown_genome/

```

:::{.callout-exercise}
#### Investigate the output of the assembly step

Go through the output files of the `spades.py` script together with the trainers. Look into the log files, discuss the files that are worth to keep long term. Investigate the level of fragmentation of the assembled novel genome.
:::

## Peforming de novo metagenomics assembly on the mixed community raw data

In this practical session we simulate a completely de novo assembly of shotgun sequencing data. We will use the same raw data we used for the alignment based metagenomics practical. As this is a synthetic community we know all the 20 genomes we mixed together in them together with their relative abundance. You can check the species names in the FASTA headers of the `sg_reference/mixed_bacterial_community_ncbi_genomes.fasta` file and their relative abundancies in the `sg_reference/mixed_bacterial_community_5M_abundance.txt` file.

### QC and pre-processing

The de novo assembly QC and pre-processing has the same first steps as any other pipeline on shotgun metagenomics data, as we did these steps already we will skip the `FastQC`, the `cutadapt` and the `trimmomatic` step and carry forward the files we generated during the [alignment based metagenomics practical](22-pract.html#standard-quality-control-and-pre-processing-of-shotgun-metagenomics-raw-data). While it is not crucial, but recommended to perform two extra pre-processing steps when we prepare the data for de novo metagenomcis assembly.

As the assembly step time and resource need is correlating significantly with the amount of input data, we can use methods to reduce the amount of raw reads without loosing important data. We will use the `clumpify.sh` script (from the `bbmap` package) to remove duplicates (PCR or optical). This algorithm removes completely matching reads or read-pairs.

```bash

# Be sure you have the metagenomics environment activated
# if not...
conda activate metagenomics

clumpify.sh

clumpify.sh in=mixedcomm_forward_paired.fq.gz in2=mixedcomm_reverse_paired.fq.gz \
out=mixedcomm_forward_paired_dedup.fq.gz out2=mixedcomm_reverse_paired_dedup.fq.gz \
dedupe=t
```
In the next step we will merge overlapping reads (like we did in the 16S pipeline).
The main aim of this is to reduce the amount of input sequence and to create longer merged reads that can help during the de novo assembly. We will use the `flash` application for this purpose.

```bash
flash --help

flash mixedcomm_forward_paired_dedup.fq.gz mixedcomm_reverse_paired_dedup.fq.gz
ls -ltr
```

The script generated 3 files, one for the merged reads (`out.extendedFrags.fastq`) and 1-1 for the non-merged paired-end reads (`out.notCombined_1.fastq` and `out.notCombined_2.fastq`). We will carry forward 3 files for the de novo assembly, two for the paired-end reads and one for unpaired and merged reads. Before we start the assembly 
we put all the unpaired reads (came from the `trimmomatic` step and the `flash` merge) into one file. Be aware that the merged file (from the `flash` step) is plain FASTQ file, while the unpaired files are gzip-ed, so needed to be unzipped.

```bash
zcat mixedcomm_forward_unpaired.fq.gz mixedcomm_reverse_unpaired.fq.gz >> out.extendedFrags.fastq
```

### De novo assembly step

The pre-processed data now ready to be fed into the de novo assembler algorithm. We are using the same application (`SPAdes`) as we used earlier for a single genome assembly, but for this purpose the software package provides a metagenome specific assembler with certain optimisations specific to metagenomics data source. First, look at the command line help of the algorithm, than execute the script with our input data. We only define the number of CPU cores to use together with the 3 input files and the output directory for the results.

```bash
# Perform the de novo assembly step
metaspades.py

metaspades.py -t 8 \
-1 out.notCombined_1.fastq \
-2 out.notCombined_2.fastq \
-s out.extendedFrags.fastq \
-o ../results/mixed_comm_de_novo &

```

:::{.callout-note}
Please note that we use the `&` symbol at the and of the command. This symbol tells the Unix / Linux system to put the run in the background and run the process even if we log out from the server. This is particularly helpful as the de novo assembly usually runs for long time. Even our simulated training data will run for about an hour, a real shotgun metagenomics data can easily provide the amount of raw reads to make the assembly step longer than a day even on a much higher spec computer (e.g., 56 or 67CPU cores).
:::