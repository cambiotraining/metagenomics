---
title: De novo metagenomics assembly
---

## Peforming de novo metagenomics assembly on the mixed community raw data

In this practical session we simulate a completely de novo assembly of shotgun sequencing data. We will use the same raw data we used for the alignment based metagenomics practical. As this is a synthetic community we know all the 20 genomes we mixed together in them together with their relative abundance. You can check the species names in the FASTA headers of the `sg_reference/mixed_bacterial_community_ncbi_genomes.fasta` file and their relative abundancies in the `sg_reference/mixed_bacterial_community_5M_abundance.txt` file.

### QC and pre-processing

The de novo assembly QC and pre-processing has the same first steps as any other pipeline on shotgun metagenomics data, as we did these steps already we will skip the `FastQC`, the `cutadapt` and the `trimmomatic` step and carry forward the files we generated during the [alignment based metagenomics practical](22-pract.html#standard-quality-control-and-pre-processing-of-shotgun-metagenomics-raw-data). While it is not crucial, but recommended to perform two extra pre-processing steps when we prepare the data for de novo metagenomcis assembly.

As the assembly step time and resource need is correlating significantly with the amount of input data, we can use methods to reduce the amount of raw reads without loosing important data. We will use the `clumpify.sh` script (from the `bbmap` package) to remove duplicates (PCR or optical). This algorithm removes completely matching reads or read-pairs.

```bash
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