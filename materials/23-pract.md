---
title: Finding an unknown genome
---

## Discover an unknown genome in low complexity mixed community

In this practical we are simulating a virus infection that is caused by a completely unknown virus. In this simulated dataset a novel genome is hidden in the raw data with a huge background of human genome sequence. This simulates a real scenario when a pathogen is in the blood or other (otherwise sterile) body fluid. To make the bioinformatics step faster we generated the human "background" from chromosome 22, so the database will be relatively small.

:::{.callout-important}
#### Activate your software environment

For this practical we need to activate the software environment called `assembly`:

```bash
mamba activate assembly
```
:::


### QC and Pre-processing

The raw date quality control and pre-processing is going the same way as we did with the mixed community data, for the details on these steps, please refer to [the appropriate practical material](22-pract.html#standard-quality-control-and-pre-processing-of-shotgun-metagenomics-raw-data).

```bash
cd sg_raw_data/
fastqc unknown_pathogen_R1.fastq unknown_pathogen_R2.fastq

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

First let's try to align our raw data to the known genome content (in this case the human chromosome 22 sequence). This is again going similarly to our previous practical, first we create a Bowtie2 database, then perform the alignment.

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
