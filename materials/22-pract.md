---
title: Practical
---

## Finding a known genome in mixed microbial community

The following practical simulates the situation when we know what we are looking for and utilise a whole-genome alignment approach to find a particular bacterial strain or virus genome in a highly complex microbial community. This technique is especially useful if the known species (potential pathogen) is hard to detect / culture with traditional laboratory methods, amplification-based methods are not specific enough (e.g., often false positive due to a common species with highly similar genome), the laboratory method requires more time and/or more expensive.

::: {.callout-tip}
#### Key Points to consider

- The closer the genome you use to the known species / strain the higher sensitivity and specificity can be achieved. If you are monitoring an outbreak with this approach, the best results can be achieved if you can isolate the microbe at least once and perform a whole genome sequencing and de novo assembly on it.
- If you plan to quantify the tracked species, you can have a broad idea (relative abundance) by comparing the aligned read number to the total number of reads in the raw data. You can achieve much better quantification (both relative and absolute) if you spike in your sample (before DNA extraction) with a known bacteria or virus, using a well-defined amount. Ideally the spiked in species should be a distant species in terms of phylogeny and should have similarish genome size.
- If you don't have your own reference genome, try to find one in public databases that is potentially the closest to your geographical location but also a recent isolate.
:::

### Standard quality control and pre-processing of shotgun metagenomics raw data

Before we perform any analysis on the raw data it is important to perform the basic quality control checks and if needed certain pre-processing and filtering steps to ensure that we are working with high quality data. When you open a new terminal in the training environment, your working directory should be the `~/Course_Materials` folder. You can always check where you are in the filesystem by using the `pwd` command or just by checking your `bash` prompt.

Switch to the folder with the raw data, check the usage of the application and run `fastqc` on our first example raw data set (a pair of FASTQ files simulating a paired-end sequencing data).

```bash
cd sg_raw_data/
fastqc -h
fastqc mixed_bacterial_community_5M_R1.fastq.gz mixed_bacterial_community_5M_R2.fastq.gz
```
::: {.callout-exercise}
#### Inspect the output of the FastQC application

FastQC generatesgraphical output report in `.html` format. This is often placed (and archived) together with the raw data files, so the quality measures can be quickly checked in the case of future usage.

Open the Html files and go through the graphs, discuss what you see. For a future reference, and to see more examples (good and bad data), please visit the [FastQC website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
:::

The next standard pre-processing step is to remove adapter, primer and other unwanted sequences from your reads. These sequence contents are side products of the next-generation sequencing technique, and as they are not coming from the template DNA, they can interfere from many downstream pipeline steps. One of the most commonly used tool for this purpose is [cutadapt](https://cutadapt.readthedocs.io/en/stable/). Check the command line help for the application and run the filtering step on the raw sequencing data (please note, that you have to activate the `metagenomics` conda environment if it is not yet active).

```bash
conda activate metagenomics

cutadapt -h

cutadapt -a CTGTCTCTTATACACATCT -A ATGTGTATAAGAGACA \
-o mixed_bacterial_community_5M_noadapt_R1.fastq.gz -p mixed_bacterial_community_5M_noadapt_R2.fastq.gz \
mixed_bacterial_community_5M_R1.fastq.gz mixed_bacterial_community_5M_R2.fastq.gz

```

The final standard pre-processing step is to remove bad quality sequencing from our raw data. Typically the 5' and 3' ends of the sequencing reads have lover quality (so it is worth to remove these ends), and it is also common that a few reads are in general bad quality. We will use the [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) application to perform this quality filtering step, please refer to the on-line website or the command line help (`trimmomatic -h`) for more information on usage and fine-tuning.

```bash
trimmomatic PE -phred33 mixed_bacterial_community_5M_noadapt_R1.fastq.gz mixed_bacterial_community_5M_noadapt_R2.fastq.gz \
mixedcomm_forward_paired.fq.gz mixedcomm_forward_unpaired.fq.gz \
mixedcomm_reverse_paired.fq.gz mixedcomm_reverse_unpaired.fq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
::: {.callout-tip}
- Please note, that we are defining output file names for both the filtered paired-end and unpaired reads. While our input is coming from a paired-end library, it is possible that one mate of a read pair has bad quality while the other mate is good enough to carry forward. In this case the good read will be carried forward as unpaired read.
- To be able to use the `trimmomatic` tool efficiently, you have to have a good understanding of both the sequencing technique, the significance of potential errors and warnings listed in the `FastQC` report and the usage of the many different options and parameters in the `trimmomatic` application itself.
:::

::: {.callout-exercise}
#### Re-analyse the filtered raw data with FastQC and compare it to the original report
{{< level 2 >}}
Run the `fastqc` step again but this time giving the output files of the `trimmomatic` step to the application. Compare the reports and graphs and discuss the differences with the trainers.

::: {.callout-answer}
To run the `fastqc` on the `trimmomatic` output, you just simply list the new file names after the command. The new `.html` files will be generated next to the files:
```bash
fastqc mixedcomm_forward_paired.fq.gz \
mixedcomm_forward_unpaired.fq.gz \
mixedcomm_reverse_paired.fq.gz \
mixedcomm_reverse_unpaired.fq.gz
```
:::
:::

### Obtaing and preparing the reference genome for alignment

Reference genomes are commonly stored in [FASTA](21-pres.html#file-and-data-formats) formatted files. If the genome contains multiple DNA species (e.g., mutliple chromosomes or genome + plasmid) it is still stored in a single FASTA file, but in a multi-FASTA format, where each DNA species has its own FASTA header. These genome files can be generated (e.g., by bacterial colony sequencing and de novo assembly) or be downloaded from on-line genome databases. The most commonly used general databases are [Ensembl](https://www.ensembl.org/index.html), [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/), [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) and the metagenomics specific [MGnify](https://www.ebi.ac.uk/metagenomics).

There are multiple reference genome files saved in the `Course_Materials/sg_reference/` directory, for this exercise we will use the the reference genome of the _Shigella flexneri strain 2013C-3749_. We can use the `grep` command to list the FASTA headers in the file (to see if there are multiple DNA entities), and use the file size (in bytes) for a good estimation of the genome size.

```bash
cd ~/Course_Materials/sg_reference/

ls
grep ">" NZ_CP034931.fa
ls -l NZ_CP034931.fa
```

We will use one of the most popular short sequencing read aligner [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) to align the raw data to the _Shigella_ genome. The aligner requires the genome to be transformed to a binary format that makes the aligner searching significantly better and faster.

```bash
bowtie2-build -q NZ_CP034931.fa shigella_genome
ls
```

:::callout-note
Please note, that the `bowtie2-build` command resulted in six new files (with the `.bt2` extension). For the future alignment run, you will not need the original FASTA file, only these. If you are working with a large size genome (e.g., human genome), you can save space by removing the original FASTA file, in our case the original genome file is small, so it is not necessary to remove it.
:::

### Performing the reference based alignment

Our first command will run the `bowtie2` aligner with the default settings, we only need to give the path for the reference genome, the input files (we are only using the trimmed paired-end data) and the output file.

```bash
cd ..
# be sure, you are in the Course_Materials directory
pwd
mkdir results

bowtie2 -x sg_reference/shigella_genome \
-1 sg_raw_data/mixedcomm_forward_paired.fq.gz \
-2 sg_raw_data/mixedcomm_reverse_paired.fq.gz \
-S results/shigella_default_alignment.sam
```

::: {.callout-exercise}
#### Re-run the alignment step by utilising a few options
{{< level 3 >}}
Refer to the above linked website or the command line help (`bowtie2 -h`) to read about the numerous options the application provides. Re-run the analysis with the following options:
- use fast end-to-end alignment
- use all the available 8 CPU cores
- only output those reads in the final SAM file that have aligned to the reference genome
- be sure to use a different output file name

After the run, compare the alignment statistics, running time and the output file size to the previous run.

::: {.callout-answer}
```bash
bowtie2 -p 8 --fast --no-unal \
-x sg_reference/shigella_genome \
-1 sg_raw_data/mixedcomm_forward_paired.fq.gz \
-2 sg_raw_data/mixedcomm_reverse_paired.fq.gz \
-S results/shigella_fast_alignment.sam
```
:::
:::

Please also take your time to read about the different options `bowtie2` provides and discuss their potential usage with the trainers.

### Sort and compress the alignment results

As we previously mentioned, while the SAM format is text based and as a result human readable, it is not a good practice to keep you alignment data in SAM format as it can take a lot of storage space. Applications that are used to extract information from alignment data are all capable of reading in and processing the compressed binary file format for alignment data. Most of the downstream applications require the alignment data to be sorted (by alignment position), so we will do the sorting and SAM -> BAM conversion in one step.

```bash
cd results/

samtools sort -@ 4 -O BAM -o shigella_default_alignment_sorted.bam shigella_default_alignment.sam
samtools sort -@ 4 -O BAM -o shigella_fast_alignment_sorted.bam shigella_fast_alignment.sam

ls -lh
rm *.sam
```
:::callout-note
Please note the amount of reduction in file size after the SAM -> BAM conversion. Remember that the BAM format contains all the data that the SAM file contained, so the SAM files can be removed.
:::

The tool we used in this step is another good example for the 'Swiss army knife' type applications in bioinformatics. Just like the previously used `bowtie2`, `samtools` have a lot of options, running modes (commands) so it can be used for many different purposes. Please take your time to read about the functions on the [on-line manual page](http://www.htslib.org/doc/samtools.html) and discuss the potential usage with the trainers.

### Alignment pipeline using multiple reference genomes

While we can use a single genome as a reference and only aiming for detecting the presence / absence of a certain genome in our mixed microbial community, these tools have the potential to provide more information, e.g., relative or absolute abundance. If you are working with a completely unknown sample (e.g., human faecal, or environmental sample) so you don't know any specific species presence and abundance in your sample, you can 'spike-in' a well-defined amount of a lab strain during the sample preparation (ideally before DNA extraction) and compare your detected pathogen genome amounts to that. In the next exercise we are simulating a slightly different scenario, In this case, you are working with a synthetic community, you mix the different species together and would like to describe the relative abundances of the different species.

:::callout-note
Synthetic communities are often used in the laboratory for different purposes:

- Studying or optimising fermentation processes that are involving multiple species
- Simulating gut microbiome with a reduced microbial community and study the changes in response to various effects (antibiotics treatment, pathogen invasion, etc)
- Studying bacterial interactions, metabolic interplay

As these communities are put together in the lab, the composition in terms of 'members' is known, dut metagenomics analysis can be used to finely map the dynamic changes in the community composition (in terms of abundance).
:::

::: {.callout-exercise}
#### Align your data to multiple genomes
{{< level 3 >}}
Using the above single genome alignment as an example, prepare an alignment file for compositional profiling. You will use the same input data (the trimmomatic paired-end output files), and similar pipeline steps, but will use a different reference genome. You can find the reference genome collection (20 different genomes) in the `sg_reference/` directory with the file name `mixed_bacterial_community_ncbi_genomes.fasta`.

::: {.callout-answer}
```bash
cd sg_reference/

bowtie2-build -q mixed_bacterial_community_ncbi_genomes.fasta all_genomes

cd ..

bowtie2 -p 8 --fast --no-unal \
-x sg_reference/all_genomes \
-1 sg_raw_data/mixedcomm_forward_paired.fq.gz \
-2 sg_raw_data/mixedcomm_reverse_paired.fq.gz \
-S results/all_genomes_fast_alignment.sam

cd results/

samtools sort -@ 4 -O BAM -o all_genomes_fast_alignment_sorted.bam all_genomes_fast_alignment.sam

rm *.sam
```
:::
:::

### Alignment file post-processing

The BAM files you have generated are ready to use in multiple applications for further processing or visualisation. Most of these applications require the alignment to be in sorted format (be genome position), that is why we went through the `samtools sort` step. Certain applications may also require you to 'index' the BAM file (to make navigation within the BAM file faster), you can do it with a simple command:

```bash
samtools index bam_file
#E.g., in your case
samtools index shigella_fast_alignment_sorted.bam
```
This command will create a `.bai` file next to the BAM file, with the same base name, so applications will find it as the index of the BAM file.

We will extract some basic statistics and coverage information from our multiple genome aligned data using various `samtools` commands.

```bash
samtools flagstat all_genomes_fast_alignment_sorted.bam

samtools stats all_genomes_fast_alignment_sorted.bam > alignment_statistics.txt
less alignment_statistics.txt

samtools coverage all_genomes_fast_alignment_sorted.bam

# Save this coverage file for future comparison
samtools coverage all_genomes_fast_alignment_sorted.bam > mixed_community_coverage_from_alignment.txt
```
:::{.callout-tip}
BAM files can be used for various other purposes:

- Visualisation by [Integrated Genome Viewer / IGV](https://www.igv.org) or [BamView](https://www.sanger.ac.uk/tool/bamview/)
- Variation calling by using the combination of `samtools`, `bcftools` and `vcftools`
- Manipulating data in Python using the [pysam](https://pysam.readthedocs.io/en/latest/index.html) package
- Manipulating data in R by using the [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) and [GenomicAlignments](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html) packages
:::

## Compositional mapping of shutgun metagenomics data

The previously used whole-genome alignment based mapping method works very well if we know what we are looking for, e.g., we try to detect a known pathogen or follow the dynamic changes in a synthetic community. When we try to describe the composition of more complex communities the potential reference database would be so big using whole genomes (Bacteria + Archea + Fungi + Virus), so we apply different reduction strategies. The two main approaches are:(i) Finding biologically meaningful genetic signatures (e.g., genes that are only present in certain taxonomic group) and search for those in the raw data; (ii) extract numerous random short fragments from genomes and try to find exact matchings for those in the raw data.

### MetaPhlAn, a marker gene based profiling method

[MetaPhlAn](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4) is using a large set of marker genes (~5.1 million) to distinguish between different taxonomic groups (clades). The method can reliably identify genomes down to the species level, further resolution (down to strain level) is possible with the developer's other algorithm [StrainPhlAn](https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1). The application uses its own database of selected genes, and due to the large amount of information it uses, it requires a relatively powerful computer to run (at least 6-8 CPU cores and 32GB RAM is highly recommended).The installation and usage is very well documented on the GitHub site (link above).

We will do the basic profiling on our pre-processed raw data from the synthetic community and will inspect the output results.

```bash
metaphlan sg_raw_data/mixedcomm_forward_paired.fq.gz,sg_raw_data/mixedcomm_reverse_paired.fq.gz \
--bowtie2out results/metaphlan.bowtie2.bz2 --nproc 5 --input_type fastq -o results/profiled_metagenome.txt

less results/profiled_metagenome.txt

# let's save the species level taxonomic and abundance information for future comparison
grep "s__" results/profiled_metagenome.txt | grep -v "t__" > results/profiled_metagenome_sp_only.txt

```

### MASH, a random k-mer based profiling method

There are several bioinformatics applications based on the extraction of multiple k-mers from genomes and then use those to detect the presence of these genomes in raw sequencing or assembled data from mixed microbial communities. The two most popular methods are [MASH](https://mash.readthedocs.io/en/latest/) and [Kraken2](https://ccb.jhu.edu/software/kraken2/). K-mers are short (typically 14-30 nucleotides) DNA fragment that are randomly cut out from genomes. The number of k-mer per genome is between a few hundred and a few thousand depending on the required final sensitivity and specificity. If we take an average setting, where we extract 1000 random 21-mer fragments from an average bacterial genome (~5 million basepairs in length), we only cover less than 5% of the genome. While this is a significant reduction in term of information content, it can still be enough to identify species with good sensitivity and specificity. In the following practical, we will demonstrate the usage of `mash`, we will create our own databases of random k-mers (with different length and number of k-mers) and will also use some generic databases (originally created by using the NCBI RefSeq database).

First, let's check the different functions of `mash`, look at their usage and inspect the RefSeq databases.

```bash
# Exit from the metagenomics environment and go back to the base
conda deactivate

mash

mash info sg_reference/RefSeqSketches.msh | head -n 20
mash info sg_reference/RefSeqSketchesDefaults.msh | head -n 20

mash sketch -h
```

We will create a new sketch file from the genomes we used to build up our synthetic community. While this sketch database will have very limited usage, we use this small set of genomes to demonstrate the process of sketching (random k-mer extraction). Doing the same process on meaningful (and potentially very large) databases would take much more time. First we examine the input data, then create a small and 'light-weight' database. 

```bash

cd sg_reference/

less mixed_bacterial_community_ncbi_genomes.fasta
grep ">" mixed_bacterial_community_ncbi_genomes.fasta
grep ">" mixed_bacterial_community_ncbi_genomes.fasta | wc -l

mash sketch -i -s 500 -k 16 -o mixed_community mixed_bacterial_community_ncbi_genomes.fasta

# Compare the original fasta file size with the sketch size
ls -lh

mash info mixed_community.msh
```

::: {.callout-exercise}
#### Create a high resolution sketch database from the same input data
{{< level 2 >}}
Modifying the options of `mash sketch` by extracting __5000 k-mers__ with the __size of 21__ nucleotides, and save it in an output file named `mixed_community_hr` (referring to high resolution). Compare the size of the two different sketch databases and the original FASTA file, print out the basic information about the new Sketch database. Refer to `mash sketch -h` if you need information on the different options.

::: {.callout-answer}
```bash
mash sketch -i -s 5000 -k 21 -o mixed_community_hr mixed_bacterial_community_ncbi_genomes.fasta

ls -lh

mash info mixed_community_hr.msh
```
:::
:::

In the next part, you will learn how to profile a mixed community using `mash screen`. After looking at the general help for the command, you will test a scenario when you have either a full genome or just fragments of it. The second reference sequence (`NZ_CP038419_1.fasta`) you will screen with `mash` is not the full genome of the bacteria but only the collection of its genes. While bacterial genomes have low amount of non-coding DNA, the genome is still highly fragmented if each gene is separated to individual fasta records. Not exactly the same way but ending up with highly fragmented genomes in metagenomics studies is very common (especially in de novo assembly). It can happen, that the randomly extracted k-mer in the original genome was in a region that is split into two in the fragmenteg genome. In this case the search algorithm will not find a good match for that k-mer.

:::{.callout-note}
The output columns for `mash screen` are the following:

- identity: level of similarity between the query and the database reference sequence
- shared-hashes: number of matching hashes (k-mers) between the query and the database reference sequence
- median-multiplicity
- p-value of false detection
- query-ID of the database entry
- query-comment of the database entry
:::

```bash
mash screen -h

less NZ_CP034931.fa
grep ">" NZ_CP034931.fa

less NZ_CP038419_1.fasta
grep -c ">" _NZ_CP038419_1.fasta

mash screen mixed_community.msh NZ_CP034931.fa
mash screen mixed_community.msh NZ_CP038419_1.fasta
```

Notice, that in both cases not only the query genome came up with high number of shared hashes. This is caused by high similarity between certain genomes but can be avoided by using the `-w` option. Please read more about this option in the application's help `mash screen -h` or on the [application's website](https://mash.readthedocs.io/en/latest/tutorials.html#screening-a-read-set-for-containment-of-refseq-genomes). Let's try running the same commands with the `-w` option and compare the results

```bash
mash screen -w mixed_community.msh NZ_CP034931.fa
mash screen mixed_community.msh NZ_CP038419_1.fasta

mash screen -w mixed_community.msh NZ_CP038419_1.fasta
mash screen mixed_community_hr.msh NZ_CP038419_1.fasta
```
In our previous examples, we used a database that was generated from the genomes we mixed together in the synthetic community against a single genome that we knew was in the synthetic community. This is rarely the case in real life settings, let's try to screen our raw sequencing data with both our own database and the general RefSeq databases. For the bigger databases we can also try to speed up the process by using multiple CPU cores. Please note that we only print out the 50 best (based on identity level) hist from the big database screens, otherwise the terminal would be flooded with the results.

```bash
mash screen -w mixed_community_hr.msh \
../sg_raw_data/mixedcomm_forward_paired.fq.gz \
../sg_raw_data/mixedcomm_reverse_paired.fq.gz

mash screen -w RefSeqSketchesDefaults.msh \
../sg_raw_data/mixedcomm_forward_paired.fq.gz \
../sg_raw_data/mixedcomm_reverse_paired.fq.gz | sort -gr -k 1 | head -n 50

mash screen -w RefSeqSketches.msh \
../sg_raw_data/mixedcomm_forward_paired.fq.gz \
../sg_raw_data/mixedcomm_reverse_paired.fq.gz | sort -gr -k 1 | head -n 50

```