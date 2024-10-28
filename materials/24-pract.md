---
title: Alignment based compositional profiling
---

## Compositional mapping of shotgun metagenomics data

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