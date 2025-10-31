---
title: Whole genome alignment based methods
---

## Finding a known genome in mixed microbial community

The following practical simulates the situation when we know what we are looking for and utilise a whole-genome alignment approach to find a particular bacterial strain or virus genome in a highly complex microbial community. This technique is especially useful if the known species (potential pathogen) is hard to detect / culture with traditional laboratory methods, amplification-based methods are not specific enough (e.g., often false positive due to a common species with highly similar genome), the laboratory method requires more time and/or more expensive.

::: {.callout-tip}
#### Key Points to consider

- The closer the genome you use to the known species / strain the higher sensitivity and specificity can be achieved. If you are monitoring an outbreak with this approach, the best results can be achieved if you can isolate the microbe at least once and perform a whole genome sequencing and de novo assembly on it.
- If you plan to quantify the tracked species, you can have a broad idea (relative abundance) by comparing the aligned read number to the total number of reads in the raw data. You can achieve much better quantification (both relative and absolute) if you spike in your sample (before DNA extraction) with a known bacteria or virus, using a well-defined amount. Ideally the spiked in species should be a distant species in terms of phylogeny and should have similarish genome size.
- If you don't have your own reference genome, try to find one in public databases that is potentially the closest to your geographical location but also a recent isolate.
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

## Alignment pipeline using multiple reference genomes

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

