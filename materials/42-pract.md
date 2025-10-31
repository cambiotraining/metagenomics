---
title: Working with MAGs
---

## Binning, assigning quality and taxonomy to MAGs, Annotation

During this practical session we will bin the contigs to multiple clusters forming metagenome assembled genomes (MAGs). These MAGs ideally represent individual genomes, we will assign quality measures for these and do de novo gene annotation and specific gene discovery.

:::{.callout-important}
#### Activate your software environment

For this practical we need to activate the software environment called `mags`:

```bash
mamba activate mags
```
:::


### Reference independent binning

We will use the `maxbin2` algorithm to reconstruct the individual genomes of those bacteria that were sequenced in the artificial mixed community. The method uses genomic properties and contig coverage values to find clusters of contigs that are potentially coming from the same source (same genome). 

:::{.callout-note}
Due to the time limitations on this training course we only analyse a single sample, `maxbi2` performs significantly better if the binning step is performed on multiple samples in the same time. The multiple samples should come from different sources (e.g., faecal microbiomes from multiple patients), but not from very different sources (e.g., not recommended to bin gut and skin microbiome samples together).
:::

First let's just inspect the results of the overnight de novo assembly run.

```bash
cd ~/Course_Materials/

cd results/mixed_comm_de_novo/
ls -lh
less contigs.fasta
grep -c ">" contigs.fasta
grep ">" contigs.fasta | head -n 50
```
:::{.callout-note}
Please note, that while this is a small simulated dataset with only 20 different species mixed together, we still end up with very high number of contigs.
:::

Let's see the application options, create the output folder and run the binning step by giving the input `contigs.fasta` file, the output directory name, and the raw sequencing data files. The raw data files will be re-aligned to the contigs to calculate precise coverage values.

```bash
run_MaxBin.pl
mkdir MAXBIN

run_MaxBin.pl -thread 8 -contig contigs.fasta -out MAXBIN/mixed_comm \
-reads ../../sg_raw_data/out.notCombined_1.fastq \
-reads2 ../../sg_raw_data/out.notCombined_2.fastq \
-reads3 ../../sg_raw_data/out.extendedFrags.fastq
```
:::{.callout-exercise}
#### Inspect the output folder and result files of the `maxbin2` run
{{< level 2 >}}
Switch to the result directory, check the file sizes and look into the different application result files (not the `.fasta` files). Discuss the results with the trainers.

::: {.callout-answer}
```bash
cd MAXBIN
ls -lh

less mixed_comm.abundance
less mixed_comm.summary 
less mixed_comm.abund1 
less mixed_comm.marker
```
:::
:::

### Assigning quality and taxonomy to the MAGs

We will use two different applications (developed by the [same research group](https://github.com/Ecogenomics)) `checkm` and `gtdbtk` to assign taxonomy and quality to the MAGs. Both methods are using single-copy gene sets to assign quality and their built-in database to assign taxonomy. CheckM is the older but accepted "gold standard" method to assign quality to MAGs. As it was developed using a smaller set of genomes, it is significantly lighter (both in terms of storage space, time and resource need) compared to GTDB-tk.

We first run the taxonomy workflow of `checkm`, followed by the detailed quality assessment.

```bash
checkm taxonomy_wf domain Bacteria -x fasta MAXBIN/ CHECKM/

checkm qa CHECKM/Bacteria.ms CHECKM/ --file CHECKM/quality_mixed_comm.tsv --tab_table -o 2
```
:::{.callout-exercise}
#### CheckM result
Inspect the result files in the `CHECKM/` directory and discuss the assigned MAG quality measures with the trainers. Refer to the [software documentation](https://github.com/Ecogenomics/CheckM/wiki) for more information on the usage and output information.

You can load in the result `.tsv` file into a spreadsheet application to make it easier to investigate.
:::

The next step is to run the `gtdbtk` on the same MAGs. As it was mentioned before, this application uses a much larger reference database, so it also requires more time to run. The next command will run for ~25mins, so it is probably a good time to have a coffee and some biscuits. 

```bash

gtdbtk classify_wf --out_dir GTDBTK/ --genome_dir MAXBIN/ -x fasta --skip_ani_screen --cpus 8
```
:::{.callout-exercise}
#### GTDB-tk result
Inspect the result files in the `GTDBTK/` directory and discuss the assigned MAG quality measures with the trainers. Refer to the [software documentation](https://ecogenomics.github.io/GTDBTk/) for more information on the usage and output information.

You can load in the result `.tsv` file into a spreadsheet application to make it easier to investigate.
:::

### De novo gene annotation and screening for AMR genes

We will use the `prokka` application to annotate genes in one of our MAGs. While we will use the most basic way of running the script, `prokka` has a [lot of options](https://github.com/tseemann/prokka) for customising the process and results.

```bash
prokka mixed_comm.003.fasta
```
In default `prokka` outputs the results (annotation files, translated sequences, log files) in a directory named `PROKKA_yyyymmdd` (with the today's date). List an check the output files that were generated from the input MAG.

In our last example we will use the tool `abricate` and its built in databases to screen for antimicrobial resistance (AMR) genes. Abricate is a versatile tool that can also search for putative plasmid sequences, virulence genes and AMR genes using different databases. Please refer to the [github website](https://github.com/tseemann/abricate) to read about all the options, built-in databases and custom database creation.

```bash
abricate mixed_comm.002.fasta
abricate mixed_comm.006.fasta
```

