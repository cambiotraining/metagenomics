---
title: Presentation
---

## Working with metagenome assembled genomes (MAGs)

After performing de novo assembly and binning on a shotgun metagenomics data, the results are individual FASTA files with sets of contigs, representing a potential genome. There are multiple inherent errors in binning that is not always easy to avoid or repaer, but at least it is important to discover.

::: {.callout-tip}
#### Potential errors and problems with MAGs

- **Incomplete genome** - most often coming from low abundance and sequence coverage. Due to the lack of appropriate raw sequence coverage the full genome cannot be reconstructed. Very low abundance often leads to no binning at all and completely missing MAG for that species.
- **Contaminated genome** - The MAG contains genomic information from more than one origin (species). The most common cause is the presence of phylogenetically close species with comparable abundance.
- **Misplaced mobile genetic elements** - traditional shotgun metagenomics sequencing separates the genomic DNA from plasmids at the lysis step, even if the plasmids are perfectly assembled, it is highly unreliable to predict which host they were originally from. Hi-C technique can provide a solution for this problem.
:::

The slides for this presentation section can be found on google slides on the [following link](https://docs.google.com/presentation/d/1_b3ufd5Uk11vlefH6ifwy7svAPUgHomXbhXTgYYIi9w/edit?usp=sharing)

### Quality control, quality comparison

There are multiple solutions to assess the quality of the MAGs after binning. These tools can be used as general QC measures, but can also provide a good platform to compare different metagenomics approaches (e.g., long vs short read sequencing, different DNA extraction methods, different amount of sequencing raw data). While basic information (genome size, number of contigs, longest contigs, coverage info, etc.) can be extracted from the final MAGs byt using simple command line commands, more sophisticated measures can be obtained in two fundamentally different ways:

- **Reference genome based quality assessment**. In this case we use reference genome databases and try to find the matching reference genomes for our MAGs. After finding the closest reference genome, we can assess the identity and coverage percentage between the genome and the MAG, find missing or surplus genetic material. While this method is highly reliable and have very good specificity for known genomes, it performs poorly in the assessment of novel species. One of the most popular method for reference based assessment is the `QUAST` package (including `MetaQUAST`).
- **Single copy gene based quality assessment**. Thanks to our existing knowledge on hundreds of thousands prokaryotic genomes, researchers were able to find sets of genes (separately for bacteria and archea) that are existing in all species and existing only in one copy in the genome. While these gene sets are only having around 100 genes, the presence of these genes provide a good proxy to predict the MAG completeness and contamination. If we find all the marker genes in a MAG and we find it only in one copy, we can assume that all the other genes are present in similarly high ratio, so we can assess high completeness and low contamination for the MAG. If we detect less marker genes in the MAG, we can assume a proportionate incompleteness for the whole genome. If we find multiple copies in our MAG for multiple marker genes (that should only be present in single copy) we can assume that the MAG is contaminated (multiple species were binned into one MAG). The most popular applications for single-copy gene based assessment are `CheckM` and `GTDB-tk`.

:::{.callout-important}
It is not recommended to take contaminated MAGs forward to further analysis (annotation, taxonomic identification, gene search) as the genomic information is likely originated from multiple sources. Incomplete genomes are suitable for further analysis, but the absence of certain genes can be the consequence of the incompleteness.
:::

### Identification of genomes and certain gene groups

Each individual multi-FASTA file (MAG) can be screened in multiple ways to obtain information about its potential taxonomy, specific gene content, or even strain level information.

- We can use the k-mer based methods we used for raw data screening (e.g., `mash`), those methods work on both raw reads and contigs.
- We can simply use the NCBI's `BLASTn` algorithm, installed locally on our computer or through the [web surface](https://blast.ncbi.nlm.nih.gov).
- For well-known species (e.g., pathogene bacteria) we can use specialised tools to identify the MAG up to strain level (`mlst` application).
- To identify plasmid sequences, we can either use reference based methods (e.g., `ABRicate`) or sequence composition based methods (e.g., `PlasFlow`, `PlasForest`)
- We can find species or genus specific tools to discover more about the MAG (e.g., `legsta` for *Legionella pneumophila*)
- The [GitHub website](https://github.com/tseemann) of Torsten Seemann is an excellent resource for applications in this category.

### De novo genome annotation

Even if we were able to find out the exact species information for a MAG, but especially if we have find novel genomes (cannot really find a good matching reference genome for it), we can perform de novo gene annotation. During this process we screen for open reading frames (ORFs) and annotate those that are looking like genes. Often, simply based on nucleitide or amino acid homology, we can also annotate the function and the name for the predicted gene. De novo gene annotation tools are generally reliable for archeal and bacterial genomes, it is much more challenging to annotate eukaryotic genomes. The two most popular de novo annotation software for prokaryotic genomes are `prokka` and `bakta`. De novo annotation pipelines result in specific file types (typically GFF3 and/or GeneBank format) containing the names, positions, DNA and translated protein sequence of the annotated genes. These files can be taken forward to do comparisons and visualisations.

### Comparison to other genomes

Nucleotide level whole-genome comparisons can be performed between a MAG and a known genome, or between MAGs (e.g., from different samples or different time points). These comparisons provide high resolution information, even single nucleotide changes can be tracked (e.g., during an outbreak). While we often use known reference genomes for comparisons (e.g., downloaded from NCBI GeneBank or RefSeq), using the laboratory's own isolates or earlier genome assemblies from the same species often provide higher similarity. These comparisons can be performed with basic tools like `BLASTn` or high resolution k-mer (high k-mer count and long k-mers) based methods (e.g., using `mash`).

We often like to see our novel MAGs in a broader context. Comparing a genome with multiple isolates/strains from the same species or even with multiple species from the same genus can provide additional and highly valuable information (pan-genome analysis). In this cases the input data is often not the FASTA file with only nucleotide data but the annotated genome(s). After comparing multiple (even hundreds of) genomes we can identify genes that are present in all isolates / MAGs representing the "core genome" in the comparison, and genes that are isolate / MAG / subgroup specific representing the "accessory genome". These categories can provide information on both isolate / species / genus specific functions and properties (e.g., metabolic capabilities, virulence, etc.). The multi genome comparisons are usually time and resource intensive steps, the two popular applications to perform it are `roary` and `panaroo`.

:::{.callout-tip}
While using specific options a good level of compatibility can be achieved between de novo annotation and multiple comparison tools, in general `prokka` works the best in combination with `roary` while `bakta` works better with `panaroo`.
:::

### Visualisation tools

There is a continuously growing list of visualisation tools for multiple parts of the de novo metagenomics assembly post processing. A few examples from the wide-range of visualisation tools:

- Genome annotation visualisation: `Artemis`, `Icarus`
- Visualise pan-genome analysis: `Phandangoo`, `PanExplorer`, `PanACEA`
- Integrated analysis for multiple purposes: `MEGAN`
- Phylogenic tree visualisation and annotation: `iTOL`, R packages
- Quality control visualisation: `MultiQC`
- MAG quality visualisation: `QUAST`, R/ggplot2