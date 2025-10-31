---
title: 16S sequencing
---

## Presentation on the 16S sequencing technique

The most commonly used amplicon based technique in metagenomics is the targeted sequencing of certain variable regions of the 16S ribosomal RNA gene in bacterial genomes. The presentation will give a short introduction to the technology, while the practical session will go through each step of the analysis process. During the practicals, you will learn about the analysis principles and will have deep understanding on the benefits and caveats of the 16S-sequencing technique.

The presentation is accessible on Google Slides on the [following link](https://docs.google.com/presentation/d/1Hel01_xsi1ZH-HbSw7jxEGvy4RDpaefFEykcaaygilc/edit?usp=sharing).

::: {.callout-tip}
#### Key points

- 16S amplicon sequencing is the quickest and cheapest method to get community composition information from complex microbial samples.
- Sequence context (primer recognition sequence and sequence composition) may introduce abundance bias to the analysis.
- Sensitivity is usually better than shotgun metagenomics
- Resolution can reach species level but rarely to strain level
- Analysis can be done in `R` / `RStudio` using gold standard pipelines and analysis tools, complemented by extensive reference databases.
:::



## 16S Sequencing practical

During the practical you will analyse the data coming from a pilot study that was aimed to compare two animal facilities (mice faecal samples). Using the results you will be able to investigate the within and between group variance and will be able to infer biologically relevant decisions. The laboratory protocol used PCR primers to amplify the V4 region of the 16S ribosomal RNA gene, the expected length of this region is about 252 base pairs.

### Preparing the environment

The base `R` installation is not providing the necessary functions for 16S amplicon sequencing and data analysis, we need to load in the required libraries for our analysis. The `dada2` package provides the functions for the data pre-processing and analysis, the `Biostrings` package helps in sequence data manipulation, the `phyloseq` package provides functions for post-processing (e.g., composition comparison) and finally we create our plots and graphs with the help of the `ggplot2` package. The last line sets the default plot theme to a simple black and white style.

```r
library("dada2")
library("phyloseq")
library("Biostrings")
library("tidyverse")
library("readxl")
theme_set(theme_bw())
```

### Data import and sample organisation

In the first line we define the path to our data folder, followed by reading the list of files for forward and reverse sequencing reads. Finally we extract the sample names from the forward file names.

```r
path <- "amplicon_16S"
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
sample.names <- str_remove(basename(fnFs), "_1.fastq.gz")
```

### Quality control and filtering

The first step during the 16S data pre-processing is the inspection of the raw sequencing read qualities and perform trimming and filtering accordingly. The `plotQualityProfile()`
function (from the `dada2` package) provides an overview of the read qualities in single or multiple input file(s).

```r
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])
```

::: {.callout-exercise}
#### Plot 3 random input files
{{< level 2 >}}
The provided code will always read in the first three files (as were listed alphabetically in the data folder), selecting the first 3 elements of the `fnFs` and `fnRs` vector. Modify the code by utilising the `sample()` function to select 3 random input file names from the same two vectors. This way you are not constrain your inspection to the same three files, if you re-run the command, you can check a diiferent set of 3 random samples.

::: {.callout-answer}
The `sample()` function needs two arguments, a vector (in our case the `fnFs` and `fnRs` variables with file names) and the number of random sampled elements.
```r
plotQualityProfile(sample(fnFs, 3))
```
:::
:::

::: {.callout-tip}
Use the `?plotQualityProfile` command in the R console to read more about the function. You will find the detailed description of the output plot in the "Details" section of the help page.
:::

The next commands will define new file names for the filtered raw data and finally perform the data filtering itself. Please use the `?filterAndTrim` command to get more information on the filtering function arguments.

```r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt_1.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_filt_2.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     maxN=0, maxEE=c(2,2), truncQ=15, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
```

::: {.callout-exercise}
#### Set the truncation length
{{< level 2 >}}
The `truncLen` argument of the function defines the lengths (first number for forward, second number for reverse reads) we plan to keep from raw reads. The length has to be set by checking the previously generated quality plots by aiming to keep only high quality regions.

::: {.callout-answer}
The forward reads are usually high quality along the full sequence, except the last few nucleotides, so we can keep almost the full length (usually we cut 10 nucleotides from the end). The reverse sequencing reads very often have a point where the quality distribution crashes (common in illumina sequencing), you need to cut more off from the end. As an example:
```r
filterAndTrim(..., truncLen=c(220,160), ...)
```
:::
:::

You can take a quick look at the amount of filtering has been done. If you notice significant drop in the number of reads carried forward, you may want to check the QC graphs and parameters given to the `filterAndTrim()` function.

```r
head(out)
```

### Error correction and read merging

The DADA2 package uses an adaptive error correction algorithm to learns the error rates introduced by PCR and sequencing. In default DADA2 learns the error rate for each samples individually but to increase sensitivity, pooling and pseudo-pooling of samples is also possible during this step. More information on the pooling techniques can be found on the [DADA2 website](https://benjjneb.github.io/dada2/pool.html).

```r
# this takes ~5 minutes to complete
errF <- learnErrors(filtFs, multithread=TRUE)  
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

In the next step we are merging the filtered and denoised forward and reverse reads to construct the full length amplicons. The function merges forward and reverse reads if those are overlapping by at least 12 bases, this length can be fine-tuned using the function's arguments (please refer to the function's help page `?mergePairs`).

```r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

### Creating and cleaning the ASV table

The final step in the DADA2 pipeline is to create the amplicon sequence variant (ASV) table that can be taken forward to more sophisticates analysis. In this table each sequence variant (V4 region with unique sequence) will be counted in each analysed sample, resulting a `n x m` dimensional table where `n` is the number of samples and `m` is the number of observed unique sequence variants.

```r
seqtab <- makeSequenceTable(mergers)

# Inspect the dimensions of the sequence table
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

::: {.callout-tip}
The sequence length table shows the lengths as headers and the counts for these lengths in the first row. Ideally most of your ASVs should be in the expected range (in the case of the V4 amplicon it is around 252bp), if you see a significant shift to both larger or smaller fragments, you may want to revisit your pre-processing steps.
:::

The ASV table need a final cleanup before it is ready to taken forward to compositional / phylogenic analysis. During the PCR amplification extended primer fragments where the polymerase enzyme didn't finish the extension on the full length sequence could re-anneal to a different template and result in chimeric sequence variants (the left and the right side of the sequence is coming from different templates). Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.

```r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#Inspect the dimensions of the final ASV table
dim(seqtab.nochim)

#Inspect the ratio of merged reads that passed the chimera filtering step
sum(seqtab.nochim)/sum(seqtab)
```

::: {.callout-tip}
The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on on factors including experimental procedures and sample complexity. While it is common that substantial amount of unique ASVs are removed (20-30%), these are usually low abundance ASVs and rarely add up to more than a few percent in terms of the amount of removed merged reads. 
:::

### Post-pipeline quality control

Our original raw data went through several pipeline steps where we applied some type of filtering, each time resulting in potential reduction in the number of data points (reads, ASVs). It is a good practice to inspect these reductions across all samples at the end of the pipeline. If we find a significant drop in numbers from one step to the other, that either highlights a general problem with our raw data or wrong choice of parameters (when we see the same drop for all samples) or signals for bad quality on just certain samples (when we only see drops only in a few samples). In both cases it is important to go back either to our laboratory notes or our pipeline and check the steps / samples.

```r
getN <- function(x) sum(getUniques(x))
track <- data.frame(
  sample = sample.names,
  input = out[, 1],
  filtered = out[, 2],
  denoisedF = sapply(dadaFs, getN), 
  denoisedR = sapply(dadaRs, getN), 
  merged = sapply(mergers, getN), 
  nonchim = rowSums(seqtab.nochim)
)
head(track)
```

The tabular data can also be visualised by `ggplot`.

```r
track %>% 
  pivot_longer(cols = -sample, names_to = "metric", values_to = "number_reads") %>% 
  mutate(metric = factor(metric, levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>% 
  ggplot(aes(metric, number_reads, group = sample)) + 
  geom_line() + 
  xlab("Pipeline step")
```

### Assign taxonomy

The discovered 16S V4 fragment sequences are still in plain sequence format in our dataset, we need to identify the ASVs and assign proper taxonomy to those. We are using the [Silva](https://www.arb-silva.de) database for this purpose.

```r
taxa <- assignTaxonomy(seqtab.nochim, "amplicon_16S/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
```


### Visualising results with Phyloseq

We will use the [Phyloseq](https://joey711.github.io/phyloseq/) package to perform further post processing and visualising differences between the faecal microbiomes from the two animal facilities. 

We are loading in the sample metadata from an external excel file and will use the `Sample_ID` column to match our DADA2 results with the metadata.

```r
metadata <- read_excel("amplicon_16S/Fecal_Sample_Collection.xlsx", 
                       sheet = "Metadata")
metadata <- metadata %>% 
  # add "R" prefix to match IDs in our sequencing data
  mutate(Sample_ID = paste0("R", Sample_ID)) %>% 
  # move Sample_ID as rownames
  column_to_rownames(var = "Sample_ID")

# make sure samples in metadata table are in the same order as the seqtab
metadata <- metadata[rownames(seqtab.nochim), ]

# confirm that sample names match between objects
all(rownames(metadata) == rownames(seqtab.nochim))
```

In the next step we will create the `phyloseq` object, this is a complex variable in R containing all the information needed for further analysis and visualisation. To construct the phyloseq object we need 3 input data: (i) the read count table of the ASVs (abundance data on different species in different samples); (ii) sample metadata (e.g., groups, treatments, disease status); (iii) taxonomic assignments for all ASVs. After creating the object we can use built-in commands to examine the content (e.g., the sample data table).

```r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

View(sample_data(ps))
```

We can use the created `phyloseq` object for many different analysis and visualisations, here we will demonstrate two often used exploration tool, dimensional reduction to see the main microbiome composition structure in all of our samples, and compositional plots for different taxonomic levels.

```r
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Facility", title="Bray NMDS")
```
The NMDS (Non-metric Multidimensional Scaling) method reduces the high dimensional complex data (each ASV is a different dimension) to two dimensions. This plot can be used to see the overall similarity / difference between samples, and especially useful to explore within and between groups varibility.

:::{.callout-exercise}
#### Investigate the plot

- Discuss the results with the trainers, try to find biological reasons for the difference between the two facilities.
- Discuss if certain animal experiments would fit better to Facility_1 or Facility_2.
- Explore the distribution of other metadata categories by replacing the `Facility` word for other column headers in `View(sample_data(ps))`.
- You can also try other methods for ordination, refer to the help page of the `ordinate()` function to see the list of available methods.
:::

Finally we can plot the composition (on different taxonomical levels) of each sample using bar plots.

```r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps, fill="Phylum") + facet_grid(~Facility, scales="free_x", space = 'free')
plot_bar(ps.top20, fill="Class") + facet_grid(~Facility, scales="free_x", space = 'free')
plot_bar(ps.top20, fill="Order") + facet_grid(~Facility, scales="free_x", space = 'free')
```
