# ----- Installation + Set up ----------------- # -----

library("dada2")
library("phyloseq")
library("Biostrings")
library("tidyverse")
library("readxl")
library("factoextra")
theme_set(theme_bw())

# ----- Input and organisation of Samples ----- # -----

# Give path to the samples
path <- "amplicon_16S"
list.files(path)

# Give forward and reverse files separately
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))

# Extract sample names from sample names by removing the `_1.fastq.gz` suffix.
sample.names <- str_remove(basename(fnFs), "_1.fastq.gz")

# ----- Sample pre-processing ----------------- # ------

# Inspect the quality of each sample in forward and reverse
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])

# In gray-scale is a heat map of the frequency of each quality score at each base position.
# The mean quality score at each position is shown by the green line, and the quartiles of the
# quality score distribution by the orange lines. The red line shows the scaled proportion of
# reads that extend to at least that position

# Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ sub-directory
filtFs <- file.path(
  path,
  "filtered",
  paste0(sample.names, "_filt_1.fastq.gz")
)
filtRs <- file.path(
  path,
  "filtered",
  paste0(sample.names, "_filt_2.fastq.gz")
)
names(filtFs) <- sample.names
names(filtRs) <- sample.names


# Using standard filtering parameters and the defined trimming based on the plot quality
out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE, # On Windows set multithread=FALSE
)
head(out)


# The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different set of error rates.
# The learnErrors method learns this error model from the data, by alternating estimation of the error rates and inference of
# sample composition until they converge on a jointly consistent solution.

# Learn Error Rates - takes ~5 minutes to complete
# The error rates for each possible transition (A→C, A→G, …) are shown. Points are the observed error rates for
# each consensus quality score. The black line shows the estimated error rates after convergence of the machine-learning algorithm.
# The red line shows the error rates expected under the nominal definition of the Q-score.
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)


# Sample Interference
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

#dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
#dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

#dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
#dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")

# Inspecting the returned dada-class object
dadaFs[[1]]
dadaRs[[1]]

### dada function processes each sample independently. However, pooling information across samples can increase sensitivity to
# sequence variants that may be present at very low frequencies in multiple samples. The dada2 package offers two types of pooling.
# dada(..., pool=TRUE) performs standard pooled processing, in which all samples are pooled together for sample inference.
# dada(..., pool="pseudo") performs pseudo-pooling, in which samples are processed independently after sharing information between
#                          samples, approximating pooled sample inference in linear time

# ----- Merge Paired reads -------------------- # -----

## Merging is performed by aligning the denoised forward reads with the
## reverse-complement of the corresponding denoised reverse reads,
## and then constructing the merged “contig” sequences

# By default, merged sequences are only output if the forward and reverse reads
# overlap by at least 12 bases, and are identical to each other in the
# overlap region (but these conditions can be changed via function arguments)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# contains the merged $sequence, its $abundance, and the indices of the $forward
# and $reverse sequence variants that were merged. Paired reads that did not exactly
# overlap were removed by mergePairs, further reducing spurious output.

# ----- Construct sequence table -------------- # -----

# Construct an amplicon sequence variant table (ASV) table,
# a higher-resolution version of the OTU table produced by traditional methods.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##   The sequence table is a matrix with rows corresponding to (and named by)
##   the samples, and columns corresponding to (and named by) the sequence variants.

##   This table contains x number of ASVs (when x is the 2nd number returned
## using the dim(seqtab) command ), and the lengths of our merged sequences
## all fall within the expected range for this V4 amplicon.

# ----- Remove chimeras ----------------------- # -----

# The core dada method corrects substitution and indel errors, but chimeras remain. Fortunately, the accuracy of sequence
# variants after denoising makes identifying chimeric ASVs simpler than when dealing with fuzzy OTUs. Chimeric sequences are
# identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant
# “parent” sequences.
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
dim(seqtab.nochim)

# gives % of reads kept after removing chimeras
(1 - (sum(seqtab.nochim) / sum(seqtab))) * 100
# ^ % lost

# ----- Track reads through the pipeline ------ # -----

# As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:

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

# table of what made it through each step

# visualise the number of reads through each step
track %>%
  # convert table to long format
  pivot_longer(
    cols = -sample,
    names_to = "metric",
    values_to = "number_reads"
  ) %>%
  # reorder the metrics by a more logical order
  mutate(
    metric = factor(
      metric,
      levels = c(
        "input",
        "filtered",
        "denoisedF",
        "denoisedR",
        "merged",
        "nonchim"
      )
    )
  ) %>%
  # calculate percentage of reads at each step
  group_by(sample) %>%
  mutate(pct = number_reads / max(number_reads) * 100) %>%
  # plot
  ggplot(aes(metric, pct, group = sample)) +
  geom_line() +
  xlab("Pipeline step")


# ----- Assign Taxonomy ----------------------- # -----

# The dada2 package also implements a method to make species level assignments
# based on exact matching between ASVs and sequenced reference strains

# assign taxonomy based on a downloaded reference database
taxa <- assignTaxonomy(
  seqtab.nochim,
  "amplicon_16S/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
  multithread = TRUE
)

# ----- Phyloseq ------------------------------ # -----

# get metadata
metadata <- read_excel(
  "amplicon_16S/Fecal_Sample_Collection.xlsx",
  sheet = "Metadata"
)
metadata <- metadata %>%
  # add "R" prefix to match IDs in our sequencing data
  mutate(Sample_ID = paste0("R", Sample_ID)) %>%
  # move Sample_ID as rownames
  column_to_rownames(var = "Sample_ID")

# make sure samples in metadata table are in the same order as the seqtab
metadata <- metadata[rownames(seqtab.nochim), ]

# confirm that sample names match between objects
all(rownames(metadata) == rownames(seqtab.nochim))

# make ps object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(metadata),
  tax_table(taxa)
)

View(sample_data(ps))

# alpha div
alpha_div <- microbiome::alpha(ps, index = c("Chao1", "Shannon"))
new_Metadata <- merge(sample_data(ps), alpha_div, by = 0, all.x = TRUE)
new_Metadata <- new_Metadata %>%
  remove_rownames() %>%
  column_to_rownames(var = "Row.names")
sample_data(ps) <- new_Metadata


#transform for relative abundances
ps_composition <- microbiome::transform(ps, "compositional")

#Shannon Plots vs your chosen variables from metadata table, once you've saved the diversity measures to a regular data frame
# take the newest sample data from ps with relative abundances + subsetted
alpha_plot_meta <- sample_data(ps_composition) %>%
  as.matrix() %>%
  as.data.frame()
alpha_plot_meta_long <- pivot_longer(
  alpha_plot_meta,
  cols = chao1:diversity_shannon,
  names_to = "Diversity_Measure",
  values_to = "Value"
)

#The above makes values = characters, so move them back to numeric
alpha_plot_meta_long$Fecal_sample_weight_mg <- as.numeric(
  alpha_plot_meta_long$Fecal_sample_weight_mg
)
alpha_plot_meta_long$Value <- as.numeric(alpha_plot_meta_long$Value)


# Adding column defining which alpha diversity is which, based on another column:
alpha_plot_meta_long <- alpha_plot_meta_long %>%
  mutate(
    Diversity_Measure_Group = case_when(
      Diversity_Measure %in% c("chao1", "diversity_shannon") ~ "Diversity"
    )
  )

# standard ggplot to browse results + create figures <
ggplot(
  alpha_plot_meta_long,
  aes(
    x = Facility,
    y = Value,
    colour = Fecal_sample_weight_mg,
    #shape = Sex,
    na.rm = T
  )
) +
  geom_point(position = position_jitterdodge(jitter.width = 0), size = 4) +
  facet_wrap(~Diversity_Measure, scales = 'free_y') +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )


# PCA
pca <- subset_samples(ps_composition, diversity_shannon <= 4.7) # remove the 3 outliers
pca <- tax_glom(ps_composition, taxrank = 'Family', NArm = T)

matrix_ps <- as.data.frame(as.matrix(pca@otu_table))
matrix_ps <- matrix_ps[, colSums(matrix_ps) != 0]
res.pca <- prcomp(matrix_ps, scale = TRUE)

meta <- as.data.frame(as.matrix(pca@sam_data))
meta$Facility <- as.factor(meta$Facility)
meta$Sex <- as.factor(meta$Sex)


fviz_pca_ind(
  res.pca,
  repel = TRUE,
  axes = c(1, 2),
  label = "ind",
  habillage = meta$Facility,
  invisible = "quali",
  addEllipses = FALSE,
  ellipse.type = "t",
  ellipse.level = 0.99, # non-normal distribution of bugs
  pointsize = 4,
  legend.title = "Cohorts",
  ggtheme = theme_minimal(),
  title = " "
) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  scale_colour_brewer(palette = "Dark2")

### alt = on DADA2 website: https://benjjneb.github.io/dada2/tutorial.html
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
plot_ordination(ps.prop, ord.nmds.bray, color = "Facility", title = "Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU / sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps, fill = "Phylum") +
  facet_grid(~Facility, scales = "free_x", space = 'free')
plot_bar(ps.top20, fill = "Class") +
  facet_grid(~Facility, scales = "free_x", space = 'free')
plot_bar(ps.top20, fill = "Order") +
  facet_grid(~Facility, scales = "free_x", space = 'free')


tmp <- tax_glom(ps_composition, taxrank = 'Family', NArm = F) # merges the blocks
plot_bar(tmp, fill = "Family") +
  facet_grid(~Facility, scales = "free_x", space = 'free')
