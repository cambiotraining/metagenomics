library(ggplot2)
library(stringr)
library(gridExtra)

fasta_headers <- read.delim("~/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents/Cam_mg_course/fasta_headers.txt", header=FALSE)
fasta_headers$V1 <- str_replace(fasta_headers$V1, ">", "")

abundance_table <- as.data.frame(str_split_fixed(fasta_headers$V1, " ", 4)[,1:3])
colnames(abundance_table) <- c("ID", "Genus", "Species")
abundance_table$fullname <- str_c(abundance_table$Genus, abundance_table$Species, sep = " ")

mixed_bacterial_community_5M_abundance <- read.delim("~/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents/Cam_mg_course/mixed_bacterial_community_5M_abundance.txt", header=FALSE)
colnames(mixed_bacterial_community_5M_abundance) <- c("ID", "Orig_abundance")

abundance_table <- merge(abundance_table, mixed_bacterial_community_5M_abundance, by = "ID")

mixed_community_coverage_from_alignment <- read.delim("~/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents/Cam_mg_course/mixed_community_coverage_from_alignment.txt")
colnames(mixed_community_coverage_from_alignment)[1] <- "ID"
mixed_community_coverage_from_alignment <- mixed_community_coverage_from_alignment[,c("ID", "meandepth")]

abundance_table <- merge(abundance_table, mixed_community_coverage_from_alignment, by = "ID")
abundance_table$Align_abundance <- abundance_table$meandepth / sum(abundance_table$meandepth)

abundance_table_to_plot <- melt(abundance_table[,c("fullname", "Orig_abundance", "Align_abundance")])
colnames(abundance_table_to_plot) <- c("Species", "Abundance_source", "Relative_abundance")
p1 <- ggplot(abundance_table_to_plot, aes(Relative_abundance, Species, fill = Abundance_source)) + geom_col(position = "dodge")

profiled_metagenome_sp_only <- read.delim("~/Library/CloudStorage/OneDrive-UniversityofCambridge/Documents/Cam_mg_course/profiled_metagenome_merged_sp_only.txt", header=FALSE)
profiled_metagenome_sp_only <- profiled_metagenome_sp_only[profiled_metagenome_sp_only$V4 == "", 1:3]

profiled_metagenome_sp_only$V1 <- str_split_i(profiled_metagenome_sp_only$V1, "\\|", 7)
profiled_metagenome_sp_only$V1 <- str_replace(profiled_metagenome_sp_only$V1, "s__", "")
profiled_metagenome_sp_only$V1 <- str_replace(profiled_metagenome_sp_only$V1, "_", " ")

colnames(profiled_metagenome_sp_only) <- c("Species", "Codes", "Abundance")
profiled_metagenome_sp_only$Rel_abundance <- profiled_metagenome_sp_only$Abundance / sum(profiled_metagenome_sp_only$Abundance)

p2 <- ggplot(profiled_metagenome_sp_only, aes(Rel_abundance, Species)) + geom_col()

ggarrange(p1, p2, ncol = 2)
