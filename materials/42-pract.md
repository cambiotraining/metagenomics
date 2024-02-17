---
title: Practical
---

```bash
########  Day 3 practical 2: De novo metagenomics practical continued  #######

# Let's onspect the output folder from the de novo assembly

cd results/mixed_comm_de_novo/
ls -lh
less contigs.fasta
grep -c ">" contigs.fasta
grep ">" contigs.fasta | head -n 50

# Let's use maxbin2 to bin the contigs

#### if not installed
micromamba deactivate
micromamba create -n maxbin2
micromamba activate maxbin2
micromamba install -c bioconda -c conda-forge maxbin2
#### end of install

# Let's see the application usage and options

run_MaxBin.pl
mkdir MAXBIN
run_MaxBin.pl -thread 8 -contig contigs.fasta -out MAXBIN/mixed_comm \
-reads ../../sg_raw_data/out.notCombined_1.fastq -reads2 ../../sg_raw_data/out.notCombined_2.fastq -reads3 ../../sg_raw_data/out.extendedFrags.fastq

# Let's inspect the output folder and some files
cd MAXBIN
ls -lh

less mixed_comm.abundance
less mixed_comm.summary 
less mixed_comm.abund1 
less mixed_comm.marker

# Check the quality of our binned assemblies
micromamba deactivate
micromamba activate metagenomics

checkm taxonomy_wf domain Bacteria -x fasta MAXBIN/ CHECKM/

# Inspect the output table and the files that were created in the CHECKM folder

# Create a more detailed quality assessment
checkm qa CHECKM/Bacteria.ms CHECKM/ --file CHECKM/quality_mixed_comm.tsv --tab_table -o 2

# Load in the created .tsv file into a spreadsheet application and discuss the content and significanc of individual columns

# Let's see if we can find out something about the taxonomic backgrount of our MAGs
gtdbtk classify_wf --out_dir GTDBTK/ --genome_dir MAXBIN/ -x fasta --skip_ani_screen --cpus 8
# Runs for ~25ins

# Inspect the output files
cd GTDBTK
less gtdbtk.bac120.summary.tsv

# Visualise one of the trees in iTol 
cd classify
gtdbtk convert_to_itol --input_tree gtdbtk.backbone.bac120.classify.tree --output_tree gtdbtk.backbone.bac120.classify.tree.itol

# Sadly these trees are super huge so really hard to visualise them in iTol

# Choose our best MAG and annotate it with Prokka
micromamba deactivate
micromamba create -n prokka
micromamba activate prokka
micromamba install -c conda-forge -c bioconda -c defaults prokka

prokka mixed_comm.003.fasta

# We can also use abricate to search for AMR and virulence genes
abricate mixed_comm.002.fasta
abricate mixed_comm.006.fasta


```

