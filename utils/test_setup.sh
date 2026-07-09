#!/usr/bin/env bash

# Stage dependencies in the steps below:
# Stage 16  ←  Stage 15  ←  Stage 14  ←  (raw input files)
# Stage 17  ←  Stage 7   ←  Stage 6   ←  (raw input files)
# Stage 18  ←  Stage 17
# Stages 19–23  ←  Stage 18 (metaspades contigs)

set -euxo pipefail

# Mamba environment names
ALIGN_ENV="alignment"
ASSEMBLY_ENV="assembly"
MAGS_ENV="mags"

# Adapter sequences used throughout the course
ADAPTER_FWD="CTGTCTCTTATACACATCT"
ADAPTER_REV="ATGTGTATAAGAGACA"

# Input data paths (relative to Course_Materials directory)
RAW_R1="sg_raw_data/mixed_bacterial_community_5M_R1.fastq.gz"
RAW_R2="sg_raw_data/mixed_bacterial_community_5M_R2.fastq.gz"
REFERENCE_SHIGELLA="sg_reference/NZ_CP034931.fa"
REFERENCE_ALL="sg_reference/mixed_bacterial_community_ncbi_genomes.fasta"
REFERENCE_HUMAN_CHR22="sg_reference/Homo_sapiens.GRCh38.dna.chromosome.22.fa"
MASH_DB_DEFAULT="sg_reference/RefSeqSketchesDefaults.msh"
UNKNOWN_R1="sg_raw_data/unknown_pathogen_R1.fastq"
UNKNOWN_R2="sg_raw_data/unknown_pathogen_R2.fastq"

mkdir -p results
WORKDIR="$(mktemp -d "results/training_machine_test.XXXXXX")"

cleanup() {
  local status=$?
  if [[ ${status} -eq 0 ]]; then
    rm -rf "${WORKDIR}"
  else
    echo "Test failed; keeping artifacts at: ${WORKDIR}"
  fi
}
trap cleanup EXIT

run_in_alignment() {
  mamba run -n "${ALIGN_ENV}" "$@"
}

run_in_assembly() {
  mamba run -n "${ASSEMBLY_ENV}" "$@"
}

run_in_mags() {
  mamba run -n "${MAGS_ENV}" "$@"
}

# Run a command with a timeout; succeed if the command exits 0 or is killed
# by the timeout (exit 124); fail on any other non-zero exit code
run_with_timeout() {
  local secs=$1
  shift
  local ec=0
  timeout "${secs}s" "$@" || ec=$?
  [[ ${ec} -eq 0 || ${ec} -eq 124 ]]
}

echo "=== Stage 1: Checking mamba, environments, and required input files ==="
command -v mamba
mamba env list | grep -E "^[[:space:]]*${ALIGN_ENV}[[:space:]]"
mamba env list | grep -E "^[[:space:]]*${ASSEMBLY_ENV}[[:space:]]"
mamba env list | grep -E "^[[:space:]]*${MAGS_ENV}[[:space:]]"
test -f "${RAW_R1}"
test -f "${RAW_R2}"
test -f "${REFERENCE_SHIGELLA}"
test -f "${REFERENCE_ALL}"
test -f "${REFERENCE_HUMAN_CHR22}"
test -f "${MASH_DB_DEFAULT}"
test -f "${UNKNOWN_R1}"
test -f "${UNKNOWN_R2}"

echo "=== Stage 2: Checking required tools in ${ALIGN_ENV} ==="
run_in_alignment fastqc --version
run_in_alignment cutadapt --version
run_in_alignment trimmomatic -version
run_in_alignment bowtie2 --version
run_in_alignment samtools --version
run_in_alignment mash --version
run_in_alignment metaphlan --version
run_in_alignment multiqc --version

echo "=== Stage 3: Checking required tools in ${ASSEMBLY_ENV} ==="
run_in_assembly fastqc --version
run_in_assembly cutadapt --version
run_in_assembly trimmomatic -version
run_in_assembly bowtie2 --version
run_in_assembly samtools --version
run_in_assembly spades.py --version
run_in_assembly clumpify.sh --version 2>&1 | head -5 || true
run_in_assembly flash --version
run_in_assembly multiqc --version

echo "=== Stage 4: Checking required tools in ${MAGS_ENV} ==="
run_in_mags run_MaxBin.pl 2>&1 | head -5 || true
run_in_mags checkm -h 2>&1 | head -5
run_in_mags gtdbtk -h 2>&1 | head -5
run_in_mags prokka --version
run_in_mags abricate --version

echo "=== Stage 5: FastQC on raw reads (practical 22) ==="
mkdir -p "${WORKDIR}/qc"
run_in_alignment fastqc "${RAW_R1}" "${RAW_R2}" -o "${WORKDIR}/qc"
test -f "${WORKDIR}/qc/mixed_bacterial_community_5M_R1_fastqc.html"
test -f "${WORKDIR}/qc/mixed_bacterial_community_5M_R2_fastqc.html"

echo "=== Stage 6: Adapter trimming with cutadapt (practical 22) ==="
mkdir -p "${WORKDIR}/cutadapt"
run_in_alignment cutadapt \
  -a "${ADAPTER_FWD}" -A "${ADAPTER_REV}" \
  -o "${WORKDIR}/cutadapt/mixed_R1_noadapt.fastq.gz" \
  -p "${WORKDIR}/cutadapt/mixed_R2_noadapt.fastq.gz" \
  "${RAW_R1}" "${RAW_R2}"
test -s "${WORKDIR}/cutadapt/mixed_R1_noadapt.fastq.gz"
test -s "${WORKDIR}/cutadapt/mixed_R2_noadapt.fastq.gz"

echo "=== Stage 7: Quality filtering with trimmomatic (practical 22) ==="
mkdir -p "${WORKDIR}/trimmed"
run_in_alignment trimmomatic PE -phred33 \
  "${WORKDIR}/cutadapt/mixed_R1_noadapt.fastq.gz" \
  "${WORKDIR}/cutadapt/mixed_R2_noadapt.fastq.gz" \
  "${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_fwd_unpaired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_rev_paired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_rev_unpaired.fq.gz" \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
test -s "${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz"
test -s "${WORKDIR}/trimmed/mixed_fwd_unpaired.fq.gz"
test -s "${WORKDIR}/trimmed/mixed_rev_paired.fq.gz"
test -s "${WORKDIR}/trimmed/mixed_rev_unpaired.fq.gz"

echo "=== Stage 8: FastQC on trimmed reads (practical 22) ==="
run_in_alignment fastqc \
  "${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_fwd_unpaired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_rev_paired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_rev_unpaired.fq.gz" \
  -o "${WORKDIR}/qc"

echo "=== Stage 9: MASH sketching and screening (practical 22) ==="
mkdir -p "${WORKDIR}/mash"
run_in_alignment mash sketch -i -s 500 -k 16 \
  -o "${WORKDIR}/mash/mixed_community" \
  "${REFERENCE_ALL}"
test -f "${WORKDIR}/mash/mixed_community.msh"
run_in_alignment mash info "${WORKDIR}/mash/mixed_community.msh"

run_in_alignment mash sketch -i -s 5000 -k 21 \
  -o "${WORKDIR}/mash/mixed_community_hr" \
  "${REFERENCE_ALL}"
test -f "${WORKDIR}/mash/mixed_community_hr.msh"
run_in_alignment mash info "${WORKDIR}/mash/mixed_community_hr.msh"

# Screen trimmed reads against the custom sketch database
run_in_alignment mash screen -w \
  "${WORKDIR}/mash/mixed_community_hr.msh" \
  "${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_rev_paired.fq.gz"

# Screen against the RefSeq default database (large)
run_in_alignment mash screen -w \
  "${MASH_DB_DEFAULT}" \
  "${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_rev_paired.fq.gz"

echo "=== Stage 10: MetaPhlAn profiling (practical 22) ==="
mkdir -p "${WORKDIR}/metaphlan"
run_in_alignment metaphlan \
  "${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz,${WORKDIR}/trimmed/mixed_rev_paired.fq.gz" \
  --bowtie2out "${WORKDIR}/metaphlan/metaphlan.bowtie2.bz2" \
  --nproc 5 \
  --input_type fastq \
  -o "${WORKDIR}/metaphlan/profiled_metagenome.txt"

echo "=== Stage 11: Bowtie2 index and alignment to Shigella genome (practical 24) ==="
mkdir -p "${WORKDIR}/shigella_ref"
run_in_alignment bowtie2-build -q "${REFERENCE_SHIGELLA}" "${WORKDIR}/shigella_ref/shigella_genome"
run_in_alignment bowtie2 \
  -x "${WORKDIR}/shigella_ref/shigella_genome" \
  -1 "${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz" \
  -2 "${WORKDIR}/trimmed/mixed_rev_paired.fq.gz" \
  -S "${WORKDIR}/shigella_default_alignment.sam"
test -s "${WORKDIR}/shigella_default_alignment.sam"

echo "=== Stage 12: SAM to BAM conversion, sorting, and indexing (practical 24) ==="
run_in_alignment samtools sort -@ 4 -O BAM \
  -o "${WORKDIR}/shigella_default_alignment_sorted.bam" \
  "${WORKDIR}/shigella_default_alignment.sam"
run_in_alignment samtools index "${WORKDIR}/shigella_default_alignment_sorted.bam"
test -s "${WORKDIR}/shigella_default_alignment_sorted.bam"
test -f "${WORKDIR}/shigella_default_alignment_sorted.bam.bai"
run_in_alignment samtools flagstat "${WORKDIR}/shigella_default_alignment_sorted.bam"

echo "=== Stage 13: Alignment to all community genomes (practical 24) ==="
mkdir -p "${WORKDIR}/all_genomes_ref"
run_in_alignment bowtie2-build -q "${REFERENCE_ALL}" "${WORKDIR}/all_genomes_ref/all_genomes"
run_in_alignment bowtie2 -p 8 --fast --no-unal \
  -x "${WORKDIR}/all_genomes_ref/all_genomes" \
  -1 "${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz" \
  -2 "${WORKDIR}/trimmed/mixed_rev_paired.fq.gz" \
  -S "${WORKDIR}/all_genomes_alignment.sam"
run_in_alignment samtools sort -@ 4 -O BAM \
  -o "${WORKDIR}/all_genomes_alignment_sorted.bam" \
  "${WORKDIR}/all_genomes_alignment.sam"
run_in_alignment samtools index "${WORKDIR}/all_genomes_alignment_sorted.bam"
test -s "${WORKDIR}/all_genomes_alignment_sorted.bam"
test -f "${WORKDIR}/all_genomes_alignment_sorted.bam.bai"
run_in_alignment samtools flagstat "${WORKDIR}/all_genomes_alignment_sorted.bam"
run_in_alignment samtools stats "${WORKDIR}/all_genomes_alignment_sorted.bam" \
  > "${WORKDIR}/alignment_statistics.txt"
run_in_alignment samtools coverage "${WORKDIR}/all_genomes_alignment_sorted.bam"

echo "=== Stage 14: QC and pre-processing of unknown pathogen reads (practical 23) ==="
mkdir -p "${WORKDIR}/qc"
run_in_assembly fastqc "${UNKNOWN_R1}" "${UNKNOWN_R2}" -o "${WORKDIR}/qc"
run_in_assembly cutadapt \
  -a "${ADAPTER_FWD}" -A "${ADAPTER_REV}" \
  -o "${WORKDIR}/unknown_noadapt_R1.fastq" \
  -p "${WORKDIR}/unknown_noadapt_R2.fastq" \
  "${UNKNOWN_R1}" "${UNKNOWN_R2}"
test -s "${WORKDIR}/unknown_noadapt_R1.fastq"
test -s "${WORKDIR}/unknown_noadapt_R2.fastq"
run_in_assembly trimmomatic PE -phred33 \
  "${WORKDIR}/unknown_noadapt_R1.fastq" \
  "${WORKDIR}/unknown_noadapt_R2.fastq" \
  "${WORKDIR}/unknown_fwd_paired.fq.gz" \
  "${WORKDIR}/unknown_fwd_unpaired.fq.gz" \
  "${WORKDIR}/unknown_rev_paired.fq.gz" \
  "${WORKDIR}/unknown_rev_unpaired.fq.gz" \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
test -s "${WORKDIR}/unknown_fwd_paired.fq.gz"
test -s "${WORKDIR}/unknown_rev_paired.fq.gz"

echo "=== Stage 15: Align unknown pathogen reads to human chr22, extract unmapped (practical 23) ==="
mkdir -p "${WORKDIR}/human_ref"
run_in_assembly bowtie2-build -q "${REFERENCE_HUMAN_CHR22}" "${WORKDIR}/human_ref/human_chr22"
run_in_assembly bowtie2 --qc-filter -p 8 --local \
  -x "${WORKDIR}/human_ref/human_chr22" \
  -U "${WORKDIR}/unknown_fwd_paired.fq.gz,${WORKDIR}/unknown_rev_paired.fq.gz,${WORKDIR}/unknown_fwd_unpaired.fq.gz,${WORKDIR}/unknown_rev_unpaired.fq.gz" \
  --un "${WORKDIR}/unknown_enriched.fastq" \
  > /dev/null
test -s "${WORKDIR}/unknown_enriched.fastq"

echo "=== Stage 16: SPAdes single genome assembly (practical 23) ==="
run_in_assembly spades.py \
  -s "${WORKDIR}/unknown_enriched.fastq" \
  -o "${WORKDIR}/unknown_genome"

echo "=== Stage 17: Deduplication and read merging for de novo assembly (practical 32) ==="
run_in_assembly clumpify.sh \
  in="${WORKDIR}/trimmed/mixed_fwd_paired.fq.gz" \
  in2="${WORKDIR}/trimmed/mixed_rev_paired.fq.gz" \
  out="${WORKDIR}/mixed_fwd_dedup.fq.gz" \
  out2="${WORKDIR}/mixed_rev_dedup.fq.gz" \
  dedupe=t
test -s "${WORKDIR}/mixed_fwd_dedup.fq.gz"
test -s "${WORKDIR}/mixed_rev_dedup.fq.gz"

run_in_assembly flash \
  -d "${WORKDIR}" -o "flash_output" \
  "${WORKDIR}/mixed_fwd_dedup.fq.gz" \
  "${WORKDIR}/mixed_rev_dedup.fq.gz"
test -s "${WORKDIR}/flash_output.extendedFrags.fastq"
test -s "${WORKDIR}/flash_output.notCombined_1.fastq"
test -s "${WORKDIR}/flash_output.notCombined_2.fastq"

# Append trimmomatic unpaired reads to the merged reads file for assembly input
zcat \
  "${WORKDIR}/trimmed/mixed_fwd_unpaired.fq.gz" \
  "${WORKDIR}/trimmed/mixed_rev_unpaired.fq.gz" \
  >> "${WORKDIR}/flash_output.extendedFrags.fastq"

echo "=== Stage 18: Metaspades metagenomics assembly (practical 32) ==="
run_in_assembly metaspades.py -t 8 \
  -1 "${WORKDIR}/flash_output.notCombined_1.fastq" \
  -2 "${WORKDIR}/flash_output.notCombined_2.fastq" \
  -s "${WORKDIR}/flash_output.extendedFrags.fastq" \
  -o "${WORKDIR}/mixed_comm_de_novo"

echo "=== Stage 19: MAG binning with MaxBin2 (practical 42) ==="
mkdir -p "${WORKDIR}/maxbin"
run_in_mags run_MaxBin.pl \
  -thread 8 \
  -contig "${WORKDIR}/mixed_comm_de_novo/contigs.fasta" \
  -out "${WORKDIR}/maxbin/mixed_comm" \
  -reads "${WORKDIR}/flash_output.notCombined_1.fastq" \
  -reads2 "${WORKDIR}/flash_output.notCombined_2.fastq" \
  -reads3 "${WORKDIR}/flash_output.extendedFrags.fastq"

echo "=== Stage 20: CheckM quality assessment (practical 42) ==="
run_in_mags checkm taxonomy_wf \
  domain Bacteria -x fasta \
  "${WORKDIR}/maxbin/" "${WORKDIR}/checkm/"

echo "=== Stage 21: GTDB-tk taxonomy classification (practical 42, with timeout) ==="
run_with_timeout 120 mamba run -n "${MAGS_ENV}" gtdbtk classify_wf \
  --out_dir "${WORKDIR}/gtdbtk/" \
  --genome_dir "${WORKDIR}/maxbin/" \
  -x fasta --skip_ani_screen --cpus 8

echo "=== Stage 22: Prokka gene annotation (practical 42) ==="
MAG_FILE=$(ls "${WORKDIR}/maxbin/"*.fasta | head -n 1)
run_in_mags prokka \
  --outdir "${WORKDIR}/prokka" \
  "${MAG_FILE}"

echo "=== Stage 23: Abricate AMR gene screening (practical 42) ==="
MAG_FILE=$(ls "${WORKDIR}/maxbin/"*.fasta | head -n 1)
run_in_mags abricate "${MAG_FILE}"

echo "All training environment checks passed."
