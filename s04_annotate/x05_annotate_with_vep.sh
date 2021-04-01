#!/bin/bash
# x05_annotate_with_vep.sh

# Alexey Larionov 30Mar2021

# Intended use:
# ./x05_annotate_with_vep.sh &> x05_annotate_with_vep.log

# Notes
# - Requires at least 16GB of RAM
#   Runs out of memory if used with m4.large (2 cores and 8 GB RAM)
# - Took ~12hrs per WGS (~7M variants) when using 4 cores
#   (m4.xlarge: 4 cores, 16GB RAM, Hight network speed)
# - Took ~10 min per VCF with ~30k variants
# - Dont edit script during its execution: somehow it may confuse bash

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
#base_name="IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag.MA-split.ID.ClinVar"

scripts_folder="${base_folder}/scripts/s04_annotate"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s04_annotate"

source_vcf="${data_folder}/thinned.recode.vcf"
output_vcf="${data_folder}/thinned.recode.vep.vcf.gz"
vep_report="${data_folder}/thinned.recode.vep.html"

#source_vcf="${data_folder}/${base_name}.vcf.gz"
#output_vcf="${data_folder}/${base_name}.VEP.vcf.gz"
#vep_report="${data_folder}/${base_name}.VEP.html"

vep_cache_info="${data_folder}/vep_cache_info.txt"

# VEP
vep="${base_folder}/tools/ensembl/vep-102-b38/vep"

# VEP cache
# cache_folder isn't needed because cache was installed to the default location in ~/.vep folder
# cache_folder="/home/ec2-user/.vep
cache_version="102"
cache_assembly="GRCh38"

# VEP plugins folder
# plugins_folder isn't needed because plugins were installed to the default location in ~/.vep folder
# plugins_folder="/home/ec2-user/.vep/Plugins"

# Data for plugins
#cadd_data_folder="${base_folder}/resources/cadd/v1.6_b38"
#cadd_snv_data="${cadd_data_folder}/whole_genome_SNVs.tsv.gz"
#cadd_indels_data="${cadd_data_folder}/gnomad.genomes.r3.0.indel.tsv.gz"

# Reference genome
# A reference genome is provided in cache, so it is, most likely, is used by VEP by default
# The exact reference used by Illumina BaseSpace Dragen is not yet clear: it yet needs to be asked from Illumina
#b38_fasta="/home/ec2-user/.vep/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"

# Num of threads: update to fit the hardware
#n_threads="16"
n_threads="4"

# PATH and PERL5LIB were configured during the AMI setup, so no updates are needed here
#export PATH="${tools_folder}/htslib/htslib-1.10.2/bin:$PATH"
#export PERL5LIB="${tools_folder}/vep/v101_b38/cpanm/lib/perl5:$PERL5LIB"

# Progress report
echo "--- Files ---"
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo "vep_report: ${vep_report}"
echo ""
echo "--- VEP cache ---"
echo ""
#echo "cache_folder: ${cache_folder}" # not needed because cache is in the default location in ~/.vep
echo "cache_version: ${cache_version}"
echo "cache_assembly: ${cache_assembly}"
echo ""
"${vep}" --show_cache_info
echo ""
cp "/home/ec2-user/.vep/homo_sapiens/102_GRCh38/info.txt" "${vep_cache_info}"
echo "See more information about used cache in the following file:"
echo "${vep_cache_info}"
echo ""
#echo "--- VEP plugins ---"
#echo ""
#echo "plugins_folder: ${plugins_folder}"
#echo ""
#echo "CADD annotation files:"
#echo "${cadd_snv_data}"
#echo "${cadd_indels_data}"
#echo ""
#echo "--- Reference ---"
#echo ""
#echo "b38_fasta: ${b38_fasta}"
#echo ""
echo "--- Other settings ---"
echo ""
echo "n_threads: ${n_threads}"
echo ""

# Annotate VCF
echo "Annotating ..."
"${vep}" \
--input_file "${source_vcf}" \
--output_file "${output_vcf}" \
--vcf \
--compress_output bgzip \
--force_overwrite \
--fork "${n_threads}" \
--stats_file "${vep_report}" \
--offline \
--cache \
--species homo_sapiens \
--cache_version "${cache_version}" \
--assembly "${cache_assembly}" \
--pick \
--gencode_basic \
--check_ref \
--everything \
--nearest symbol \
--total_length \
--check_existing \
--exclude_null_alleles

echo ""

#  --fasta "${b38_fasta}" \
#  --dir_cache "${cache_folder}" \

#  --dir_plugins "${plugins_folder}" \
#  --plugin CADD,"${cadd_snv_data}","${cadd_indels_data}" \

# Index annotated vcf
echo "Indexing ..."
bcftools index "${output_vcf}"
echo ""

# Added VEP annotations
echo "Added VEP annotations:"
echo ""
bcftools +split-vep -l "${output_vcf}"
echo ""

# Completion message
echo "Done"
date
