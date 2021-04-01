#!/bin/bash
# s03_split_multiallelic_sites.sh

# Alexey Larionov 01Apr2021

# Use:
# ./s03_split_multiallelic_sites.sh &> s03_split_multiallelic_sites.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="zhang_hg38.bwa.QC.MA-flag"

scripts_folder="${base_folder}/scripts/s02_split_multiallelic_sites"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s02_split_multiallelic_sites"
source_vcf="${data_folder}/${base_name}.vcf.gz"
output_vcf="${data_folder}/${base_name}.MA-split.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Split multiallelic sites
echo "Splitting multiallelic sites ..."
bcftools norm "${source_vcf}" \
--multiallelics -any \
--do-not-normalize \
--output "${output_vcf}" \
--output-type z \
--threads 2
echo ""

# Indexing
echo "Indexing ..."
echo ""
bcftools index "${output_vcf}"

# Completion mesage
echo "Done"
date
