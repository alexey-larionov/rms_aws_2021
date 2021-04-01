#!/bin/bash
# s01_flag_multiallelic_sites.sh

# Alexey Larionov 01Apr2021

# Use:
# ./s01_flag_multiallelic_sites.sh &> s01_flag_multiallelic_sites.log

# References & examples

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="zhang_hg38.bwa.QC"

scripts_folder="${base_folder}/scripts/s02_split_multiallelic_sites"
cd "${scripts_folder}"

source_folder="${base_folder}/data/s01_copy_vcf"
source_vcf="${source_folder}/${base_name}.vcf.gz"

output_folder="${base_folder}/data/s02_split_multiallelic_sites"
rm -fr "${output_folder}"
mkdir -p "${output_folder}"
multiallelic_vcf="${output_folder}/${base_name}.MA.vcf.gz"
output_vcf="${output_folder}/${base_name}.MA-flag.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "multiallelic_vcf: ${multiallelic_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Select multiallelic sites
echo "Selecting multiallelic sites ..."
bcftools view "${source_vcf}" \
--min-alleles 3 \
--output "${multiallelic_vcf}" \
--output-type z \
--threads 2

echo ""

# Indexing
echo "Indexing ..."
echo ""
bcftools index "${multiallelic_vcf}"

# Flag multiallelic sites
echo "Adding multiallelic flag ..."
bcftools annotate "${source_vcf}" \
--annotations "${multiallelic_vcf}" \
--mark-sites MULTIALLELIC \
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
