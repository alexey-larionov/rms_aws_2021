#!/bin/bash
# s01_add_ids.sh

# Alexey Larionov 29Mar2021

# Use:
# ./s01_add_ids.sh &> s01_add_ids.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag.MA-split"

scripts_folder="${base_folder}/scripts/s04_annotate"
cd "${scripts_folder}"

source_folder="${base_folder}/data/s03_split_multiallelic_sites"
source_vcf="${source_folder}/${base_name}.vcf.gz"
output_folder="${base_folder}/data/s04_annotate"
rm -fr "${output_folder}"
mkdir -p "${output_folder}"
output_vcf="${output_folder}/${base_name}.ID.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Adding variant IDs
echo "Adding variant IDs ..."
bcftools annotate "${source_vcf}" \
--set-id '%CHROM\_%POS\_%REF\_%ALT' \
--output "${output_vcf}" \
--output-type z \
--threads 4
echo ""

#Indexing
echo "Indexing ..."
bcftools index "${output_vcf}"
echo ""

# Completion message
echo "Done"
date
