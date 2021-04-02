#!/bin/bash
# s01_add_ids.sh

# Alexey Larionov 02Apr2021

# Use:
# ./s01_add_ids.sh &> s01_add_ids.log

# Note
# RS numbers will be added again by VEP

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="zhang_hg38.bwa.QC.MA-flag.MA-split"

scripts_folder="${base_folder}/scripts/s03_annotate"
cd "${scripts_folder}"

source_folder="${base_folder}/data/s02_split_multiallelic_sites"
source_vcf="${source_folder}/${base_name}.vcf.gz"
output_folder="${base_folder}/data/s03_annotate"
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

# Check result
echo "Source VCF"
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL\t%INFO\n' "${source_vcf}" | head
echo ""
echo "Output VCF"
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL\t%INFO\n' "${output_vcf}" | head
echo ""

# Completion message
echo "Done"
date
