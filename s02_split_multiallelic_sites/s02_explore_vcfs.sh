#!/bin/bash

# s02_explore_vcfs.sh
# Alexey Larionov, 28Mar2021

# Intended use:
# ./s02_explore_vcfs.sh &> s02_explore_vcfs.log

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

data_folder="${base_folder}/data/s02_split_multiallelic_sites"

multiallelic_vcf="${data_folder}/${base_name}.MA.vcf.gz"
flagged_vcf="${data_folder}/${base_name}.MA-flag.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "vcf: ${vcf}"
echo "data_folder: ${data_folder}"
echo ""

# Count variants

echo "Multiallelic variant counts:"
echo ""
bcftools +counts "${multiallelic_vcf}"
echo ""

echo "Variant counts in vcf with flagged MA sites:"
echo ""
bcftools +counts "${flagged_vcf}"
echo ""

# Explore vcf-s

echo "String added to the header of flagged file"
echo ""
bcftools view -h "${flagged_vcf}" | grep "##INFO=<ID=MULTIALLELIC"
echo ""

#bcftools view -H "${flagged_vcf}" | grep "MULTIALLELIC" | head -n 1 # chr1:935954
#bcftools view -H "${flagged_vcf}" | head -n 1 # chr1:13273

echo "A multiallelic site in source vcf"
echo ""
bcftools view -H -r chr1:935954 "${multiallelic_vcf}" | cut -f -8
echo ""

echo "A multiallelic site in flagged vcf"
echo ""
bcftools view -H -r chr1:935954 "${flagged_vcf}" | cut -f -8
echo ""

echo "A biallelic site in flagged vcf"
echo ""
bcftools view -H -r chr1:13273 "${flagged_vcf}" | cut -f -8
echo ""

# Completion message
echo "Done"
date
echo ""
