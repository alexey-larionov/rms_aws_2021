#!/bin/bash
# s01_compress_vcf.sh

# Alexey Larionov 01Apr2021

# Use:
# ./s01_compress_vcf.sh &> s01_compress_vcf.log

# Note
# The compression is needed not only to reduce the file size, but also because
# some bcftoolls commands require the file to be compressed and indexed

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="zhang_hg38.bwa.QC"

scripts_folder="${base_folder}/scripts/s01_copy_vcf"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s01_copy_vcf"
source_vcf="${data_folder}/${base_name}.vcf"
output_vcf="${data_folder}/${base_name}.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Compressing VCF file
echo "Compressing ..."
bcftools view "${source_vcf}" \
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
