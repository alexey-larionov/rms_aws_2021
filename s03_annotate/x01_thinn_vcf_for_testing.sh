#!/bin/bash
# x01_thinn_vcf_for_testing.sh

# Alexey Larionov 02Apr2021

# Intended use:
# ./x01_thinn_vcf_for_testing.sh &> x01_thinn_vcf_for_testing.log

# Note
# This script is kept here just in case: if a VCF file needs to be
# thinned for VEP debugging: its hard to debug, when each run takes ~30 min

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="data.name" # update the name

scripts_folder="${base_folder}/scripts/s03_annotate"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s03_annotate" # update the folder if needed
source_vcf="${data_folder}/${base_name}.vcf.gz"

output_vcf="${data_folder}/thinned.recode.vcf"
# the prefix ("${data_folder}/thinned") will be specifyed in the --out option (see below)
# the default suffix ".recode.vcf" will be added by vcftools

# Progress report
vcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Thinn vcf using vcftools
# - this will thinn the file to keep minimal distance of 100000 bases between variants
# - the --out option specifyes "prefix" of the output file
# - the default suffix ".recode.vcf" will be added to the output file name by vcftools
echo "Thinning ..."
vcftools \
--gzvcf "${source_vcf}" \
--thin 100000 \
--recode \
--recode-INFO-all \
--out "${data_folder}/thinned"
echo ""

# Explore result

echo "Number of contigs in source vcf header:"
bcftools view -h "${source_vcf}" | grep ^##contig | wc -l
echo ""

echo "Number of contigs in output vcf header:"
bcftools view -h "${output_vcf}" | grep ^##contig | wc -l
echo ""

echo "Counts in source vcf:"
bcftools +counts "${source_vcf}"
echo ""

echo "Counts in output vcf:"
bcftools +counts "${output_vcf}"
echo ""

# Completion message
echo "Done"
date
echo ""
