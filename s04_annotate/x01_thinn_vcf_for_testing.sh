#!/bin/bash
# x01_thinn_vcf_for_testing.sh

# Alexey Larionov 30Mar2021

# Use:
# ./x01_thinn_vcf_for_testing.sh &> x01_thinn_vcf_for_testing.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag.MA-split.ID.ClinVar.std-Chr"

scripts_folder="${base_folder}/scripts/s04_annotate"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s04_annotate"
source_vcf="${data_folder}/${base_name}.vcf.gz"

output_vcf="${data_folder}/thinned.recode.vcf"
# the default suffix ".recode.vcf" will be added by vcftools
# to the prefix, specifyed in the --out option (see below)

# Progress report
vcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Thinn vcf using vcftools
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
