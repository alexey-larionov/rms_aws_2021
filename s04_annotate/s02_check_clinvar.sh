#!/bin/bash
# s02_check_clinvar.sh

# Alexey Larionov 29Mar2021

# Use:
# ./s02_check_clinvar.sh &> s02_check_clinvar.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag.MA-split.ID"

scripts_folder="${base_folder}/scripts/s04_annotate"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s04_annotate"
data_vcf="${data_folder}/${base_name}.vcf.gz"

clinvar_folder="${base_folder}/resources/ClinVar/GRCh38/v_20210110"
clinvar_vcf="${clinvar_folder}/clinvar_20210110.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "data_vcf: ${data_vcf}"
echo "clinvar_vcf: ${clinvar_vcf}"
echo ""

# Look at contig names

echo "Contig names in data:"
echo ""
bcftools view -h "${data_vcf}" | grep ^##contig | wc -l
bcftools view -h "${data_vcf}" | grep ^##contig | head -n 30
echo "..."
bcftools view -h "${data_vcf}" | grep ^##contig | tail
echo ""

echo "Contig names in ClinVar:"
echo ""
bcftools view -h "${clinvar_vcf}" | grep ^##contig | wc -l
bcftools view -h "${clinvar_vcf}" | grep ^##contig
echo ""

# Look at reference

echo "Reference in data:"
echo ""
bcftools view -h "${data_vcf}" | grep ^##reference
echo ""

echo "Reference in ClinVar:"
echo ""
bcftools view -h "${clinvar_vcf}" | grep ^##reference
echo ""

# Completion message
echo "Done"
date
