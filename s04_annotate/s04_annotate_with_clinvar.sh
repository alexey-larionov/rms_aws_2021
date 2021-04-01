#!/bin/bash
# s04_annotate_with_clinvar.sh

# Alexey Larionov 29Mar2021

# Use:
# ./s04_annotate_with_clinvar.sh &> s04_annotate_with_clinvar.log

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
source_vcf="${data_folder}/${base_name}.vcf.gz"
output_vcf="${data_folder}/${base_name}.ClinVar.vcf.gz"
clinvar_vcf="${data_folder}/clinvar_20210110.chr.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo "clinvar_vcf: ${clinvar_vcf}"
echo ""

# Annotate using bcftools
echo "Annotating ..."
bcftools annotate "${source_vcf}" \
--annotations "${clinvar_vcf}" \
--columns INFO \
--output "${output_vcf}" \
--output-type z \
--threads 4
echo ""

# Index output vcf
echo "Indexing ..."
bcftools index "${output_vcf}"
echo ""

# Summary of INFO fields
echo "Number of INFO fields in source vcf:"
bcftools view -h "${source_vcf}" | grep ^##INFO | wc -l
echo ""

echo "Number of INFO fields in ClinVar-annotated vcf:"
bcftools view -h "${output_vcf}" | grep ^##INFO | wc -l
echo ""

echo "List of INFO fields in ClinVar-annotated vcf:"
bcftools view -h "${output_vcf}" | grep ^##INFO
echo ""

# Progress report
echo "Selected lines from output_vcf:"
echo ""
bcftools query \
-i 'ALLELEID != "."' \
-f '%ID\t%ALLELEID\t%CLNSIG\t%CLNDN\n' "${output_vcf}" | head
echo ""

bcftools query \
-i 'ALLELEID != "." & CLNSIG != "Benign"' \
-f '%ID\t%ALLELEID\t%CLNSIG\t%CLNDN\n' "${output_vcf}" | head
echo ""

# Completion message
echo "Done"
date
echo ""
