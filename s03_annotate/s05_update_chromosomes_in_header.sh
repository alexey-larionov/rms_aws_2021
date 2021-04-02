#!/bin/bash
# s05_update_chromosomes_in_header.sh

# Alexey Larionov 02Apr2021

# Use:
# ./s05_update_chromosomes_in_header.sh &> s05_update_chromosomes_in_header.log

# Notes:
# This step is optional: just to make some plots look nicier in VEP-html report.
# Of course, its only possible after we filtered al variants outside of
# standard chromosomes during the previous step.

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="zhang_hg38.bwa.QC.MA-flag.MA-split.ID.std-Chr"

scripts_folder="${base_folder}/scripts/s03_annotate"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s03_annotate"
source_vcf="${data_folder}/${base_name}.vcf.gz"
output_vcf="${data_folder}/${base_name}.Reheaded.vcf.gz"

header_old="${data_folder}/header.old"
header_new="${data_folder}/header.new"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Extract old header
bcftools view -h "${source_vcf}" > "${header_old}"

# Modify old header: exclude the alt etc
cat "${header_old}" | \
grep -v ^##contig=\<ID=HLA | \
grep -v ^##contig=\<ID=chrUn | \
grep -v ^##contig=\<ID=chrEBV | \
grep -v ^##contig=\<ID=chr.*_alt | \
grep -v ^##contig=\<ID=chr.*_random \
> "${header_new}"

# Reheader VCF using bcftools
echo "Reheading ..."
bcftools reheader "${source_vcf}" \
--header "${header_new}" \
--output "${output_vcf}" \
--threads 4
echo ""

# Note that bcftools reheader does not accept '--output-type' option

# Index output vcf
echo "Indexing ..."
bcftools index "${output_vcf}"
echo ""

# Explore result

echo "Number of contigs in source vcf header:"
bcftools view -h "${source_vcf}" | grep ^##contig | wc -l
echo ""

echo "Number of contigs in output vcf header:"
bcftools view -h "${output_vcf}" | grep ^##contig | wc -l
echo ""

echo "Contigs in source vcf:"
bcftools view -h "${source_vcf}" | grep ^##contig | head -n 30
echo "..."
echo ""

echo "Contigs in output vcf:"
bcftools view -h "${output_vcf}" | grep ^##contig
echo ""

# Completion message
echo "Done"
date
