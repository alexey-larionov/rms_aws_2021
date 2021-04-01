#!/bin/bash
# s06_update_chromosomes_in_header.sh

# Alexey Larionov 30Mar2021

# Use:
# ./s06_update_chromosomes_in_header.sh &> s06_update_chromosomes_in_header.log

# Note:
# The new header was made manually:
# - First the previous header was obtained with
#   bcftools -h view IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag.MA-split.ID.ClinVar.std-Chr.vcf.gz
# - Then all contigs except for chr1-22,X,Y were manually removed

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
output_vcf="${data_folder}/${base_name}.std-Chr.Reheaded.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Extract old header


# Modify old header


# Reheader VCF using bcftools
echo "Reheader ..."
bcftools ... "${source_vcf}" \

--output "${output_vcf}" \
--output-type z \
--threads 4
echo ""

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
