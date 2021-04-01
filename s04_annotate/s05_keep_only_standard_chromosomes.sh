#!/bin/bash
# s05_keep_only_standard_chromosomes.sh

# Alexey Larionov 30Mar2021

# Use:
# ./s05_keep_only_standard_chromosomes.sh &> s05_keep_only_standard_chromosomes.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="IHCAPX8_dragen_joint.hard-filtered.PF.MA-flag.MA-split.ID.ClinVar"

scripts_folder="${base_folder}/scripts/s04_annotate"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s04_annotate"
source_vcf="${data_folder}/${base_name}.vcf.gz"
output_vcf="${data_folder}/${base_name}.std-Chr.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo ""

# Annotate using bcftools
echo "Select variants ..."
bcftools view "${source_vcf}" \
--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \
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
