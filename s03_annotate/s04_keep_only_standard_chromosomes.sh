#!/bin/bash
# s04_keep_only_standard_chromosomes.sh

# Alexey Larionov 02Apr2021

# Use:
# ./s04_keep_only_standard_chromosomes.sh &> s04_keep_only_standard_chromosomes.log

# Stop at runtime errors
set -e

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"
base_name="zhang_hg38.bwa.QC.MA-flag.MA-split.ID"

scripts_folder="${base_folder}/scripts/s03_annotate"
cd "${scripts_folder}"

data_folder="${base_folder}/data/s03_annotate"
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
--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM \
--output "${output_vcf}" \
--output-type z \
--threads 4
echo ""

# Index output vcf
echo "Indexing ..."
bcftools index "${output_vcf}"
echo ""

# Explore result

echo "Counts in source vcf:"
bcftools +counts "${source_vcf}"
echo ""

echo "Counts in output vcf:"
bcftools +counts "${output_vcf}"
echo ""

echo "Number of contigs in source vcf header:"
bcftools view -h "${source_vcf}" | grep ^##contig | wc -l
echo ""

echo "Number of contigs in output vcf header:"
bcftools view -h "${output_vcf}" | grep ^##contig | wc -l
echo ""

# Completion message
echo "Done"
date
