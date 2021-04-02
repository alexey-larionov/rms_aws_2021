#!/bin/bash
# s03_rename_clinvar_chr.sh

# Alexey Larionov 02Apr2021

# Use:
# ./s03_rename_clinvar_chr.sh &> s03_rename_clinvar_chr.log

# Start message
echo $0
date
echo ""

# Files and folders
base_folder="/home/share"

scripts_folder="${base_folder}/scripts/s03_annotate"
cd "${scripts_folder}"

clinvar_folder="${base_folder}/resources/ClinVar/GRCh38/v_20210110"
initial_clinvar_vcf="${clinvar_folder}/clinvar_20210110.vcf.gz"

data_folder="${base_folder}/data/s03_annotate"
updated_clinvar_vcf="${data_folder}/clinvar_20210110.chr.vcf.gz"

# This file is manually prepared using results of the previous script
chromosomes_translation_file="${scripts_folder}/s03_clinvar_chr_translation.txt"

# This file is used just for comparison
data_vcf="${data_folder}/zhang_hg38.bwa.QC.MA-flag.MA-split.ID.vcf.gz"

# Progress report
bcftools --version
echo ""
echo "initial_clinvar_vcf: ${initial_clinvar_vcf}"
echo "updated_clinvar_vcf: ${updated_clinvar_vcf}"
echo "chromosomes_translation_file: ${chromosomes_translation_file}"
echo ""

# Check headers from clinvar and gel vcf-s

echo "--- Contigs in public ClinVar VCF header ---"
echo "" # h = only header
bcftools view -h "${initial_clinvar_vcf}" | grep "^##contig"
echo ""

echo "--- Number of records in public ClinVar VCF ---"
echo "" # H = no header
bcftools view -H "${initial_clinvar_vcf}" | wc -l
echo ""

echo "--- Contigs in GEL VCF ---"
echo ""
bcftools view -h "${data_vcf}" | grep "^##contig" | head -n 27
echo "..."
echo ""

echo "--- Translation file ---"
echo ""
cat "${chromosomes_translation_file}"
echo ""

# Update ClinVar VCF
echo "Updating ClinVar VCF ..."
bcftools annotate "${initial_clinvar_vcf}" \
--rename-chrs "${chromosomes_translation_file}" \
--output "${updated_clinvar_vcf}" \
--output-type z \
--threads 4
echo ""

echo "Indexing ..."
bcftools index "${updated_clinvar_vcf}"
echo ""

# Check result
echo "--- Contigs in the updated ClinVar VCF header ---"
echo ""
bcftools view -h "${updated_clinvar_vcf}" | grep "^##contig"
echo ""
echo "--- Number of records in updated ClinVar VCF ---"
echo ""
bcftools view -H "${updated_clinvar_vcf}" | wc -l
echo ""

# Progress report
echo "Done"
date
