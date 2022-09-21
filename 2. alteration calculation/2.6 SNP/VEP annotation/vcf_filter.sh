# !/bin/bash

dataset=$1
input_dir="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP"
output_dir="/data/taoyuhuan/projects/exOmics_RNA/level_3_Editing_SNP_ASE/${dataset}/SNP_filtered"
mkdir -p ${output_dir}

for vcf in $(ls ${input_dir} | grep "rmEDIT.SNP.gz")
do
echo "Start filter ${vcf} at `date`."
zcat ${input_dir}/${vcf} > ${input_dir}/${vcf%.*}.vcf
vcftools --vcf ${input_dir}/${vcf%.*}.vcf --remove-filtered-all --min-meanDP 10 --recode --stdout > ${output_dir}/${vcf%.*}.filtered.vcf
raw_vcf=$(less ${input_dir}/${vcf} | wc -l)
filtered_vcf=$(less ${output_dir}/${vcf%.*}.filtered.vcf | wc -l)
echo "Raw: "${raw_vcf}
echo "Filtered: "${filtered_vcf}
done

