#!/bin/bash
#SBATCH -J Alt.promoter
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/cutadapt"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_Alt.promoter/${dataset}"
mkdir ${output}
mkdir ${output}/salmon
ref="/data/taoyuhuan/projects/exOmics_RNA/bin"
threads=16
libtype=A

echo "Start quantify Alternative promoter for ${dataset} at `date`."
for sample in $(ls ${input} | grep fastq.gz | cut -d "." -f 1)
do
echo "Start quantify ${sample} at `date`"
mkdir ${output}/salmon/${sample}
salmon quant --threads ${threads} --libType A -i ${ref}/genome/salmon-index -r ${input}/${sample}.fastq.gz --validateMappings --gcBias -o ${output}/salmon/${sample}
echo "End quantify ${sample} at `date`"
done


/data/taoyuhuan/projects/exOmics_RNA/bin/scripts/prepareTxQuantMat.py -i ${output}/salmon -o ${output}/TPM-by-tx_${dataset}.txt -t ${ref}/genome/tx-info.txt -q TPM

/data/taoyuhuan/projects/exOmics_RNA/bin/scripts/getPromoterActivity.py -i ${output}/TPM-by-tx_${dataset}.txt -o ${output}/TPM-by-promoter_${dataset}.txt -p ${ref}/genome/promoter/tx2tss.10.txt

echo "Finish quantify Alternative promoter for ${dataset} at `date`."
