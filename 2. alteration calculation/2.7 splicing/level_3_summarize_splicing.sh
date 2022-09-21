#!/bin/bash
#SBATCH -J summary_rMATs
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output=%j.out
#SBATCH --error=%j.err


dataset=$1
positive=$2
negative=$3
rMATS_input="/data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/$1"
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/$1/${positive}vs${negative}"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/$1/${positive}vs${negative}/summary"
mkdir -p ${output}
mkdir -p ${input}/matrix
mkdir -p ${input}/differential_splicing

cat ${rMATS_input}/${positive}.txt | tr "," "\n" > ${output}/${positive}.txt
cat ${rMATS_input}/${negative}.txt | tr "," "\n" > ${output}/${negative}.txt

echo "Start summarizing ${dataset} at `date`"
for splicing_type in MXE A3SS A5SS SE RI
do
echo "Processing ${splicing_type}:"
/data/taoyuhuan/projects/exOmics_RNA/bin/scripts/summarize-splicing.py --input ${input}/${splicing_type}.MATS.JC.txt --outdir ${output} --type ${splicing_type} --method JC --pos ${output}/${positive}.txt --neg ${output}/${negative}.txt
done
echo "Finish summarizing ${dataset} at `date`"

Rscript /data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/splicing_aggregate.R -i ${input}/summary -m ${input}/matrix/${dataset}_${positive}_vs_${negative}_Inclevel_matrix.txt -s ${input}/differential_splicing/${dataset}_${positive}_vs_${negative}.txt
