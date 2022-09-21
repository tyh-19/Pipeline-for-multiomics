#!/bin/bash
#SBATCH -J rMATs
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

dataset=$1
positive=$2
negative=$3
input="/data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/$1"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_splicing/$1/${positive}vs${negative}"
mkdir -p ${output}
gtf="/data/taoyuhuan/projects/exOmics_RNA/bin/genome/gtf/gencode.v27.annotation.gtf"

echo "Start call alternative splicing for $1: ${positive} vs ${negative} at `date`"
rmats.py --b1 ${input}/${positive}.txt --b2 ${input}/${negative}.txt --gtf ${gtf} --od ${output} -t single  --libType fr-unstranded --readLength 100

