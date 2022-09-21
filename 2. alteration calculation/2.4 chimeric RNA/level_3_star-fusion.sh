#!/bin/bash
#SBATCH -J star-fusion
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#conda activate ChimericRNA

dataset=$1
strandness="reverse" # forward, reverse or no
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/unmapped"
output="/data/taoyuhuan/projects/exOmics_RNA/level_3_chimeric/${dataset}"
mkdir -p ${output}
mkdir -p ${output}/star-fusion

echo "Start counting ${dataset} chimeric RNA at `date`."
for sample in $(ls ${input})
#for sample in $(cat /data/taoyuhuan/projects/exOmics_RNA/GSE133684/raw/sample_ids_reverse.txt)
do
echo "Mapping ${sample} chimeric RNA at `date`."
mkdir -p ${output}/${sample}
STAR --genomeDir /data/taoyuhuan/reference/ChimericRNA/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx \
          --readFilesIn ${input}/${sample}/circRNA_1.fastq.gz ${input}/${sample}/circRNA_2.fastq.gz \
          --twopassMode Basic \
          --runThreadN 16 \
          --outFileNamePrefix ${output}/${sample}/${sample} \
          --outSAMtype BAM Unsorted \
          --outReadsUnmapped Fastx \
          --readFilesCommand gzip -d -c \
          --outSAMmultNmax 1 \
          --seedPerWindowNmax 20 \
	  --outSAMstrandField intronMotif \
          --outSAMunmapped Within \
          --chimSegmentMin 12 \
          --chimJunctionOverhangMin 8 \
          --chimOutJunctionFormat 1 \
          --alignSJDBoverhangMin 10 \
          --alignMatesGapMax 100000 \
          --alignIntronMax 100000 \
          --alignSJstitchMismatchNmax 5 -1 5 5 \
          --outSAMattrRGline ID:GRPundef \
          --chimMultimapScoreRange 3 \
          --chimScoreJunctionNonGTAG -4 \
          --chimMultimapNmax 20 \
          --chimNonchimScoreDropMin 10 \
          --peOverlapNbasesMin 12 \
          --peOverlapMMp 0.1 \
          --alignInsertionFlush Right \
          --alignSplicedMateMapLminOverLmate 0 \
          --alignSplicedMateMapLmin 30

echo "Repairing ${sample} at `date`"
# To fix the bug of STAR_2.5.3a_modified 
/data/taoyuhuan/tools/bbmap/bbmap/repair.sh in=${output}/${sample}/${sample}Unmapped.out.mate1 in2=${output}/${sample}/${sample}Unmapped.out.mate2 out=${output}/${sample}/unmapped_1.fastq.gz out2=${output}/${sample}/unmapped_2.fastq.gz overwrite=t

echo "Fusion predicting ${sample} at `date`"
mkdir -p ${output}/star-fusion/${sample}
STAR-Fusion --genome_lib_dir /data/taoyuhuan/reference/ChimericRNA/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir \
            --CPU 16 \
            -J ${output}/${sample}/${sample}Chimeric.out.junction \
            --output_dir ${output}/star-fusion/${sample}

cat ${output}/star-fusion/${sample}/star-fusion.fusion_predictions.tsv | cut -f 1 | sed "1d" > ${output}/star-fusion/${sample}/fusion_genes.txt

FusionInspector --fusions ${output}/star-fusion/${sample}/fusion_genes.txt \
                --genome_lib /data/taoyuhuan/reference/ChimericRNA/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir \
                --left_fq ${input}/${sample}/circRNA_1.fastq.gz --right_fq ${input}/${sample}/circRNA_2.fastq.gz \
                --output_dir ${output}/star-fusion/${sample} \
                --out_prefix finspector \
                --CPU 16 --annotate \
                --vis
done

Rscript /data/taoyuhuan/projects/exOmics_RNA/level_3_chimeric/STARFusion_matrix.R -i ${output}/star-fusion -o ${output}/star-fusion_matrix

