###README for small_RNAseq_pipeline
-------------------------------------------------------------------------
##Before start

0.Background:
QIAseq_pipeline: for QIAseq method, a commerical small RNA library kit with UMIs(unique molecular identifiers)[http://www.ebiotrade.com/custom/qiagen/170329/index.html]
small_pipeline: a more common pipeline for small RNA seq without UMIs

1.Software requirements:
Python 3.7.4
R 3.6.1
snakemake 5.5.4
trim_galore 0.6.7
FastQC 0.11.8
STAR 2.5.3a_modified [note: v2.5.3a is also accepted]
samtools 1.10 [note: using htslib 1.10]
UMI-tools 1.1.0
bowtie2 2.3.5.1

--------------------------------------------------------------------------
##Usage
0. Setup working environment with the following directories:
	./data
	./sample_id
	./pipeline
	./scripts
	./ref
	./output

1. Prepare rawdata file in ./data
	all fastq files(gz compressed) should store in a subfolder, such as ./data/GSE104251,
	all fastq files should end with "fastq.gz".
2. Provide sample ids in ./sample_id
	sample ids should be provided in a txt file named consistant with you raw file subfolder, such as GSE104251.txt,
	sample ids are strings before "fastq.gz",
	sample ids should be "\n" seprated,
	make sure all sample ids exist in you raw file subfolder.
3. Copy pipelines to your working environment directory
	copy from: /data/taoyuhuan/projects/small_RNAseq_pipeline/pipeline
4. Copy scripts to your working environment directory
        copy from: /data/taoyuhuan/projects/small_RNAseq_pipeline/scripts
	scripts are some tools call by pipeline
5. Copy ref to your working environment directory
        copy from: /data/taoyuhuan/projects/small_RNAseq_pipeline/ref
6. Create an "output" directory

7. Initialize your pipeline by vim
	write dataset name to the 1st line, such as: dataset = "GSE104251",
	set your mapping sequence for different RNA biptypes, such as: sequences = ["spikein_small","univec","rRNA","mtRNA","miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","genome","circRNA"].

8. Copy config and scripts to your working environment directory
	copy from: /data/taoyuhuan/projects/small_RNAseq_pipeline/Snakemake_cluster.json, /data/taoyuhuan/projects/small_RNAseq_pipeline/Snakemake_cluster.json

9. Start pipeline by the following command:
	sbatch Snakemake_cluster.sh {dataset}
	{dataset} can be any with rawfile and sample id provided, such as GSE104251.

------------------------------------------------------------------------------
