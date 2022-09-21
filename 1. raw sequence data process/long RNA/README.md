# cfRNAseq


### Preparation environment
1. Install conda
```
wget -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh
bash Miniconda2-4.5.11-Linux-x86_64.sh
```
2. Create an environemt for cfRNAseq, cfRNAseq.yaml can be find in ./cfRNAseq
```
conda env create -f cfRNAseq.yaml
```
if not success, try create environment by yourself
```
#python==3.9.6
#R=3.6
#Java version "11.0.9.1-internal"
conda create --name cfRNAseq
conda activate cfRNAseq
conda install snakemake
conda install openjdk
conda install -c bioconda snakemake
conda install cutadapt=3.4
conda install star=2.5.3a
conda install samtools
conda install bedtools
conda install subread
```
3. Prepera tools and reference
```
cp /BioII/lulab_b/taoyuhuan/cfRNAseq/tools /{your directory}
cp /BioII/lulab_b/taoyuhuan/cfRNAseq/reference /{your directory}
```
4. Set parameters in cfRNAseq.snakefile
```
#config
tool_dir="cfRNAseq/tools"

#reference
genome_dir="cfRNAseq/reference/genome_index/star"
fasta_dir="cfRNAseq/reference/fasta"
bed_dir="cfRNAseq/reference/bed"
gtf_dir="cfRNAseq/reference/gtf"

#input
data_dir="cfRNAseq/test_RNAseq/fastq"
sample_ids=["CRC-2418488_test_neg","CRC-2276341_test_neg"]

#output
output_dir="cfRNAseq/output"

#cutadapt
thread_cutadapt=4
clean_levels=['cutadapt','trimGC']
adaptor1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adaptor2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

#mapping and operate on bam file
thread_mapping=16
threads_decompress=4
map_steps=['spikein_long','univec','rRNA','hg38_long','circRNA']
map_steps_sortbyName=['rRNA','hg38_long','circRNA']
map_steps_rmdup=['rRNA','hg38_long','circRNA']

#count parameter
strandness="reverse" #{"no","forward","reverse"}
paired_end="True" #{"True","False"}
count_multimap_reads="True" #{"True","False"}
min_mapping_quality=0
count_overlapping_features="True" #{"True","False"}
count_levels=['hg38_long','hg38_long_rmdup','intron_spanning']

#temp
temp_dir="cfRNAseq/temp"
```
5. Run example data
```
snakemake --snakefile RNAseq.snakefile --cores 16
```
