# Analysis-pipeline-for-multiomics
### 1.raw sequence data process

**1.1 Long RNA**

**Introduction:** This pipeline is designed for SMARTer-seq pico input v2, which utilizing template switch strategy to add adaptor sequence at 5' of RNA fragments. This strategy introduces ~3 GC pairs at the 5' end of read. We trimmed these artificial GC pairs to allow more precise read assignment. We also applied sequential mapping strategy, sequentially mapped to ERCC spike-in, univec, rRNA, genome and circular RNA. Sequential mapping strategy can reduce chances of multi-mapping, result in higher mapping ratio.

**Requirements：**

Environments and tools:

```
python=3.9.6
R=3.6
snakemake=5.5.4
cutadapt=3.4
star=2.5.3a
samtools=1.10
bedtools=2.30.0
GATK==4.1.9.0
BBMap==38.90
```

Reference and index:

The version of genome sequence is GRCh38. Annotation is GENCODE release 27. More detailed information please refer to [GENCODE website](https://www.gencodegenes.org/human/stats_27.html).

Due to the large size of reference and index, we suggest beginers to download from [GENCODE download page](https://www.gencodegenes.org/human/release_27.html).

**Demostration：**

Before start, we need to set: 

1. all the tools;
2. reference used in this pipeline;
3. path to input files and their identifiers;
4. output directory;
5. basic parameters: strandness of seqencing methods; single_end or paired_end, etc.;
6. threads for cutadaptor, mapping and compress.

After all parameters set up, raw sequence processing can be simlified to 1 command using the snakemake pipeline.

```
snakemake --snakefile RNAseq.snakefile
```

*More details on configuration settings or advanced usage in high-performance cluster, please refer to ./1. raw sequence data process/long RNA/README.md*



**1.2 small RNA**

**Introduction：**

We provide 2 customized pipeline for small RNA processing: small_pipeline and QIAseq_pipeline. small_pipeline is designed for common small RNA sequencing without UMIs (unique molecular identifiers). QIAseq_pipeline is designed for a commerical small RNA library kit with UMIs. 

**Requirements：**

```
Python=3.7.4 
R=3.6
snakemake=5.5.4
trim_galore=0.6.7
FastQC=0.11.8
bowtie2=2.3.5.1
STAR=2.5.3a_modified [note: v2.5.3a is also accepted]
samtools=1.10 [note: using htslib 1.10] 
UMI-tools=1.1.0 
```

**Demostration：**

Before start, we need to prepare: 

1. all raw data(gz compressed fastq files, end with fastq.gz) in a folder;
2. a list of all samples required to be processed, sample identifiers should be seperated by "\n";
3. Path to indexed reference;
4. output directory;
5. mappin sequence: such as ["spikein_small","univec","rRNA","miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","genome","circRNA"];

After all preparation done, raw sequence processing can be simlified to 1 command using the snakemake pipeline.

```
#for library with UMI, such as QIAseq data
snakemake --snakefile QIAseq_pipeline

#for library without UMI
snakemake --snakefile small_pipeline
```

*More details on configuration settings or advanced usage in high-performance cluster, please refer to ./1. raw sequence data process/small RNA/README.md*



**1.3 DNA**

Introduction：

Requirements：

Demostration：



**1.4 DNA methylation**

Introduction：

Requirements：

Demostration：



-----



### 2. alteration calculation

Alteration calculation is based on middle or downstream output of **1.raw sequence data process**

**2.1 RNA abundance** [based on count of reads in **1.raw sequence data process**]

Introduction：RNA abundance is counting by read assigned to the all exons. To comparison betweeb samples, usually raw counts should be normalized by transcript length and library size. TPM (trancript per million read) is the common normalized form of RNA abundance. If length of transcripts are similar, such as small RNAs. Transcript length can be ignored, then CPM (count per million read) is enough for normalization.

Requirements：

```
R=3.6
R packages:
	argparse
	Biobase
	matrixStats
	dplyr
	tidyr
	biomaRt
	EnsDb.Hsapiens.v86
```

Demostration：

Before normalized abundance, raw count matrix should be prepare in the form of tab seperated file. Row is  gene identifiers and colum is sample identifiers. The following Rscript support 2 normalization method: TPM and CPM. The gene identifiers can be ensembl_gene_id, hgnc_symbol, external_gene_name and gene_names in EnsDb.Hsapiens.v86.

```
Rscript TPM.R -m raw_count_matrix.txt -n TPM -ID ensembl_gene_id -o ./output_directory
```

*More detailed information can refer to the help of TPM.R*



2.2 alternative promoter

Introduction：

Requirements：

Demostration：



2.3 allelic expression

Introduction：

Requirements：

Demostration：

2.4 chimeric RNA

Introduction：

Requirements：

Demostration：

2.5 editing

Introduction：

Requirements：

Demostration：

2.6 SNP

Introduction：

Requirements：

Demostration：

2.7 splicing

Introduction：

Requirements：

Demostration：



3. detailed analysis

   3.1 differential analysis

   3.2 gini index

   3.3 PCA analysis and distance calculation

   3.4 Cancer outlier detection

   3.5 functional enrichment analysis

   3.6 multiomics integrative pathway enrichment

   3.7 Cell type abundance estimation

   3.8 Correlation analysis