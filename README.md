# Analysis-pipeline-for-multiomics
### 1.Raw sequence data process

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

**1.3 cfDNA Methylation/cfMeDIP (DIP-seq)**

**Requirements**

- smk: ./1.raw sequence data process/DIP-seq/snakemake/DIP-seq-pe.snakemake
- cfg: ./1.raw sequence data process/DIP-seq/config/test.yaml
- env:
  - ./1.raw sequence data process/DIP-seq/snakemake/envs/DIP.yml
  - ./1.raw sequence data process/DIP-seq/snakemake/envs/py27.yml
  - ./1.raw sequence data process/DIP-seq/snakemake/envs/r35.yml

**Demostration：**
```
dst="test"
snakemake --rerun-incomplete --keep-going --printshellcmds --reason \
  --use-conda --nolock --latency-wait 200 --restart-times 1 --jobs 16 \
  --snakefile snakemake/DIP-seq-pe.snakemake \
  --configfile config/${dst}.yaml \

  > run-${dst}_methylation.log 2>&1
```

**1.4 cfDNA WGS (DNA-seq)**

**Requirements**

- smk:
  - ./1.raw sequence data process/DNA-seq/snakemake/DNA-seq-common.snakemake
  - ./1.raw sequence data process/DNA-seq/snakemake/DNA-seq-pe.snakemake
- cfg: ./1.raw sequence data process/DNA-seq/config/test.yaml
- env:
  - ./1.raw sequence data process/DNA-seq/snakemake/envs/DNA.yml
    - ./1.raw sequence data process/DIP-seq/snakemake/envs/py27.yml
    - ./1.raw sequence data process/DIP-seq/snakemake/envs/r35.yml

**Demostration：**

```
dst="test"
snakemake --rerun-incomplete --keep-going --printshellcmds --reason \
  --use-conda --nolock --latency-wait 200 --restart-times 1 --jobs 16 \
  --snakefile snakemake/DNA-seq-pe.snakemake \
  --configfile config/${dst}.yaml \

  > run-${dst}_WGS.log 2>&1
```


### 2. Alteration calculation

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

Demonstration：

Before normalized abundance, raw count matrix should be prepare in the form of tab seperated file. Row is  gene identifiers and colum is sample identifiers. The following Rscript support 2 normalization method: TPM and CPM. The gene identifiers can be ensembl_gene_id, hgnc_symbol, external_gene_name and gene_names in EnsDb.Hsapiens.v86.

```
Rscript TPM.R -m raw_count_matrix.txt -n TPM -ID ensembl_gene_id -o ./output_directory
```

*More detailed information can refer to the help of TPM.R*



**2.2 alternative promoter**  [based on adaptor cutted sequences in **1.raw sequence data process**]

Introduction：transcript isoform abundance was quantified by *salmon* and normalized to TPM. TPMs of isoforms with transcript start sites within 10 bp (sharing the same promoter) were aggregated as one promoter activity. TPM < 1 promoter is defined as inactive promoter. Promoter with the highest relative promoter activity is defined as major promoter. The remaining promoters are defined as minor promoter

Requirements：

```
python=3.7.4
R=3.6
R packages:
	argparse
	dplyr
	GenomicFeatures
	doParallel
	parallel

salmon=1.4.0
```

Demonstration：

Step A. generate normalized promoter usage.

Prepare adaptor cutted fastq files in ./${dataset_name}/output/cutadapt, and reference files in ${ref}.

```
sh level_3_Alt.promoter.sh ${dataset_name}
```

*Details in help of salmon, prepareTxQuantMat.py, getPromoterActivity.py, AlternativePromoter_Normalize.R*

Step B. classify promoters

If we need to classify promoter into inactive, minor and major promoters, we can run the following Rscript.

```
Rscript Promoter_usage.R -n ${normalized promoter usage matrix file from step A}
```

*Details in help of Promoter_usage.R*



**2.3 allelic expression** [based on non-duplicated sequences mapped to genome in **1.raw sequence data process**]

Introduction：*GATK ASEReadCounter* were used to identify allele specific expression gene site based on SNP sites. ﻿For each individual, ﻿Allelic expression (AE, AE = |0.5 − Reference ratio|, Reference ratio = Reference reads/Total reads) was calculated for all sites with ≥16 reads

Requirements：

```
python=3.7.4
R=3.6
R packages:
	argparse
	dplyr
	doParallel

gatk-4.1.9.0
```



Demonstration：

Step A. generate per-sample allelic alter count at each SNP site

Prepare non-duplicated sequences mapped to genome by 2passMode in STAR in ./${dataset_name}/output/bam, and reference files in ${ref}.

```
sh level_3_Editing_SNP_ASE.sh ${dataset_name}
```

Step B. Summarize per-sample allelic count and reference count to Allelic expression matrix

 ```
Rscript ASE_AE.R -i ${output}/ASE/COSMIC -o ${output}/matrix/ASE_COSMIC_${dataset}.txt
 ```

*Details in help of ASE_AE.R*



**2.4 chimeric RNA** [based on sequences unmapped to genome in **1.raw sequence data process**]

Introduction：Reads unaligned to genome were remapped to chimeric junctions by *STAR-fusion* to identify chimeric RNA.

Requirements：

```
Python=3.7.0
R=3.5.1
R packages:
	argparse
	dplyr
	doParallel
STAR=2.5.3a_modified [2.5.3a is ok]
STAR-Fusion=1.10.0
```

*Detailed environment please refer to ChimericRNA.yaml*


Demonstration：

Prepare sequences unmapped to genome in ./${dataset_name}/output/unmapped, and reference files in ${ref}.

```
sh level_3_star-fusion.sh ${dataset_name}
```

*Details in help of STAR-Fusion and STARFusion_matrix.R*



**2.5 editing** [based on non-duplicated sequences mapped to genome in **1.raw sequence data process**]

Introduction：*GATK ASEReadCounter* were used to identify editing sites based on REDIportal. The editing ratio was defined as allele count divided by total count.

Requirements：

```
python=3.7.4
R=3.6
R packages:
	argparse
	dplyr
	doParallel

gatk-4.1.9.0
```



Demostration：

Step A. generate per-sample allelic alter count at each RNA editing site

Prepare non-duplicated sequences mapped to genome by 2passMode in STAR in ./${dataset_name}/output/bam, and reference files in ${ref}.

```
sh level_3_Editing_SNP_ASE.sh ${dataset_name}
```

Step B. Summarize per-sample allelic count and reference count to RNA editing site matrix

 ```
Rscript Editing_ratio.R -i ${output}/REDIportal -o ${output}/matrix/Editing_ratio.txt
 ```

*Details in help of Editing_ratio.R*



**2.6 SNP** [based on non-duplicated sequences mapped to genome in **1.raw sequence data process**]

Introduction：intron-spanning reads were splited by *GATK SplitNCigarReads* for confident SNP calling at RNA level. Then, alterations were identified by *GATK HaplotypeCaller* and filtered by *GATK VariantFilteration* with the following 4 criteria: strand bias defined by fisher exact test phred-scaled p value (FS) < 20, variant confidence (QUAL) divided by the unfiltered depth (QD) > 2, total number of reads at the variant site (DP) > 10, SNP quality (QUAL) > 20. Allele fraction was defined as allele count divided by total count (reference count and allele count).

Requirements：

```
python=3.7.4
R=3.6
R packages:
	argparse
	dplyr
	doParallel

gatk-4.1.9.0
ensembl-vep=104.3
```



Demonstration：

Step A. generate per-sample allelic alter count at each RNA editing site

Prepare non-duplicated sequences mapped to genome by 2passMode in STAR in ./${dataset_name}/output/bam, and reference files in ${ref}.

```
sh level_3_Editing_SNP_ASE.sh ${dataset_name}
```

Step B. summarize per-sample allelic count and reference count to RNA editing site matrix

 ```
Rscript SNP_AF.R -i ${output}/SNP -o ${output}/matrix
 ```

*Details in help of SNP_AF.R*

Step C. annotate SNPs by VEP

```
#Firstly filter low quality SNP in Step.B
sh vcf_filter.sh ${dataset_name}

#Then annotate high qualilty SNP by VEP
sh VEP.sh ${dataset_name}

#mutation classification can be summaried
Summary_mutation_class_forplot.R  -i ${directory with VEP annotated files} -o ${output}
```

*Details in help of VEP, SNP_AF.R and Summaru_mutation_class_forplot.R*



**2.7 splicing** [based on non-duplicated sequences mapped to genome in **1.raw sequence data process**]

Introduction：The percent spliced-in (PSI) score of each alternative splicing event was calculated using *rMATs-turbo*.

Requirements：

```
python=2.7.12
R=3.6
rmats=4.1.2
```

Demonstration：

Prepare positive and negative group of sample identifers and their full path to non-duplicated sequences mapped to genome by 2passMode in STAR.

```
sh level3_splicing.sh ${dataset_name} ${positive} ${negative}
```

*Details in help of rmats.py*

**2.8 methylation level**

included within ./1.raw sequence data process/DIP-seq/snakemake/DIP-seq-pe.snakemake (rule: count_matrix_gene)

**2.9 cfDNA WGS CNV**

Requirements：

- smk: ./1.raw sequence data process/DNA-seq/snakemake/DNA-seq-pe-CNV.snakemake
- cfg: ./1.raw sequence data process/DNA-seq/config/test.yaml
- env:
  - ./1.raw sequence data process/DNA-seq/snakemake/envs/DNA_CNV_WisecondorX.yml
  - ./1.raw sequence data process/DNA-seq/snakemake/envs/DNA_CNV_CNVkit.yml

Demonstration：
included within pre-process (snakemake/DNA-seq-pe-CNV.snakemake)


**2.10 cfDNA WGS FragSize**

Requirements：

- smk: ./1.raw sequence data process/DNA-seq/snakemake/DNA-seq-pe-FragSize.snakemake
- cfg: ./1.raw sequence data process/DNA-seq/config/test.yaml
- env:
  - ./1.raw sequence data process/DNA-seq/snakemake/envs/DNA_FragSize.yml

Demonstration：

included within pre-process (./1.raw sequence data process/DNA-seq/snakemake/DNA-seq-pe-FragSize.snakemake)


**2.11 cfDNA WGS WPS**

(modified from https://github.com/kircherlab/cfDNA)

Requirements：

- smk: ./1.raw sequence data process/DNA-seq/snakemake/WPS/snakefile_WPS.smk
- cfg: ./1.raw sequence data process/DNA-seq/snakemake/WPS/config/config.yml
- env: ./1.raw sequence data process/DNA-seq/snakemake/WPS/workflow/envs/cfDNA.yml

Demonstration：

included within pre-process (./1.raw sequence data process/DNA-seq/WPS/snakefile_WPS.smk)

### 3. Detailed analysis

Documented in ./3.detailed analysis/Multiomics_final.R
