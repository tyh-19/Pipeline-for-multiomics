######################################################################
######################################################################
#config
tool_dir="/data/taoyuhuan/RNAseq/tools"

#reference
genome_dir="/data/taoyuhuan/RNAseq/reference/genome_index/star"
fasta_dir="/data/taoyuhuan/RNAseq/reference/fasta"
bed_dir="/data/taoyuhuan/RNAseq/reference/bed"
gtf_dir="/data/taoyuhuan/RNAseq/reference/gtf"

#input
data_dir="/data/taoyuhuan/RNAseq/test_RNAseq/fastq"
sample_ids=["CRC-2418488_test_neg","CRC-2276341_test_neg"]

#output
output_dir="/data/taoyuhuan/RNAseq/test_RNAseq/output"

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
temp_dir="/data/taoyuhuan/RNAseq/test_RNAseq/temp"

######################################################################
######################################################################

#pipeline
def get_all_inputs(wildcards):
	available_inputs = dict(
		cutadapt=expand('{output_dir}/cutadapt/{sample_id}_{mate_index}.fastq.gz',
			output_dir=output_dir, sample_id=sample_ids, mate_index=[1, 2]),
		trimGC=expand('{output_dir}/trimGC/{sample_id}_{mate_index}.fastq.gz',
			output_dir=output_dir, sample_id=sample_ids, mate_index=[1, 2]),
		bam=expand('{output_dir}/bam/{sample_id}/{map_step}.bam',
			output_dir=output_dir, sample_id=sample_ids, map_step=map_steps),
		bam_sort_by_name=expand('{output_dir}/bam/{sample_id}/{map_step_sortbyName}_sortbyName.bam',
			output_dir=output_dir, sample_id=sample_ids, map_step_sortbyName=map_steps_sortbyName),
		bam_rmdup=expand('{output_dir}/bam/{sample_id}/{map_step_rmdup}_rmdup.bam',
			output_dir=output_dir, sample_id=sample_ids, map_step_rmdup=map_steps_rmdup),
		bam_rumdup_RG=expand('{output_dir}/bam/{sample_id}/hg38_long_rmdup_RG.bam',
			output_dir=output_dir, sample_id=sample_ids),
		bam_intron_spanning=expand('{output_dir}/bam/{sample_id}/hg38_long_intron_spanning.bam',
			output_dir=output_dir, sample_id=sample_ids),
		unmapped=expand('{output_dir}/unmapped/{sample_id}/{map_step}_{mate_index}.fastq.gz',
			output_dir=output_dir, sample_id=sample_ids, map_step=map_steps, mate_index=[1, 2]),
		count=expand('{output_dir}/counts/hg38_long/{sample_id}/featurecount',
			output_dir=output_dir, sample_id=sample_ids),
		count_rmdup=expand('{output_dir}/counts/hg38_long_rmdup/{sample_id}/featurecount',
			output_dir=output_dir, sample_id=sample_ids),
		count_intron_spanning=expand('{output_dir}/counts/intron_spanning/{sample_id}/featurecount',
			output_dir=output_dir, sample_id=sample_ids),
		count_summary=expand('{output_dir}/counts/hg38_long/{sample_id}/featurecount.summary',
			output_dir=output_dir, sample_id=sample_ids),
		count_rmdup_summary=expand('{output_dir}/counts/hg38_long_rmdup/{sample_id}/featurecount.summary',
			output_dir=output_dir, sample_id=sample_ids),
		count_intron_spanning_summary=expand('{output_dir}/counts/intron_spanning/{sample_id}/featurecount.summary',
			output_dir=output_dir, sample_id=sample_ids),
#		ref_gene_table=expand('{output_dir}/log/count_matrix/reference_genes.txt',
#			output_dir=output_dir),
		count_matrix=expand('{output_dir}/count_matrix/{count_level}',
			output_dir=output_dir,count_level=count_levels),
		fastq_read_pair=expand('{output_dir}/log/{sample_id}/read_pair/{clean_level}.txt',
			output_dir=output_dir, sample_id=sample_ids,clean_level=clean_levels),
		raw_read_pair=expand('{output_dir}/log/{sample_id}/read_pair/raw.txt',
			output_dir=output_dir, sample_id=sample_ids),
		samtools_stats=expand('{output_dir}/log/{sample_id}/samtool_stats/{map_step}.txt',
			output_dir=output_dir, sample_id=sample_ids, map_step=map_steps),
		samtools_stats_rmdup=expand('{output_dir}/log/{sample_id}/samtool_stats/{map_step_rmdup}_rmdup.txt',
			output_dir=output_dir, sample_id=sample_ids, map_step_rmdup=map_steps_rmdup),
		samtools_stats_intron_spanning=expand('{output_dir}/log/{sample_id}/samtool_stats/hg38_long_intron_spanning.txt',
			output_dir=output_dir, sample_id=sample_ids),
		bam_read_pair=expand('{output_dir}/log/{sample_id}/read_pair/{map_step}.txt',
			output_dir=output_dir, sample_id=sample_ids, map_step=map_steps),
		bam_rmdup_read_pair=expand('{output_dir}/log/{sample_id}/read_pair/{map_step_rmdup}_rmdup.txt',
			output_dir=output_dir, sample_id=sample_ids, map_step_rmdup=map_steps_rmdup),
		bam_intron_spanning_read_pair=expand('{output_dir}/log/{sample_id}/read_pair/hg38_long_intron_spanning.txt',
			output_dir=output_dir, sample_id=sample_ids),
		genomic_region=expand("{output_dir}/log/{sample_id}/genomic_region",
			output_dir=output_dir, sample_id=sample_ids)
	)
	enabled_inputs = list(available_inputs.keys())
	inputs = []
	for key, l in available_inputs.items():
		if key in enabled_inputs:
			inputs += l
	return inputs

rule all:
	input:
		get_all_inputs

rule cutadapt_pe:
	input:
		fastq1=data_dir+'/{sample_id}_1.fastq.gz',
		fastq2=data_dir+'/{sample_id}_2.fastq.gz'
	output:
		fastq1='{output_dir}/cutadapt/{sample_id}_1.fastq.gz',
		fastq2='{output_dir}/cutadapt/{sample_id}_2.fastq.gz'
	log:
		'{output_dir}/log/{sample_id}/cutadaptor.log'
	shell:
		'''
		cutadapt --pair-filter any -j {thread_cutadapt} -q 30,30 -a {adaptor1} -A {adaptor2} \
		--trim-n -m 16 -o {output.fastq1} -p {output.fastq2} \
		{input.fastq1} {input.fastq2} > {log} 2>&1
		'''

rule trimGC:
	input:
		fastq1='{output_dir}/cutadapt/{sample_id}_1.fastq.gz',
		fastq2='{output_dir}/cutadapt/{sample_id}_2.fastq.gz'
	output:
		fastq1='{output_dir}/trimGC/{sample_id}_1.fastq.gz',
		fastq2='{output_dir}/trimGC/{sample_id}_2.fastq.gz'
	log:
		'{output_dir}/log/{sample_id}/trimGC.log'
	shell:
		'''
		python {tool_dir}/trimGC.py -m 16 -s {strandness} -o {output_dir}/trimGC/{wildcards.sample_id} -i {output_dir}/cutadapt/{wildcards.sample_id} > {log} 2>&1
		'''

map_command_pe = '''STAR --genomeDir {params.index} \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --runThreadN {thread_mapping} \
            --outFileNamePrefix {params.output_prefix} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand gzip -d -c \
            --outSAMmultNmax 1 \
            --seedPerWindowNmax {params.seedPerWindowNmax} > {log} 2>&1
        mv {params.output_prefix}Aligned.out.bam {output.bam}
        {tool_dir}/bbmap/repair.sh in={params.output_prefix}Unmapped.out.mate1 in2={params.output_prefix}Unmapped.out.mate2 out={output.unmapped1} out2={output.unmapped2} >> {log} 2>&1
        rm -f {params.output_prefix}Unmapped.out.mate1 {params.output_prefix}Unmapped.out.mate2
        '''


rule map_spikein_long:
	input:
		fastq1='{output_dir}/trimGC/{sample_id}_1.fastq.gz',
		fastq2='{output_dir}/trimGC/{sample_id}_2.fastq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/spikein_long.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/spikein_long_1.fastq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/spikein_long_2.fastq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/spikein_long/mapping.log'
	params:
		index=genome_dir+'/spikein_long',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/spikein_long/',
		seedPerWindowNmax=20
	run:
		shell(map_command_pe)

rule map_univec:
	input:
		fastq1='{output_dir}/unmapped/{sample_id}/spikein_long_1.fastq.gz',
		fastq2='{output_dir}/unmapped/{sample_id}/spikein_long_2.fastq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/univec.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/univec_1.fastq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/univec_2.fastq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/univec/mapping.log'
	params:
		index=genome_dir+'/univec',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/univec/',
		seedPerWindowNmax=20
	run:
                shell(map_command_pe)

rule map_rRNA:
	input:
		fastq1='{output_dir}/unmapped/{sample_id}/univec_1.fastq.gz',
		fastq2='{output_dir}/unmapped/{sample_id}/univec_2.fastq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/rRNA.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/rRNA_1.fastq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/rRNA_2.fastq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/rRNA/mapping.out'
	params:
		index=genome_dir+'/rRNA',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/rRNA/',
		seedPerWindowNmax=20
	run:
		shell(map_command_pe)

rule map_hg38_long:
	input:
		fastq1='{output_dir}/unmapped/{sample_id}/rRNA_1.fastq.gz',
		fastq2='{output_dir}/unmapped/{sample_id}/rRNA_2.fastq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/hg38_long.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/hg38_long_1.fastq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/hg38_long_2.fastq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/hg38_long/mapping.log'
	params:
		index=genome_dir+'/hg38_long',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/hg38_long/',
		seedPerWindowNmax=50
	run:
		shell(map_command_pe)

rule map_circRNA:
	input:
		fastq1='{output_dir}/unmapped/{sample_id}/hg38_long_1.fastq.gz',
		fastq2='{output_dir}/unmapped/{sample_id}/hg38_long_2.fastq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/circRNA.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/circRNA_1.fastq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/circRNA_2.fastq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/circRNA/mapping.log'
	params:
		index=genome_dir+'/circRNA',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/circRNA/',
		seedPerWindowNmax=20
	run:
		shell(map_command_pe)


rule sort_bam_by_name:
	input:
		'{output_dir}/bam/{sample_id}/{map_step_sortbyName}.bam'
	output:
		'{output_dir}/bam/{sample_id}/{map_step_sortbyName}_sortbyName.bam'
	shell:
		'''
		samtools sort -n -T {temp_dir} -o {output} {input}
		'''

rule remove_duplicates:
	input:
		'{output_dir}/bam/{sample_id}/{map_step_sortbyName}_sortbyName.bam'
	output:
		bam_rmdup='{output_dir}/bam/{sample_id}/{map_step_sortbyName}_rmdup.bam',
		metrics='{output_dir}/log/{sample_id}/remove_duplicate/{map_step_sortbyName}_rmdup'
	log:
		'{output_dir}/log/{sample_id}/picard_rmdup/{map_step_sortbyName}_remove_duplicates.log'
	shell:
		'''
		java -jar {tool_dir}/picard/picard.jar \
		MarkDuplicates REMOVE_DUPLICATES=true \
		ASSUME_SORT_ORDER=queryname \
		I={input} \
		O={output.bam_rmdup} \
		M={output.metrics} \
		READ_NAME_REGEX=null > {log} 2>&1
		'''

rule add_ReadGroup:
	input:
		'{output_dir}/bam/{sample_id}/hg38_long_rmdup.bam'
	output:
		'{output_dir}/bam/{sample_id}/hg38_long_rmdup_RG.bam'
	log:
		'{output_dir}/log/{sample_id}/gatk/add_readgroup.log'
	params:
		sample='{sample_id}'
	shell:
		'''
		{tool_dir}/gatk-4.1.9.0/gatk AddOrReplaceReadGroups --java-options -Xmx4G \
		--INPUT {input} \
		--OUTPUT {output} \
		--TMP_DIR {temp_dir} \
		-SO queryname \
		--RGLB library \
		--RGPL illumina \
		--RGPU HiSeq2000 \
		--RGSM {params.sample} > {log} 2>&1
		'''

rule get_Intron_spanning:
	input:
		'{output_dir}/bam/{sample_id}/hg38_long_rmdup_RG.bam'
	output:
		'{output_dir}/bam/{sample_id}/hg38_long_intron_spanning.bam'
	log:
		'{output_dir}/log/{sample_id}/gatk/get_intron_spanning.log'
	shell:
		'''
		{tool_dir}/gatk-4.1.9.0/gatk SplitNCigarReads --java-options -Xmx4G \
		--input {input} \
		--output {output} \
		-R {fasta_dir}/hg38.fa > {log} 2>&1
		'''

rule featurecounts:
	input:
		bam='{output_dir}/bam/{sample_id}/hg38_long_sortbyName.bam',
		gtf=gtf_dir+'/long_RNA.gtf'
	output:
		counts='{output_dir}/counts/hg38_long/{sample_id}/featurecount',
		summary='{output_dir}/counts/hg38_long/{sample_id}/featurecount.summary'
	params:
		strandness={'no': 0, 'forward': 1, 'reverse': 2}[strandness],
		paired_end={'True': '-p', 'False': ''}[paired_end],
		min_mapping_quality={min_mapping_quality},
		count_multimap_reads={'True': '-M','False': ''}[count_multimap_reads],
		count_overlapping_features={'True': '-O','False': ''}[count_overlapping_features]
	log:
		'{output_dir}/log/{sample_id}/featurecount/hg38_long.log'
	shell:
		'''
		featureCounts {params.count_overlapping_features} -t exon -g gene_id {params.count_multimap_reads} \
		-s {params.strandness} -Q {params.min_mapping_quality} \
		{params.paired_end} -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
		'''

rule featurecounts_rmdup:
	input:
		bam='{output_dir}/bam/{sample_id}/hg38_long_rmdup.bam',
		gtf=gtf_dir+'/long_RNA.gtf'
	output:
		counts='{output_dir}/counts/hg38_long_rmdup/{sample_id}/featurecount',
		summary='{output_dir}/counts/hg38_long_rmdup/{sample_id}/featurecount.summary'
	params:
		strandness={'no': 0, 'forward': 1, 'reverse': 2}[strandness],
		paired_end={'True': '-p', 'False': ''}[paired_end],
		min_mapping_quality={min_mapping_quality},
		count_multimap_reads={'True': '-M','False': ''}[count_multimap_reads],
		count_overlapping_features={'True': '-O','False': ''}[count_overlapping_features]
	log:
		'{output_dir}/log/{sample_id}/featurecount/hg38_long_rmdup.log'
	shell:
		'''
		featureCounts {params.count_overlapping_features} -t exon -g gene_id {params.count_multimap_reads} \
		-s {params.strandness} -Q {params.min_mapping_quality} \
		{params.paired_end} -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
		'''

rule featurecounts_intron_spanning:
	input:
		bam='{output_dir}/bam/{sample_id}/hg38_long_intron_spanning.bam',
		gtf=gtf_dir+'/long_RNA.gtf'
	output:
		counts='{output_dir}/counts/intron_spanning/{sample_id}/featurecount',
		summary='{output_dir}/counts/intron_spanning/{sample_id}/featurecount.summary'
	params:
		strandness={'no': 0, 'forward': 1, 'reverse': 2}[strandness],
		paired_end={'True': '-p', 'False': ''}[paired_end],
		min_mapping_quality={min_mapping_quality},
		count_multimap_reads={'True': '-M','False': ''}[count_multimap_reads],
		count_overlapping_features={'True': '-O','False': ''}[count_overlapping_features]
	log:
		'{output_dir}/log/{sample_id}/featurecount/hg38_long_intron_spanning.log'
	shell:
		'''
		featureCounts {params.count_overlapping_features} -t exon -g gene_id {params.count_multimap_reads} \
		-s {params.strandness} -Q {params.min_mapping_quality} \
		{params.paired_end} -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
		'''

#rule ref_gene_id_type_name:
#	input:
#		gtf=gtf_dir+'/long_RNA.gtf'
#	output:
#		'/data/taoyuhuan/RNAseq/reference/gtf/reference_genes.txt'
#	shell:
#		'''
#		awk '{for (i = 1; i<= NF; i++) {if ($i ~/gene_id|gene_type|gene_name/) {printf "%s ",$(i+1)}} print ""}' {input.gtf} | sed -e 's/"//g' -e 's/;//g' -e 's/ /|/g' | sort -k1,1 | uniq > {output}
#		sed -i "1d" {output}
#		sed -i "1iensembl_id|gene_type|gene_name" {output}
#		sed -i "s/|$//" {output}
#		'''

rule count_matrix:
	input:
		ref_genes='/data/taoyuhuan/RNAseq/reference/gtf/reference_genes.txt'
	output:
		count_matrix=directory('{output_dir}/count_matrix/{count_level}')
	params:
		input='{output_dir}/counts/{count_level}',
		counts='{output_dir}/log/count_matrix/counts_{count_level}',
		rownames='{output_dir}/log/count_matrix/rownames_{count_level}'
	shell:
		'''
		sh {tool_dir}/count_matrix.sh {params.input} {params.rownames} {params.counts} {input.ref_genes} {output.count_matrix}
		'''


rule samtools_stats:
	input:
		'{output_dir}/bam/{sample_id}/{map_step}.bam'
	output:
		'{output_dir}/log/{sample_id}/samtool_stats/{map_step}.txt'
	wildcard_constraints:
		map_step='(spikein_long)|(univec)|(rRNA)|(hg38_long)|(circRNA)'
	shell:
		'''
		samtools stats {input} > {output}
		'''

rule samtools_stats_rmdup:
	input:
		'{output_dir}/bam/{sample_id}/{map_step_rmdup}_rmdup.bam'
	output:
		'{output_dir}/log/{sample_id}/samtool_stats/{map_step_rmdup}_rmdup.txt'
	shell:
		'''
		samtools stats {input} > {output}
		'''

rule samtools_stats_intron_spanning:
	input:
		'{output_dir}/bam/{sample_id}/hg38_long_intron_spanning.bam'
	output:
		'{output_dir}/log/{sample_id}/samtool_stats/hg38_long_intron_spanning.txt'
	shell:
		'''
		samtools stats {input} > {output}
		'''

rule count_raw_fastq_read_pairs:
	input:
		raw=data_dir+'/{sample_id}_1.fastq.gz'
	output:
		'{output_dir}/log/{sample_id}/read_pair/raw.txt'
	shell:
		'''
		pigz -p {threads_decompress} -d -c {input.raw} | wc -l | awk '{{print int($0/4)}}' > {output}
		'''

rule count_clean_fastq_read_pairs:
	input:
		'{output_dir}/{clean_level}/{sample_id}_1.fastq.gz'
	output:
		'{output_dir}/log/{sample_id}/read_pair/{clean_level}.txt'
	wildcard_constraints:
		clean_level='(cutadapt)|(trimGC)'
	shell:
		'''
		pigz -p {threads_decompress} -d -c {input} | wc -l | awk '{{print int($0/4)}}' > {output}
		'''

rule count_bam_read_pairs:
	input:
		'{output_dir}/log/{sample_id}/samtool_stats/{map_step}.txt'
	output:
		'{output_dir}/log/{sample_id}/read_pair/{map_step}.txt'
	wildcard_constraints:
		map_step='(spikein_long)|(univec)|(rRNA)|(hg38_long)|(circRNA)'
	shell:
		'''
		awk 'BEGIN{{OFS="\t";FS="\t"}}/^SN/{{if($2 == "reads mapped and paired:") print int($3/2)}}' {input} > {output}
		'''

rule count_bam_rmdup_read_pairs:
	input:
		'{output_dir}/log/{sample_id}/samtool_stats/{map_step_rmdup}_rmdup.txt'
	output:
		'{output_dir}/log/{sample_id}/read_pair/{map_step_rmdup}_rmdup.txt'
	shell:
		'''
		awk 'BEGIN{{OFS="\t";FS="\t"}}/^SN/{{if($2 == "reads mapped and paired:") print int($3/2)}}' {input} > {output}
		'''

rule count_bam_intron_spanning_read_pairs:
	input:
		'{output_dir}/log/{sample_id}/samtool_stats/hg38_long_intron_spanning.txt'
	output:
		'{output_dir}/log/{sample_id}/read_pair/hg38_long_intron_spanning.txt'
	shell:
		'''
		awk 'BEGIN{{OFS="\t";FS="\t"}}/^SN/{{if($2 == "reads mapped and paired:") print int($3/2)}}' {input} > {output}
		'''

rule genomic_region:
	input:
		'{output_dir}/bam/{sample_id}/hg38_long.bam'
	output:
		directory('{output_dir}/log/{sample_id}/genomic_region')
	log:
		'{output_dir}/log/{sample_id}/genomic_region/genomic_region.log'
	shell:
		'''
		sh {tool_dir}/sequential.assign.long.sh {input} {output} {bed_dir} > {log} 2>&1
		'''


