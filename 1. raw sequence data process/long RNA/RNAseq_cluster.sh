#!/bin/bash
mkdir -p cluster_log

echo "start at `date`"
snakemake --rerun-incomplete --nolock \
	-s RNAseq.snakefile -j 2 --latency-wait 60 \
	--cluster-config RNAseq_cluster.json \
	--cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.ntasks} --job-name {cluster.jobname} --partition {cluster.partition} --output {cluster.output} --error {cluster.error}" > cluster_log/0.RNAseq.log 2> cluster_log/0.RNAseq.err
echo "end at `date`"
