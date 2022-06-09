#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y


if [[ ( $@ == "--help") ||  $@ == "-h" ]]; then
    echo "Usage: source submit.sh SMK_NAME RUN_NAME"
    echo "SMK_NAME - Name of snakemake file (without .smk extension) under single_steps/ to run on cluster"
    echo "RUN_NAME - Optional argument to name run. Config file for run will be copied to folder containing cluster log files (.submissions/<date><time>/) with run name prefixed"
    echo "-h/--help - print this help message and exit"
    exit 0
fi

WORKFLOW="single_steps/${1}.smk"

if [ "$2" != "" ]; then
    RUN_NAME="$1"_"$2"
else
    RUN_NAME=$1
fi

FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}
cp config/${1}_config.yaml ${FOLDER}/${RUN_NAME}_${1}_config.yaml

snakemake -s ${WORKFLOW} \
--conda-prefix "/SAN/vyplab/vyplab_reference_genomes/conda_envs/" \
--use-conda \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster/${1}.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
--keep-going
