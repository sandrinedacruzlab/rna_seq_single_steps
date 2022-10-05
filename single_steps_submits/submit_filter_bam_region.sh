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

if [ "$1" != "" ]; then
    RUN_NAME=$1
else
    RUN_NAME=$""
fi

FOLDER=$(date +"%Y%m%d%H%M")
WRITEFOLDER=../submissions/$FOLDER
mkdir -p $WRITEFOLDER
cp single_steps/filter_bam_region.smk $WRITEFOLDER/filter_bam_region.smk

snakemake -s single_steps/filter_bam_region.smk \
--jobscript cluster_qsub.sh \
--cluster-config config/cluster/filter_bam_region.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $WRITEFOLDER {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
--use-conda
