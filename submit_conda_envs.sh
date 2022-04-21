#!/bin/bash
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=6G,h_rt=20:00:00,tmem=6G

# join stdout and stderr output
#$ -j y
#$ -R y
#$ -N single_steps_conda_envs

snakemake \
-p \
-s $1 \
--conda-prefix "/SAN/vyplab/vyplab_reference_genomes/conda_envs/" \
--use-conda \
--conda-create-envs-only \
--conda-frontend mamba \
--cores 1
