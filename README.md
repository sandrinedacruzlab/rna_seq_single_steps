# Simple snakemake pipelines for common/standard bioinformatics analyses


This repo contains a collection of ad-hoc and flexible pipelines for running common bioinformatics tools on a directory of input files.

See [Available Pipelines](#Available-pipelines) section for a list & description of available workflows

# Running instructions

Perform a dry-run to check run declaration:
```
snakemake -n -p -s single_steps/<workflow>.smk
```

# Available pipelines

### Salmon (single_steps/salmon.smk)

Run [Salmon](https://github.com/COMBINE-lab/salmon) transcript quantification on a directory of input FASTQ files (paired-end or single-end). The workflow can use a pre-provided index or generate an full index with the genome sequence as decoys from provided reference files.

If generating a custom index, you have the following options:
- Provide a GTF file of transcript models - transcript sequences will be extracted using [gffread](https://github.com/gpertea/gffread)
- Provide a FASTA file of transcript sequences - will be passed straight to `salmon index`

If generating a custom index, you **must provide a genome sequence FASTA file**. To use a pre-computed index, you just need to provide the workflow with the path to the directory containing the index.

Snakefile: `single_steps/salmon.smk`

Config file: `config/salmon_config.yaml`

Cluster config file: `config/cluster/salmon.yaml`


### Pull FASTQs from BAM files (single_steps/sort_pull.smk)

Uses [samtools](https://github.com/samtools/samtools) to sort **input coordinate-sorted BAM** files by read name and extract reads to FASTQ files. Supports single-end or paired-end reads.

Note that name sorting is performed using `samtools collate` (standard mode). This ensures that reads of the same name are grouped together in the BAM file, but does not perform a full alphabetical sort (and so is quicker than `samtools sort`). This means that read pairs will be ordered randomly in the output BAM file. See [documentation](http://www.htslib.org/doc/samtools-collate.html) for a full breakdown and description of potential ramifications (for our typical RNA-seq analyses this shouldn't be an issue).

Snakefile: `single_steps/sort_pull.smk`

Config file: `config/sort_pull_config.yaml`

Cluster config file: `config/cluster/sort_pull.yaml`
