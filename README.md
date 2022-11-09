# Simple snakemake pipelines for common/standard bioinformatics analyses


This repo contains a collection of ad-hoc and flexible pipelines for running common bioinformatics tools on a directory of input files.

See [Available Pipelines](#Available-pipelines) section for a list & description of available workflows. Well defined/flexible pipelines include:
- [Salmon transcript quantification](#Salmon)
- [Pull FASTQs from coordinate-sorted BAM files](#Pull-FASTQs-from-BAM-files)
- [Generate summary statistics from BAM files with samtools stats](#Generate-summary-statistics-from-BAM-files-with-samtools-stats)
- [Infer strandedness of RNA-seq reads from BAM files with infer_experiment.py](#Infer-strandedness-of-an-RNA-seq-experiment-from-BAM-files-with-infer_experiment.py)

# Running instructions

Perform a dry-run to check run declaration:
```
snakemake -n -p -s single_steps/<workflow>.smk
```

# Available pipelines

### Salmon

Run [Salmon](https://github.com/COMBINE-lab/salmon) transcript quantification on a directory of input FASTQ files (paired-end or single-end). The workflow can use a pre-provided index or generate an full index with the genome sequence as decoys from provided reference files.

If generating a custom index, you have the following options:
- Provide a GTF file of transcript models - transcript sequences will be extracted using [gffread](https://github.com/gpertea/gffread)
- Provide a FASTA file of transcript sequences - will be passed straight to `salmon index`

If generating a custom index, you **must provide a genome sequence FASTA file**. To use a pre-computed index, you just need to provide the workflow with the path to the directory containing the index.

Snakefile: `single_steps/salmon.smk`

Config file: `config/salmon_config.yaml`

Cluster config file: `config/cluster/salmon.yaml`


### Pull FASTQs from BAM files

Uses [samtools](https://github.com/samtools/samtools) to sort **input coordinate-sorted BAM** files by read name and extract reads to FASTQ files. Supports single-end or paired-end reads.

Note that name sorting is performed using `samtools collate` (standard mode). This ensures that reads of the same name are grouped together in the BAM file, but does not perform a full alphabetical sort (and so is quicker than `samtools sort`). This means that read pairs will be ordered randomly in the output BAM file. See [documentation](http://www.htslib.org/doc/samtools-collate.html) for a full breakdown and description of potential ramifications (for our typical RNA-seq analyses this shouldn't be an issue).

Snakefile: `single_steps/sort_pull.smk`

Config file: `config/sort_pull_config.yaml`

Cluster config file: `config/cluster/sort_pull.yaml`


### Generate summary statistics from BAM files with samtools stats

Runs `samtools stats` over a set of input BAM files. Also collapses the 'summary/SN section' into a single summary table for all samples. See [documentation](http://www.htslib.org/doc/samtools-stats.html) for a full breakdown/description of calculated metrics.

Snakefile: `single_steps/samtools_stats.smk`

Config file: `config/samtools_stats_config.yaml`

Cluster config file: `config/cluster/samtools_stats.yaml`

Note: Requires pandas to be installed outside of pipeline, which is usually satisfied by a standard snakemake installation.


### Infer strandedness of an RNA-seq experiment from BAM files with infer_experiment.py

Runs RSeQC's `infer_experiment.py` over a set of input BAM files to infer the 'strandedness' of the input RNA reads. Also requires transcript annotation in BED12 format, which can be generated from a GTF file using the recipe in `single_steps/gtf_to_bed12.smk`.

See [Documentation](https://rseqc.sourceforge.net/#infer-experiment-py) for interpretation of the output files. The following blog posts are also handy for translating the definition to the correct parameter for popular RNA-seq tools - [1](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/), [2](https://littlebitofdata.com/en/2017/08/strandness_in_rnaseq/).

Snakefile: `single_steps/infer_experiment.smk`

Config file: `config/infer_experiment_config.yaml`

Cluster config file: `config/cluster/infer_experiment.yaml`

Command to submit to UCL cluster: `source submit.sh infer_experiment <optional_run_name>`
