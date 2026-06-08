configfile: "config/salmon_config.yaml"

import os


fastq_dir = config["fastq_dir"]
fastq1_suffix = config["fastq1_suffix"]
fastq2_suffix = config["fastq2_suffix"]
salmon_index_dir = config["salmon_index_dir"]
salmon_index_name = config["salmon_index_name"]
out_dir = config["out_dir"]
log_subdir = os.path.join(out_dir, config["logs_dir"])

SAMPLES = [f.replace(fastq1_suffix, "") for f in os.listdir(fastq_dir) if f.endswith(fastq1_suffix)]

sys.stderr.write(f"Basenames for input FASTQ files - {', '.join(SAMPLES)}\n")


assert isinstance(config["generate_fasta"], bool), f"'generate_fasta' must be True/False boolean, {config['generate_fasta']} (type {type(config['generate_fasta'])}) was provided"
assert isinstance(config["salmon_quant_flags"], list), f"'salmon_quant_flags' must be a list, {config['salmon_quant_flags']} (type {type(config['salmon_quant_flags'])}) was provided"
assert all((isinstance(x, str) for x in config["salmon_quant_flags"])), f"All elements of 'salmon_quant_flags' must be strings"

assert config["end_type"] in ["pe", "se"], f"'end_type' must be one of 'pe' (paired-end) or 'se' (single_end), {config['end_type']} was passed"

# Double check can find corresponding mate files for fastq2_suffix
if config["end_type"] == "pe":
    fq2_samples = [f.replace(fastq2_suffix, "") for f in os.listdir(fastq_dir) if f.endswith(fastq2_suffix)]

    assert sorted(SAMPLES) == sorted(fq2_samples), f"Inconsistent sample names between 'fastq1_suffix' & 'fastq2_suffix' (could not map mate files)"


if not os.path.exists(log_subdir):
    os.system(f"mkdir -p {log_subdir}")

if not os.path.exists(os.path.join(salmon_index_dir, salmon_index_name)):
    os.system(f"mkdir -p {os.path.join(salmon_index_dir, salmon_index_name)}")


if config['end_type'] == "pe":
    ruleorder: salmon_quant_pe > salmon_quant_se
else:
    ruleorder: salmon_quant_se > salmon_quant_pe

# # Set global conda environment, avoids specifying for each rule
# conda: "../envs/single_steps.yaml"

rule all:
    input:
        expand(os.path.join(out_dir, "{sample}", "quant.sf"), sample=SAMPLES)


rule custom_txome_fasta:
    '''
    Generate FASTA file of input transcripts for use with Salmon
    '''
    input:
        gtf = config["gtf"],
        genoma_fa = config["genome_fasta"]

    output:
        fa = os.path.join(salmon_index_dir, "index_input_files", salmon_index_name, "transcripts.fa")

    log:
        stdout = os.path.join(log_subdir,
                     "custom_txome_fasta.stdout.log"),
        stderr = os.path.join(log_subdir,
                     "custom_txome_fasta.stderr.log")                  

    # conda: "../envs/single_steps.yaml"
    container: 
        "https://depot.galaxyproject.org/singularity/gffread%3A0.12.9--hf426362_0"
    # container: "docker://quay.io/biocontainers/gffread:0.9.12--0"

    shell:
        """
        gffread \
        -w {output.fa} \
        -g {input.genome_fa} \
        {input.gtf} \
        1> {log.stdout} \
        2> {log.stderr}
        """

# def target_txome_fasta(make_fasta, custom_path):
#     '''
#     '''

#     assert isinstance(make_fasta, bool)

#     if make_fasta:
#         # Need to return output of snakemake rule
#         return custom_path
#     else:
#         # Return path to provided file
#         return config["transcripts_fasta"]


rule generate_full_decoys:
    '''
    Generate combined FASTA of target transcripts and rest of genome
    Used to generate selective-alignment compatible FASTA file for Salmon index
    https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    '''
    input:
        genome_fa = config["genome_fasta"],
        txome_fa = lambda wildcards: rules.custom_txome_fasta.output.fa if config["generate_fasta"] else config["transcripts_fasta"]
        # txome_fa = target_txome_fasta(config["generate_fasta"], )
        # os.path.join(SALMON_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}","papa.transcripts.fa"),

    output:
        gentrome_fa = os.path.join(salmon_index_dir, "index_input_files", salmon_index_name,"gentrome.fa"),
        decoys = os.path.join(salmon_index_dir, "index_input_files", salmon_index_name, "decoys.txt")

    log:
        os.path.join(log_subdir,
                     "generate_full_decoys.stderr.log")

    shell:
        """
        grep "^>" {input.genome_fa} | cut -d " " -f 1 > {output.decoys} && \
        sed -i.bak -e 's/>//g' {output.decoys} && \
        cat {input.txome_fa} {input.genome_fa} > {output.gentrome_fa} \
        2> {log}
        """


rule salmon_index:
    input:
        gentrome_fa = rules.generate_full_decoys.output.gentrome_fa,
        decoys = rules.generate_full_decoys.output.decoys

    output:
        dir = directory(os.path.join(salmon_index_dir, salmon_index_name)),

    params:
        k = config["salmon_kmer_size"],
        # outdir = subpath(output.seq, parent=True)

    threads:
        1
        # config["index_threads"]

    log:
        stdout = os.path.join(log_subdir,
                     "salmon_index.stdout.log"),
                 stderr = os.path.join(log_subdir,
                     "salmon_index.stderr.log"),            

    # conda: "../envs/single_steps.yaml"

    container: "docker://quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"

    shell:
        """
        salmon index \
        -t {input.gentrome_fa} \
        -i {output.dir} \
        --decoys {input.decoys} \
        -k {params.k} \
        -p {threads} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule salmon_quant_pe:
    input:
        fast1 = os.path.join(fastq_dir, "{sample}" + fastq1_suffix),
        fast2 = os.path.join(fastq_dir, "{sample}" + fastq2_suffix),
        index = rules.salmon_index.output.dir

    output:
        sf = os.path.join(out_dir, "{sample}", "quant.sf")

    params:
        index_dir = os.path.join(salmon_index_dir, salmon_index_name),
        output_dir = subpath(output.sf, parent=True),
        libtype = config["salmon_strand_info"],
        extra_flags = " ".join(config["salmon_quant_flags"])

    threads:
        1
        # config["quant_threads"]

    log:
        stdout = os.path.join(log_subdir,
                     "salmon_quant_pe.{sample}.stdout.log"),
        stderr = os.path.join(log_subdir,
                     "salmon_quant_pe.{sample}.stderr.log")

    # conda: "../envs/single_steps.yaml"

    container: "docker://quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --threads {threads} \
        -o {params.output_dir} \
        {params.extra_flags} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule salmon_quant_se:
    input:
        fast1 = os.path.join(fastq_dir, "{sample}" + fastq1_suffix),
        index = rules.salmon_index.output.dir

    output:
        os.path.join(out_dir, "{sample}", "quant.sf")

    params:
        index_dir = os.path.join(salmon_index_dir, salmon_index_name),
        output_dir = subpath(output.sf, parent=True),
        libtype = config["salmon_strand_info"],
        extra_flags = " ".join(config["salmon_quant_flags"])

    threads:
        1
        # config["quant_threads"]

    log:
       stdout = os.path.join(log_subdir,
                     "salmon_quant_se.{sample}.stderr.log"),
       stderr = os.path.join(log_subdir,
                     "salmon_quant_se.{sample}.stderr.log")

    # conda: "../envs/single_steps.yaml"

    container: "docker://quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"


    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        -r {input.fast1} \
        --threads {threads} \
        -o {params.output_dir} \
        {params.extra_flags} \
        1> {log.stdout} \
        2> {log.stderr}
        """
