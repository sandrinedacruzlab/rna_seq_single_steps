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
assert config["end_type"] in ["pe", "se"], f"'end_type' must be one of 'pe' (paired-end) or 'se' (single_end), {config['end_type']} was passed"

# Double check can find corresponding mate files for fastq2_suffix
if config["end_type"] == "pe":
    fq2_samples = [f.replace(fastq2_suffix, "") for f in os.listdir(fastq_dir) if f.endswith(fastq2_suffix)]

    assert SAMPLES == fq2_samples, f"Inconsistent sample names between 'fastq1_suffix' & 'fastq2_suffix' (could not map mate files)"


if not os.path.exists(log_subdir):
    os.system(f"mkdir -p {log_subdir}")

if not os.path.exists(os.path.join(salmon_index_dir, salmon_index_name)):
    os.system(f"mkdir -p {os.path.join(salmon_index_dir, salmon_index_name)}")


if config['end_type'] == "pe":
    ruleorder: salmon_quant_pe > salmon_quant_se
else:
    ruleorder: salmon_quant_se > salmon_quant_pe

# Set global conda environment, avoids specifying for each rule
conda: "../envs/single_steps.yaml"

rule all:
    input:
        expand(os.path.join(out_dir, "{sample}", "quant.sf"), sample=SAMPLES)


rule custom_txome_fasta:
    '''
    Generate FASTA file of input transcripts for use with Salmon
    '''
    input:
        config["gtf"]

    output:
        os.path.join(salmon_index_dir, salmon_index_name, "transcripts.fa")

    params:
        genome_fa = config["genome_fasta"]

    log:
        os.path.join(log_subdir,
                     "custom_txome_fasta.log")

    conda: "../envs/single_steps.yaml"

    shell:
        """
        gffread \
        -w {output} \
        -g {params.genome_fa} \
        {input}
        """


def target_txome_fasta(make_fasta, custom_path):
    '''
    '''

    assert isinstance(make_fasta, bool)

    if make_fasta:
        # Need to return output of snakemake rule
        return custom_path

    else:
        return config["transcripts_fasta"]


rule generate_full_decoys:
    '''
    Generate combined FASTA of target transcripts and rest of genome
    Used to generate selective-alignment compatible FASTA file for Salmon index
    https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    '''
    input:
        genome_fa = config["genome_fasta"],
        txome_fa = target_txome_fasta(config["generate_fasta"], rules.custom_txome_fasta.output)
        # os.path.join(SALMON_SUBDIR, "min_jnc_{min_jnc}", "min_frac_{min_frac}", "min_cov_{min_cov}","papa.transcripts.fa"),

    output:
        gentrome_fa = os.path.join(salmon_index_dir, salmon_index_name, "gentrome.fa"),
        decoys = os.path.join(salmon_index_dir, salmon_index_name, "decoys.txt")

    log:
        os.path.join(log_subdir,
                     "generate_full_decoys.log")

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
        seq = os.path.join(salmon_index_dir, salmon_index_name, "seq.bin"),
        pos = os.path.join(salmon_index_dir, salmon_index_name, "pos.bin")

    params:
        k = config["salmon_kmer_size"],
        outdir = os.path.join(salmon_index_dir, salmon_index_name, "")

    threads:
        config["threads"]

    log:
        os.path.join(log_subdir,
                     "salmon_index.log")

    conda: "../envs/single_steps.yaml"

    shell:
        """
        salmon index \
        -t {input.gentrome_fa} \
        -i {params.outdir} \
        --decoys {input.decoys} \
        -k {params.k} \
        -p {threads} \
        &> {log}
        """


rule salmon_quant_pe:
    input:
        fast1 = os.path.join(fastq_dir, "{sample}" + fastq1_suffix),
        fast2 = os.path.join(fastq_dir, "{sample}" + fastq2_suffix),
        index = rules.salmon_index.output.seq

    output:
        os.path.join(out_dir, "{sample}", "quant.sf")

    params:
        index_dir = os.path.join(salmon_index_dir, salmon_index_name),
        output_dir = os.path.join(out_dir, "{sample}"),
        libtype = config["salmon_strand_info"]

    threads:
        config["threads"]

    log:
        os.path.join(log_subdir,
                     "salmon_quant_pe.{sample}.log")

    conda: "../envs/single_steps.yaml"

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --threads {threads} \
        -o {params.output_dir} \
        --seqBias \
        --posBias \
        --gcBias \
        &> {log}
        """


rule salmon_quant_se:
    input:
        fast1 = os.path.join(fastq_dir, "{sample}" + fastq1_suffix),
        index = rules.salmon_index.output.seq

    output:
        os.path.join(out_dir, "{sample}", "quant.sf")

    params:
        index_dir = os.path.join(salmon_index_dir, salmon_index_name),
        output_dir = os.path.join(out_dir, "{sample}"),
        libtype = config["salmon_strand_info"],

    threads:
        config["threads"]

    log:
        os.path.join(log_subdir,
                     "salmon_quant_se.{sample}.log")

    conda: "../envs/single_steps.yaml"

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        -r {input.fast1} \
        --threads {threads} \
        -o {params.output_dir} \
        --seqBias \
        --posBias \
        --gcBias \
        &> {log}
        """
