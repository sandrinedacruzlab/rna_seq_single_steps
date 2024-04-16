configfile: "config/ribodetector_config.yaml"

assert config["end_type"] in ["pe", "se"], f"'end_type' must be one of 'pe' (paired-end) or 'se' (single_end), {config['end_type']} was passed"

fastq_dir = config["fastq_dir"]
fastq1_suffix = config["fastq1_suffix"]
fastq2_suffix = config["fastq2_suffix"]
out_dir = config["out_dir"]
log_subdir = os.path.join(out_dir, config["logs_dir"])


SAMPLES = [f.replace(fastq1_suffix, "") for f in os.listdir(fastq_dir) if f.endswith(fastq1_suffix)]

if not os.path.exists(log_subdir):
    os.system(f"mkdir -p {log_subdir}")

# Double check can find corresponding mate files for fastq2_suffix
if config["end_type"] == "pe":
    fq2_samples = [f.replace(fastq2_suffix, "") for f in os.listdir(fastq_dir) if f.endswith(fastq2_suffix)]

    assert sorted(SAMPLES) == sorted(fq2_samples), f"Inconsistent sample names between 'fastq1_suffix' & 'fastq2_suffix' (could not map mate files)"


if config['end_type'] == "pe":
    ruleorder: ribodetector_cpu_pe > ribodetector_cpu_se
else:
    ruleorder: ribodetector_cpu_se > ribodetector_cpu_pe


def get_sample_fastqs(sample, fq_dir, end_type, fq1_suffix, fq2_suffix):
    '''
    '''

    if end_type == "pe":

        fqs = [os.path.join(fq_dir, sample + fq1_suffix), os.path.join(fq_dir, sample + fq2_suffix)]

    else:
        fqs = [os.path.join(fq_dir, sample + fq1_suffix)]

    for idx, f in enumerate(fqs):
        assert os.path.exists(f), f"inferred fq not found for sample|mate - {sample} | {idx + 1} - {f}"

    return fqs


fq1_targets = expand(os.path.join(out_dir, "nonrrna", "{sample}" + config["read1_num_suffix"] + config["out_suffix_filtered"]), sample=SAMPLES)
fq2_targets = expand(os.path.join(out_dir, "nonrrna", "{sample}" + config["read2_num_suffix"] + config["out_suffix_filtered"]), sample=SAMPLES)

rule all:
    input:
        fq1_targets + fq2_targets if config["end_type"] == "pe" else fq1_targets


rule ribodetector_cpu_pe:
    input:
        lambda wildcards: get_sample_fastqs(wildcards.sample, fastq_dir, config["end_type"], fastq1_suffix, fastq2_suffix)

    output:
        # Dynamic output files based on input sample names
        nonrna_fq1=os.path.join(out_dir, "nonrrna", "{sample}" + config["read1_num_suffix"] + config["out_suffix_filtered"]),
        nonrna_fq2=os.path.join(out_dir, "nonrrna", "{sample}" + config["read2_num_suffix"] + config["out_suffix_filtered"])
    params:
        length=config["read_length"],  # Example read length
        ensure=config["ensure"],  # Classification mode
        threads=config["threads"],  # Number of threads
        chunk_size=config["chunk_size"],  # Chunk size for memory management
        log=os.path.join(log_subdir, "ribodetector.{sample}.log"),  # Log file for each sample
        rrna_output = [os.path.join(out_dir, "rrna", "{sample}" + config["read1_num_suffix"] + config["out_suffix_filtered"]), #  putting here as non-essential
                       os.path.join(out_dir, "rrna", "{sample}" + config["read2_num_suffix"] + config["out_suffix_detected"])]

    log:
        stderr=os.path.join(log_subdir, "rule_ribodetector_cpu_pe.{sample}.stderr.txt")

    container:
        "docker://quay.io/biocontainers/ribodetector:0.3.1--pyhdfd78af_0"

    shell:
        """
        ribodetector_cpu \
        -l {params.length} \
        -i {input} \
        -o {output} \
        -r {params.rrna_output} \
        -e {params.ensure} \
        -t {params.threads} \
        --chunk_size {params.chunk_size} \
        --log {params.log} \
        2> {log.stderr}
        """


rule ribodetector_cpu_se:
    input:
        lambda wildcards: get_sample_fastqs(wildcards.sample, fastq_dir, config["end_type"], fastq1_suffix, fastq2_suffix)

    output:
        nonrna_fq1=os.path.join(out_dir, "nonrrna", "{sample}" + config["read1_num_suffix"] + config["out_suffix_filtered"])
    params:
        length=config["read_length"],  # Example read length
        ensure=config["ensure"],  # Classification mode
        threads=config["threads"],  # Number of threads
        chunk_size=config["chunk_size"],  # Chunk size for memory management
        log=os.path.join(log_subdir, "ribodetector.{sample}.log"),  # Log file for each sample
        rrna_output = [os.path.join(out_dir, "rrna", "{sample}" + config["read1_num_suffix"] + config["out_suffix_filtered"]), #  putting here as non-essential
                       os.path.join(out_dir, "rrna", "{sample}" + config["read2_num_suffix"] + config["out_suffix_detected"])]

    log:
        stderr=os.path.join(log_subdir, "rule_ribodetector_cpu_se.{sample}.stderr.txt")

    container:
        "quay.io/biocontainers/ribodetector:0.3.1--pyhdfd78af_0"

    shell:
        """
        ribodetector_cpu \
        -l {params.length} \
        -i {input} \
        -o {output} \
        -r {params.rrna_output} \
        -e {params.ensure} \
        -t {params.threads} \
        --chunk_size {params.chunk_size} \
        --log {params.log} \
        2> {log.stderr}
        """

