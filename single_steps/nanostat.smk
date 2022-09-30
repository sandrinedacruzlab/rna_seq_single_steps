configfile: "config/nanostat_config.yaml"

import os
import sys

in_bam_dir = config["input_dir"]
bam_suffix = config["bam_suffix"]
out_dir = config["output_dir"]

summary_out = os.path.join(out_dir, config["summary_out_file"])
log_dir = os.path.join(out_dir, config["log_dir"])

SAMPLES = [f.replace(bam_suffix, "")
           for f in os.listdir(in_bam_dir)
           if f.endswith(bam_suffix)]


# Check that each input sample has a BAM index (.bai)
for s in SAMPLES:
    assert os.path.exists(os.path.join(in_bam_dir, s + bam_suffix + ".bai")), f".bai index file does not exist at same location as input BAM file for sample {s}"

sys.stderr.write(f"Inferred sample names - {', '.join(SAMPLES)}\n")


wildcard_constraints:
    sample = "|".join(SAMPLES)

rule all:
    input:
        summary_out

rule nanostat:
    input:
        bam = os.path.join(in_bam_dir, "{sample}" + bam_suffix),
        bai = os.path.join(in_bam_dir, "{sample}" + bam_suffix)

    output:
        os.path.join(out_dir, "{sample}.nanostat_summary.tsv")

    log:
        os.path.join(log_dir, "{sample}_nanostat.log")

    threads:
        config["nanostat_threads"]

    container:
        "docker://quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0"

    shell:
        """
        NanoStat -n {output} \
        --tsv \
        --threads {threads} \
        --bam {input.bam} &> {log}
        """


rule combine_nanostat:
    input:
        bams = expand(os.path.join(out_dir, "{sample}.nanostat_summary.tsv"), sample=SAMPLES)

    output:
        os.path.join(out_dir, "all_samples.nanostat_summary.tsv")

    params:
        script = "scripts/combine_nanostat_tables.py",
        in_dir = out_dir,
        suffix = ".nanostat_summary.tsv"

    log:
        os.path.join(log_dir, "combine_nanostat.log")

    shell:
        """
        python {params.script} \
        {params.in_dir} \
        {params.suffix} \
        {output} &> {log}
        """
