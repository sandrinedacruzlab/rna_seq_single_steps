configfile: "config/infer_experiment_config.yaml"

import os
import sys

in_bam_dir = config["input_dir"]
bam_suffix = config["bam_suffix"]
out_dir = config["output_dir"]
log_dir = os.path.join(out_dir, config["log_dir"])

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]

# Check that each input sample has a BAM index (.bai)
for s in SAMPLES:
    assert os.path.exists(os.path.join(in_bam_dir, s + bam_suffix + ".bai")), f".bai index file does not exist at same location as input BAM file for sample {s}"


wildcard_constraints:
    sample = "|".join(SAMPLES)


rule all:
    input:
        expand(os.path.join(out_dir, "{sample}.infer_experiment.txt"), sample=SAMPLES)


rule infer_experiment:
    input:
        bam = os.path.join(in_bam_dir, "{sample}" + bam_suffix),
        idx = os.path.join(in_bam_dir, "{sample}" + bam_suffix + ".bai")

    output:
        os.path.join(out_dir, "{sample}.infer_experiment.txt")

    params:
        annotation = config["bed12"],
        sample_size = config["sample_size"], # Number of reads to sample from BAM
        min_qual = config["min_qual"] # min map qual to be considered 'uniquely mapped'

    conda:
        "../envs/single_steps.yaml"

    log:
        os.path.join(log_dir, "{sample}.infer_experiment.log")

    shell:
        """
        infer_experiment.py \
        -i {input.bam} \
        -r {params.annotation} \
        -s {params.sample_size} \
        -q {params.min_qual} > {output}
        """
