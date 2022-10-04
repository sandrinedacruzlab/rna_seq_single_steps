configfile: "config/filter_bam_region_config.yaml"

import os
import sys


in_bam_dir = os.path.join(config["input_dir"], "")
out_dir = os.path.join(config["output_dir"], "")
log_dir = os.path.join(config["output_dir"], config["log_dir"], "")

bam_suffix = config["bam_suffix"]

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]

# Check that each input sample has a BAM index (.bai)
for s in SAMPLES:
    assert os.path.exists(os.path.join(in_bam_dir, s + bam_suffix + ".bai")), f".bai index file does not exist at same location as input BAM file for sample {s}"


rule all:
    input:
        expand(os.path.join(out_dir, "{sample}.regions.bam"),
               sample=SAMPLES),
        expand(os.path.join(out_dir, "{sample}.regions.bam.bai"),
               sample=SAMPLES),



rule subset_bam:
    input:
        bam = os.path.join(in_bam_dir, "{sample}" + config["bam_suffix"]),
        bed = config["regions_bed_file"]

    output:
        bam = os.path.join(out_dir, "{sample}.regions.bam")

    params:
        extra_threads = config["samtools_extra_threads"]

    threads:
        config["samtools_extra_threads"] + 1

    log:
        os.path.join(log_dir, "{sample}.subset_bam.log")

    conda:
        "../envs/single_steps.yaml"

    shell:
        """
        samtools view -h -b \
        -L {input.bed} \
        -o {output.bam} \
        --threads {params.extra_threads} \
        {input.bam} \

        """


rule index_bam:
    input:
        rules.subset_bam.output.bam

    output:
        os.path.join(out_dir, "{sample}.regions.bam.bai")

    log:
        os.path.join(log_dir, "{sample}.index_bam.log")

    conda:
        "../envs/single_steps.yaml"

    shell:
        """
        samtools index {input} 2> {log}
        """
