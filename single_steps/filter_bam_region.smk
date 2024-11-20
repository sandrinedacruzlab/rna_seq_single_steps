configfile: "config/filter_bam_region_config.yaml"

import os
import sys


in_bam_dir = os.path.join(config["input_dir"], "")
out_dir = os.path.join(config["output_dir"], "")
log_dir = os.path.join(config["output_dir"], config["log_dir"], "")

bam_suffix = config["bam_suffix"]
out_bam_suffix = config["filtered_bam_suffix"]

assert out_bam_suffix.endswith(".bam"), f"output bam suffix (config['filtered_bam_suffix']) must end in '.bam', following passed - {out_bam_suffix}"

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]

# Check that each input sample has a BAM index (.bai)
for s in SAMPLES:
    assert os.path.exists(os.path.join(in_bam_dir, s + bam_suffix + ".bai")), f".bai index file does not exist at same location as input BAM file for sample {s}"

if not os.path.exists(log_dir):
    os.makedirs(log_dir, exist_ok=True)


sys.stderr.write(f"Inferred sample names for input BAM files - {', '.join(SAMPLES)}\n")


rule all:
    input:
        expand(os.path.join(out_dir, "{sample}" + out_bam_suffix),
               sample=SAMPLES),
        expand(os.path.join(out_dir, "{sample}" + out_bam_suffix + ".bai"),
               sample=SAMPLES),



rule subset_bam:
    input:
        bam = os.path.join(in_bam_dir, "{sample}" + config["bam_suffix"]),
        bed = config["regions_bed_file"]

    output:
        bam = os.path.join(out_dir, "{sample}" + out_bam_suffix)

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
        {input.bam} 2> {log}
        """


rule index_bam:
    input:
        rules.subset_bam.output.bam

    output:
        os.path.join(out_dir, "{sample}" + out_bam_suffix + ".bai")

    log:
        os.path.join(log_dir, "{sample}.index_bam.log")

    conda:
        "../envs/single_steps.yaml"

    shell:
        """
        samtools index {input} 2> {log}
        """
