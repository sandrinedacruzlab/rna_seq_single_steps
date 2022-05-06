configfile: "config/sort_pull_config.yaml"
import os


#########---------------
## Input parameters
#########---------------
config_path = "config/sort_pull_config.yaml" # this is for copying later

in_bam_dir = config["input_dir"]
out_temp_dir = config["temp_dir"]
out_fastq_dir = config["fastq_spot"]

bam_suffix = config["bam_suffix"]
end_type = config["end_type"]
assert end_type in ["pe", "se"], f"'end_type' must be one of 'pe' (paired-end) or 'se' (single-end) - {end_type} was provided"

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]

sys.stderr.write(f"Basenames for input BAM files - {', '.join(SAMPLES)}\n")
# sys.stderr.write("".join(str(s) for s in SAMPLES) + "\n")


if not os.path.exists(out_fastq_dir):
    os.system(f"mkdir -p {out_fastq_dir}")

if not os.path.exists(out_temp_dir):
    os.system(f"mkdir -p {out_temp_dir}")

# conda: "../envs/single_steps.yaml"

localrules: copy_config

rule all:
  input:
    expand(os.path.join(out_fastq_dir, "{sample}.1.pulled.fastq.gz"),
           sample = SAMPLES),
    expand(os.path.join(out_fastq_dir, "{sample}.2.pulled.fastq.gz"),
           sample = SAMPLES if end_type == "pe" else []),
    os.path.join(out_fastq_dir, "config_sort_pull.yaml")


rule name_sort:
    '''
    samtools collate puts reads of same name/pair next to eachother in the sorted BAM, but doesn't care about sorting read names alphabetically within the whole file
    This is quicker than a full name sort and is fine as input to samtools fastq
    Unlikely to be necessary for pe read names to be alphabetically sorted (just that mates occur at same position in each file)
    '''
    input:
        aligned_bam = os.path.join(in_bam_dir, "{sample}" + bam_suffix)

    output:
        temp(os.path.join(out_temp_dir, "{sample}.namesorted.bam"))

    threads:
        config["collate_extra_threads"] + 1

    conda:
        "../envs/single_steps.yaml"

    shell:
        """
        samtools collate \
        -@ {threads} \
        -o {output} \
        {input.aligned_bam}
        """

if end_type == "pe":
    rule bam_to_fastq:
        input:
            rules.name_sort.output

        output:
            one = os.path.join(out_fastq_dir, "{sample}.1.pulled.fastq.gz"),
            two = os.path.join(out_fastq_dir, "{sample}.2.pulled.fastq.gz"),

        params:
            singletons = os.path.join(out_fastq_dir, "{sample}.singletons.fastq.gz")

        conda: "../envs/single_steps.yaml"

        threads:
            config["fastq_extra_threads"] + 1

        shell:
            """
            samtools fastq \
            -1 {output.one} \
            -2 {output.two} \
            -s {params.singletons} \
            -N \
            -@ {threads} \
            {input}
            """

else:
    rule bam_to_fastq:
        input:
            rules.name_sort.output

        output:
            one = os.path.join(out_fastq_dir, "{sample}.1.pulled.fastq.gz")

        threads:
            config["fastq_extra_threads"] + 1

        conda: "../envs/single_steps.yaml"

        shell:
            """
            samtools fastq \
            -@ {threads} \
            -0 {output.one} \
            {input}
            """


rule copy_config:
    input:
        conf = config_path,
        fqs = expand(os.path.join(out_fastq_dir, "{sample}.1.pulled.fastq.gz"),
                     sample = SAMPLES)
    output:
        os.path.join(out_fastq_dir, "config_sort_pull.yaml")

    shell:
        """
        cp {input.conf} {output}
        """
