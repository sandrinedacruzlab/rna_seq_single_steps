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

rule all:
  input:
    expand(os.path.join(out_fastq_dir, "{sample}.1.pulled.fastq.gz"),
           sample = SAMPLES),
    expand(os.path.join(out_fastq_dir, "{sample}.2.pulled.fastq.gz"),
           sample = SAMPLES if end_type == "pe" else []),

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
        1
        #config["collate_extra_threads"] + 1

    conda:
        "../envs/single_steps.yaml"

    container:
        "docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

    shell:
        """
        EXTRA_THREADS=$(({threads} -1))
        samtools collate \
        -@ $EXTRA_THREADS \
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
            singletons = os.path.join(out_fastq_dir, "{sample}.singletons.fastq.gz")

#        params:
#             singletons=subpath(output.one, strip_suffix=".1.pulled.fastq.gz") + ".singletons.fastq.gz"
#            singletons = os.path.join(out_fastq_dir, "{sample}.singletons.fastq.gz")

        conda: "../envs/single_steps.yaml"

        container:
            "docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"

        threads:
             1
#            config["fastq_extra_threads"] + 1

        shell:
            """
            EXTRA_THREADS=$(({threads} -1))
            samtools fastq \
            -1 {output.one} \
            -2 {output.two} \
            -s {output.singletons} \
            -N \
            -@ $EXTRA_THREADS \
            {input}
            """

else:
    rule bam_to_fastq:
        input:
            rules.name_sort.output

        output:
            one = os.path.join(out_fastq_dir, "{sample}.1.pulled.fastq.gz")

        threads:
            1
            #config["fastq_extra_threads"] + 1

        conda: "../envs/single_steps.yaml"
        container:
            "docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0"
        shell:
            """
            EXTRA_THREADS=$(({threads} -1))
            samtools fastq \
            -@ $EXTRA_THREADS \
            -0 {output.one} \
            {input}
            """

