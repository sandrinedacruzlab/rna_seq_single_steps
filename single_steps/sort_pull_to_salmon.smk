

#########---------------
## Input parameters
#########---------------
configfile: "config/sort_pull_to_salmon.yaml"
import os

project_dir = config["project_dir"]

salmon_index_dir = config["salmon_index_dir"]
gtf = config["gtf"]
salmon_strand_info = config["salmon_strand_info"]
end_type = config["end_type"]


in_bam_dir = config["input_dir"]
out_temp_dir = config["temp_dir"]
out_fastq_dir = config["fastq_spot"]
salmon_out_dir = config['salmon_out']


bam_suffix = config["bam_suffix"]


assert end_type in ["pe", "se"], f"'end_type' must be one of 'pe' (paired-end) or 'se' (single-end) - {end_type} was provided"
print(end_type)
SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]
print(SAMPLES)
# sys.stderr.write("".join(str(s) for s in SAMPLES) + "\n")
if not os.path.exists(out_fastq_dir):
    os.system(f"mkdir -p {out_fastq_dir}")

if not os.path.exists(salmon_out_dir):
    os.system(f"mkdir -p {salmon_out_dir}")

if not os.path.exists(out_temp_dir):
    os.system(f"mkdir -p {out_temp_dir}")

# conda: "../envs/single_steps.yaml"



rule all:
  input:
    expand(os.path.join(salmon_out_dir, "{sample}/quant.sf"), sample = SAMPLES)

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
        priority: 1
        input:
            rules.name_sort.output

        output:
            one = temp(os.path.join(out_fastq_dir, "{sample}.1.pulled.fastq.gz")),
            two = temp(os.path.join(out_fastq_dir, "{sample}.2.pulled.fastq.gz")),

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
        priority: 1
        input:
            rules.name_sort.output

        output:
            one = temp(os.path.join(fastq_spot, "{sample}.1.pulled.fastq.gz"))

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

if end_type == "pe":
    rule salmon_quant:
        priority: 2
        input:
            fast1 = os.path.join(out_fastq_dir, "{sample}.1.pulled.fastq.gz"),
            fast2 = os.path.join(out_fastq_dir, "{sample}.2.pulled.fastq.gz")

        output:
            os.path.join(salmon_out_dir, "{sample}/quant.sf")


        params:
            index_dir = salmon_index_dir,
            salmon_out_dir = os.path.join(salmon_out_dir, "{sample}"),
            libtype = config["salmon_strand_info"],
            extra_flags = " ".join(config["salmon_quant_flags"])

        threads:
            config["quant_threads"]

        conda: "../envs/single_steps.yaml"

        shell:
            """
            salmon quant \
            --index {params.index_dir} \
            --libType {params.libtype} \
            --mates1 {input.fast1} \
            --mates2 {input.fast2} \
            --threads {threads} \
            -o {params.salmon_out_dir} \
            {params.extra_flags} 
            """

else:
    rule salmon_quant:
        priority: 2
        input:
            fast1 = os.path.join(fastq_spot, "{sample}.1.pulled.fastq.gz")

        output:
            os.path.join(salmon_out_dir, "{sample}/quant.sf")

        params:
            index_dir = salmon_index_dir,
            salmon_out_dir = os.path.join(salmon_out_dir, "{sample}"),
            libtype = config["salmon_strand_info"],
            extra_flags = " ".join(config["salmon_quant_flags"])

        threads:
            config["quant_threads"]

        conda: "../envs/single_steps.yaml"

        shell:
            """
            salmon quant \
            --index {params.index_dir} \
            --libType {params.libtype} \
            -r {input.fast1} \
            --threads {threads} \
            -o {params.salmon_out_dir} \
            {params.extra_flags} 
            """
