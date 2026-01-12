configfile: "config/htseq_count_config.yaml"
import os
import sys


#########---------------
## Input parameters
#########---------------

in_bam_dir = config["input_dir"]
out_counts_dir = config["output_dir"]
gtf_file = config["gtf_file"]

bam_suffix = config["bam_suffix"]

# HTSeq parameters from config (i.e. a nested dictionary)
htseq_params = config["htseq_params"]

# Get sample names
SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]

sys.stderr.write(f"Basenames for input BAM files - {', '.join(SAMPLES)}\n")

# Create output directories if needed
if not os.path.exists(out_counts_dir):
    os.makedirs(out_counts_dir)

out_aln_dir = os.path.join(out_counts_dir, "samout")
if htseq_params["samout_enabled"]:
    if not os.path.exists(out_aln_dir):
        os.makedirs(out_aln_dir)


# def get_samout_param(samout_enabled, samout_dir, sample_names, samout_format):
#     """Generate the -o parameter string for htseq-count samout option"""
#     if not samout_enabled:
#         return ""
    
#     out_sam_files = [os.path.join(samout_dir, f"{sample}.annotated.{samout_format}") 
#                      for sample in sample_names]
#     return "-o " + " ".join(out_sam_files)


# Determine output files for rule all
# always the count matrix, but include individual BAM files if samout option is enabled
def get_all_outputs():
    outputs = [os.path.join(out_counts_dir, "counts.tsv")]
    
    # Add SAM output files if enabled
    if htseq_params.get("samout_enabled", False):
        samout_dir = os.path.join(out_counts_dir, "samout")
        outputs.extend(
            expand(os.path.join(samout_dir, "{sample}.annotated." + htseq_params["samout_format"]),
                   sample=SAMPLES)
        )

    return outputs


if htseq_params["samout_enabled"]:
    ruleorder: htseq_count_output_alignments > htseq_count
else:
    ruleorder: htseq_count > htseq_count_output_alignments

wildcard_constraints:
    sample="|".join(SAMPLES)


rule all:
    input:
        get_all_outputs()

rule htseq_count:
    input:
        bam_files = expand(os.path.join(in_bam_dir, "{sample}" + bam_suffix), sample=SAMPLES),
        gtf = gtf_file
    
    output:
        counts = os.path.join(out_counts_dir, "counts.tsv")
    
    params:
        order = htseq_params["order"],
        max_reads_in_buffer = htseq_params["max_reads_in_buffer"],
        stranded = htseq_params["stranded"],
        minaqual = htseq_params["minaqual"],
        htseq_type = htseq_params["type"],
        idattr = htseq_params["idattr"],
        additional_attr = " ".join([f"--additional-attr={attr}" for attr in htseq_params["additional_attr"]]) if htseq_params["additional_attr"] else "",
        add_chromosome_info = "--add-chromosome-info" if htseq_params["add_chromosome_info"] else "",
        feature_query = f"--feature-query='{htseq_params['feature_query']}'" if htseq_params["feature_query"] else "",
        mode = htseq_params["mode"],
        nonunique = htseq_params["nonunique"],
        secondary_alignments = htseq_params["secondary_alignments"],
        supplementary_alignments = htseq_params["supplementary_alignments"],
        nprocesses = htseq_params["nprocesses"],
        quiet = "--quiet" if htseq_params["quiet"] else ""

    log:
        stdout=os.path.join(out_counts_dir, "logs", "htseq_count.stdout.txt"),
        stderr=os.path.join(out_counts_dir, "logs", "htseq_count.stderr.txt")
    
    threads:
        1
    
    container:
        "docker://quay.io/biocontainers/htseq:2.0.9--py311h8fb3dee_0"
    
    shell:
        """
        htseq-count \
        --order={params.order} \
        --max-reads-in-buffer={params.max_reads_in_buffer} \
        --stranded={params.stranded} \
        --minaqual={params.minaqual} \
        --type={params.htseq_type} \
        --idattr={params.idattr} \
        {params.additional_attr} \
        {params.add_chromosome_info} \
        {params.feature_query} \
        --mode={params.mode} \
        --nonunique={params.nonunique} \
        --secondary-alignments={params.secondary_alignments} \
        --supplementary-alignments={params.supplementary_alignments} \
        --nprocesses={threads} \
        --counts_output={output.counts} \
        {params.quiet} \
        {input.bam_files} \
        {input.gtf} \
        1> {log.stdout} \
        2> {log.stderr}
        """

rule htseq_count_output_alignments:
    input:
        bam_files = expand(os.path.join(in_bam_dir, "{sample}" + bam_suffix), sample=SAMPLES),
        gtf = gtf_file
    
    output:
        counts = os.path.join(out_counts_dir, "counts.tsv"),
        samout = expand(os.path.join(out_aln_dir, "{sample}.annotated." + htseq_params["samout_format"]), 
                       sample=SAMPLES) if htseq_params.get("samout_enabled", False) else []
    
    params:
        order = htseq_params["order"],
        max_reads_in_buffer = htseq_params["max_reads_in_buffer"],
        stranded = htseq_params["stranded"],
        minaqual = htseq_params["minaqual"],
        htseq_type = htseq_params["type"],
        idattr = htseq_params["idattr"],
        additional_attr = " ".join([f"--additional-attr={attr}" for attr in htseq_params["additional_attr"]]) if htseq_params["additional_attr"] else "",
        add_chromosome_info = "--add-chromosome-info" if htseq_params["add_chromosome_info"] else "",
        feature_query = f"--feature-query='{htseq_params['feature_query']}'" if htseq_params["feature_query"] else "",
        mode = htseq_params["mode"],
        nonunique = htseq_params["nonunique"],
        secondary_alignments = htseq_params["secondary_alignments"],
        supplementary_alignments = htseq_params["supplementary_alignments"],
        samout = output.samout,
        # samout = "-o " + " ".join(list(output.samout)) if htseq_params["samout_enabled"] else "",
        samout_format = htseq_params["samout_format"],
        nprocesses = htseq_params["nprocesses"],
        quiet = "--quiet" if htseq_params["quiet"] else ""

    log:
        stdout=os.path.join(out_counts_dir, "logs", "htseq_count_output_alignments.stdout.txt"),
        stderr=os.path.join(out_counts_dir, "logs", "htseq_count_output_alignments.stderr.txt")
    
    threads:
        1
    
    container:
        "docker://quay.io/biocontainers/htseq:2.0.9--py311h8fb3dee_0"
    
    shell:
        """
        htseq-count \
        --order={params.order} \
        --max-reads-in-buffer={params.max_reads_in_buffer} \
        --stranded={params.stranded} \
        --minaqual={params.minaqual} \
        --type={params.htseq_type} \
        --idattr={params.idattr} \
        {params.additional_attr} \
        {params.add_chromosome_info} \
        {params.feature_query} \
        --mode={params.mode} \
        --nonunique={params.nonunique} \
        --secondary-alignments={params.secondary_alignments} \
        --supplementary-alignments={params.supplementary_alignments} \
        -o {params.samout} \
        --samout-format={params.samout_format} \
        --nprocesses={threads} \
        --counts_output={output.counts} \
        {params.quiet} \
        {input.bam_files} \
        {input.gtf} \
        1> {log.stdout} \
        2> {log.stderr}
        """

# def get_samout_param(samout_enabled):
#     """Generate the -o parameter string for htseq-count samout option"""
#     if not samout_enabled:
#         return ""
    
#     out_sam_files = rules.htseq_count.output.samout
#     return "-o " + " ".join(out_sam_files)