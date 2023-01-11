configfile: "config/feature_counts_config.yaml"
import os



config_path = "config/feature_counts_config.yaml" # this is for copying later


# project_dir = "/SAN/vyplab/alb_projects/data/linked_buratti_hnrnpk/"
out_spot = config["out_spot"]
bam_spot = config["bam_spot"]

# mouse and human gtf, comment dependening on your species
# gtf =  "/SAN/vyplab/vyplab_reference_genomes/annotation/mouse/gencode/gencode.vM22.annotation.gtf"
# gtf =  "/SAN/vyplab/vyplab_reference_genomes/annotation/human/GRCh38/gencode.v31.annotation.gtf"
gtf = config["gtf"]

feature_counts_strand_info = config["feature_counts_strand_info"]
end_type = config["end_type"]
bam_suffix = config["bam_suffix"]

# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------

output_dir = os.path.join(out_spot, "")
bam_dir = os.path.join(bam_spot, "")

SAMPLES = [f.replace(bam_suffix, "") for f in os.listdir(bam_dir) if f.endswith(bam_suffix)]

print(SAMPLES)

if not os.path.exists(output_dir):
    os.system("mkdir -p {0}".format(output_dir))


localrules: all, copy_config

rule all:
  input:
    expand(output_dir + "{sample}_featureCounts_results.txt", sample = SAMPLES),
    os.path.join(output_dir, "feature_counts_config.yaml")
if config["gtf_used"]:
    rule feature_counts:
        input:
            aligned_bam = os.path.join(bam_dir,"{sample}" + bam_suffix)

        output:
            out_name = os.path.join(output_dir, "{sample}_featureCounts_results.txt")

        params:
            ref_anno = gtf,
            stranded = feature_counts_strand_info,
            feature_type = config["feature_type"], # exons extracted for counting from ref_anno
            attr_type = config["attribute_type"], # Generate meta-features for counting via this group/id
            extra_attr = ",".join(["gene_name"]), # extra metadata to extract & report in count output
            count_level = "-f" if config["count_at"] == "feature" else "", # count at exon or gene level?
            paired_end = "-p" if end_type == "pe" else "", #count fragments if paired end
            long_reads = "-L" if config["long_reads"] else "" # turn on/off long reads setting

        container:
            "docker://quay.io/biocontainers/subread:2.0.3--h7132678_0"

        shell:
            """
            featureCounts \
            -a {params.ref_anno} \
            -t {params.feature_type} \
            -g {params.attr_type} \
            --extraAttributes {params.extra_attr} \
            {params.count_level} \
            {params.paired_end} \
            {params.stranded} \
            {params.long_reads} \
            -o {output.out_name} \
            {input.aligned_bam}
            """
else:
    rule feature_counts:
        input:
            aligned_bam = os.path.join(bam_dir,"{sample}" + bam_suffix)

        output:
            out_name = os.path.join(output_dir, "{sample}_featureCounts_results.txt")

        params:
            ref_anno = config["saf"],
            stranded = feature_counts_strand_info,
            paired_end = "-p" if end_type == "pe" else "", #count fragments if paired end
            long_reads = "-L" if config["long_reads"] else "" # turn on/off long reads setting

        container:
            "docker://quay.io/biocontainers/subread:2.0.3--h7132678_0"

        shell:
            """
            featureCounts \
            -a {params.ref_anno} \
            -F SAF \
            {params.paired_end} \
            {params.stranded} \
            {params.long_reads} \
            -o {output.out_name} \
            {input.aligned_bam}
            """

rule copy_config:
    input:
        fc_output = expand(output_dir + "{sample}_featureCounts_results.txt", sample = SAMPLES),
        conf = config_path

    output:
        os.path.join(output_dir, "feature_counts_config.yaml")

    shell:
        """
        cp {input.conf} {output}
        """
