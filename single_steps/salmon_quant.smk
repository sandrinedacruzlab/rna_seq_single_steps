configfile: "config/salmon_config.yaml"
import os

#########---------------
## Input parameters
#########---------------
config_path = "config/salmon_config.yaml" # this is for copying later

project_dir = config["project_dir"]
out_spot = config["out_spot"]
fastq_spot = config["fastq_spot"]
salmon_index_dir = config["salmon_index_dir"]

gtf = config["gtf"]
salmon_strand_info = config["salmon_strand_info"]
end_type = config["end_type"]

# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------

output_dir = os.path.join(project_dir,out_spot)
fastq_dir = os.path.join(project_dir,fastq_spot)

SAMPLES, = glob_wildcards(os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz"))
print(SAMPLES)

if not os._exists(output_dir):
    os.system("mkdir -p {0}".format(output_dir))

rule all:
  input:
    expand(os.path.join(output_dir, "{sample}_quant.sf"), sample = SAMPLES),
    os.path.join(project_dir, "config_salmon_quant.yaml")

if end_type == "pe":
    rule salmon_quant:
        input:
            fast1 = os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz"),
            fast2 = os.path.join(fastq_dir, "{sample}_2.merged.fastq.gz"),
        output:
            os.path.join(output_dir, "{sample}_quant.sf"),
        params:
            tempout = os.path.join(output_dir, "{sample}"),
            tempout2 = os.path.join(output_dir, "{sample}", "quant.sf"),
            salmon_strand_info = salmon_strand_info,
            salmon_index_dir = salmon_index_dir,
            gtf = gtf
        shell:
            """
            /SAN/vyplab/alb_projects/tools/salmon-latest_linux_x86_64/bin/salmon quant\
             -i {params.salmon_index_dir} \
             -l {params.salmon_strand_info} \
             -1 {input.fast1} \
             -2 {input.fast2} \
             -o {params.tempout}\
             --geneMap {params.gtf}
             mv {params.tempout2} {output}
            """
else:
    rule salmon_quant:
        input:
            fast1 = os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz")
        output:
            os.path.join(output_dir, "{sample}_quant.sf"),
        params:
            tempout = os.path.join(output_dir, "{sample}"),
            tempout2 = os.path.join(output_dir, "{sample}", "quant.sf"),
            salmon_strand_info = salmon_strand_info,
            salmon_index_dir = salmon_index_dir,
            gtf = gtf
        shell:
            """
            /SAN/vyplab/alb_projects/tools/salmon-latest_linux_x86_64/bin/salmon quant\
             -i {params.salmon_index_dir} \
             -l {params.salmon_strand_info} \
             -r {input.fast1} \
             -o {params.tempout}\
             --geneMap {params.gtf}
             mv {params.tempout2} {output}
            """

rule copy_config:
    input:
        conf = config_path
    output:
        os.path.join(project_dir, "config_salmon_quant.yaml")
    shell:
        """
        cp {input.conf} {output}
        """
