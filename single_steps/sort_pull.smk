configfile: "config/sort_pull_config.yaml"
import os

# =-------DON"T TOUCH ANYTHING PAST THIS POINT ----------------------------

#########---------------
## Input parameters
#########---------------
config_path = "config/sort_pull_config.yaml" # this is for copying later

project_dir = config["project_dir"]
out_spot = config["out_spot"]
bam_spot = config["bam_spot"]
fastq_spot = config["fastq_spot"]

output_dir = os.path.join(project_dir,out_spot)
bam_dir = os.path.join(project_dir,bam_spot)
fastq_dir = os.path.join(project_dir, fastq_spot)

bam_suffix = config["bam_suffix"]
end_type = config["end_type"]
SAMPLES, = glob_wildcards(os.path.join(bam_dir, "{sample}" + bam_suffix))
print(SAMPLES)

if not os._exists(fastq_dir):
    os.system("mkdir -p {0}".format(fastq_dir))

if not os._exists(output_dir):
    os.system("mkdir -p {0}".format(output_dir))

rule all:
  input:
    expand(os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz"), sample = SAMPLES),
    os.path.join(output_dir, "config_sort_pull.yaml")

rule name_sort:
    input:
        aligned_bam = os.path.join(bam_dir, "{sample}" + bam_suffix)
    output:
       temp(os.path.join(output_dir, "{sample}_namesorted.bam"))
    shell:
        """
        samtools sort -n -@ 2 {input.aligned_bam} -o {output}
        """
rule copy_config:
    input:
        conf = config_path
    output:
        os.path.join(output_dir, "config_sort_pull.yaml")
    shell:
        """
        cp {input.conf} {output}
        """
if end_type == "pe":
  rule bam_to_fastq:
      input:
          name_sort_bam = os.path.join(output_dir, "{sample}_namesorted.bam")
      output:
          one = temp(os.path.join(fastq_dir, "{sample}_1.merged.fastq")),
          two = temp(os.path.join(fastq_dir, "{sample}_2.merged.fastq"))
      shell:
          """
          bedtools bamtofastq -i {input} \
                        -fq {output.one} \
                        -fq2 {output.two}
          """
  rule gunzip_fastq:
      input:
          one = os.path.join(fastq_dir, "{sample}_1.merged.fastq"),
          two = os.path.join(fastq_dir, "{sample}_2.merged.fastq")
      output:
          one_out = os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz"),
          two_out = os.path.join(fastq_dir, "{sample}_2.merged.fastq.gz")
      shell:
          """
          gzip {input.one}
          gzip {input.two}
          """
else:
  rule bam_to_fastq:
      input:
          name_sort_bam = os.path.join(output_dir, "{sample}_namesorted.bam")
      output:
          one = temp(os.path.join(fastq_dir, "{sample}_1.merged.fastq"))
      shell:
          """
          bedtools bamtofastq -i {input} \
                        -fq {output.one}
          """
  rule gunzip_fastq:
      input:
          one = os.path.join(fastq_dir, "{sample}_1.merged.fastq")
      output:
          one_out = os.path.join(fastq_dir, "{sample}_1.merged.fastq.gz")
      shell:
          """
          gzip {input.one}
          """
