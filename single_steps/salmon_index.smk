configfile: "config/salmon_config.yaml"

#########---------------
## Input parameters
#########---------------
config_path = "config/salmon_config.yaml" # this is for copying later

transcripts = config["transcripts"]
threads = config["threads"]
output_index_dir = config["salmon_index_dir"]

rule salmon_index:
    input:
        transcripts = transcripts
    output:
        directory(output_index_dir)
    params:
        threads = threads
    shell:
        """
        /SAN/vyplab/alb_projects/tools/salmon-latest_linux_x86_64/bin/salmon index \
        -t {input.transcripts} -i {output} --threads {threads} --gencode
        """

