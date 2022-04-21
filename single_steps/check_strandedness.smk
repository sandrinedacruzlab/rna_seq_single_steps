configfile: "config/check_strandedness_config.yaml"
include: "/SAN/vyplab/alb_projects/pipelines/rna_seq_snakemake/rules/helpers.py"
import os

#########---------------
## Input parameters
#########---------------
config_path = "config/check_strandedness_config.yaml" # this is for copying later
csv_path = config["sampleCSVpath"]

gtf = config["gtf"]
transcripts = config["transcripts"]

output_dir = config["output_dir"]

SAMPLES = pd.read_csv(csv_path, sep = ",")
SAMPLES = SAMPLES.replace(np.nan, '', regex=True)

UNITS = SAMPLES['unit'].tolist()
SAMPLE_NAMES = SAMPLES['sample_name'].tolist()

'''
#so I want this rule to be run ONCE for every fast1, so the wild cards I'm giving are the 'name' of the fastq of the first read
#name = [re.sub(".fastq.gz","",strpd.rpartition('/')[2]) for strpd in SAMPLES['fast1'].tolist()]
name1 = SAMPLES['fast1'][0]
if config['end_type'] == 'pe':
    name1 = SAMPLES['fast2'][0]
'''

if not os._exists(output_dir):
    os.system("mkdir -p {0}".format(output_dir))

########-----------------
localrules: all, copy_config

rule all:
    input:
        expand(os.path.join(output_dir, "{unit}_{name}.cntTable"),zip, unit = UNITS, name = SAMPLE_NAMES),
        os.path.join(output_dir, "config_strandedness.yaml")

rule check_strandedness:
    input:
        os.path.join(output_dir, "config_strandedness.yaml"),
        fastq_file = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = True),
        fastq_file2 = lambda wildcards: return_fastq(wildcards.name,wildcards.unit,first_pair = False),
    output:
        os.path.join(output_dir, "{unit}_{name}.cntTable")
    params:
        transcripts = transcripts,
        gtf = gtf
    shell:
        """
        set +u;
        source activate kallisto
        check_strandedness  --gtf {params.gtf} \
        --transcripts {params.transcripts}  \
        --reads_1 {input.fastq_file} \
        --reads_2 {input.fastq_file2} \
        -n 10000
        """

rule copy_config:
    input:
        conf = config_path
    output:
        os.path.join(output_dir, "config_strandedness.yaml")
    shell:
        """
        cp {input.conf} {output}
        """