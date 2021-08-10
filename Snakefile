import pandas as pd
import glob

configfile: "config_grid.yaml"

def find_path(dir, barcode):
    # path = glob.glob("{}/**/fastq_pass/barcode{}/".format(dir,barcode), recursive=True)
    path = glob.glob("{}/{}".format(dir,barcode))
    return str(path[0])

sample_table = pd.read_csv("ref/sample_table.csv").set_index("Sample")
sample_table["path"] = sample_table.apply(lambda x: find_path(x.Directory, x.Barcode), axis = 1)

SAMPLES = sample_table.index.values

TIME = [10, 15]

print(SAMPLES, TIME)

include: "rules/move_files.smk"
include: "rules/host_remove.smk"

rule all:
    input:
        expand("results/{sample}/{time}_minutes/microbial/hg38_unmapped.fastq", sample = SAMPLES, time = TIME)

# rule test:
#     input:
#         lambda wildcards: sample_table.path[wildcards.sample]
#     output:
#         "results/{sample}_merged.fastq"
#     shell:
#         "cat {input}/*.fastq > {output}"


