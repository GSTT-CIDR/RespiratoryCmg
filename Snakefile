import pandas as pd
import glob
import time

configfile: "config.yaml"

def find_path(dir, barcode):
    # Work on new method for this bit    
    path = glob.glob("/data/{}/**/fastq_pass/barcode*{}/".format(dir,barcode), recursive=True)
    #path = glob.glob("{}/{}".format(dir,barcode))
    return str(path[0])

print("waiting 1 minute before running")
time.sleep(60)

sample_table = pd.read_csv("ref/sample_table.tsv", sep="\t").set_index("Sample")
sample_table["path"] = sample_table.apply(lambda x: find_path(x.Directory, x.Barcode), axis = 1)

SAMPLES = sample_table.index.values

TIME = [120]

print(SAMPLES, TIME)

include: "rules/move_files.smk"
include: "rules/host_remove.smk"
include: "rules/centrifuge.smk"
include: "rules/amr.smk"
include: "rules/qc.smk"


rule all:
    input:
        expand("results/{sample}/{time}_minutes/amr/scagaire_gene_summary.tsv", sample = SAMPLES, time = TIME),
        expand("results/{sample}/{time}_minutes/qc/nanostat_summary.txt", sample = SAMPLES, time = TIME)    



