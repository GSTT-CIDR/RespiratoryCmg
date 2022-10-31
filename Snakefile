import pandas as pd
import glob
import time
import os
import sys
configfile: "config.yaml"


def find_path(exp, sample, barcode, data_dir):
    barcode_dict = {"1": "01", "2":"02", "3":"03", "4":"04", "5":"05", "6":"06", "7":"07", "8":"08", "9":"09", "10":"10", "11":"11", "12":"12a"}
    path_str = f"{str(data_dir)}/{exp}/{sample}/**/fastq_pass/barcode{barcode_dict[str(barcode)]}/"
    path = glob.glob(path_str, recursive=True)
    if len(path) == 0:
        sys.exit(f"Error: path error with sample {exp} barcode {barcode}. Check variables for experiment, sample or barcode")
    else:
        return path[0]   
	
#print("Waiting 1 minute before running")
#time.sleep(60)
data_dir = config["data_dir"]
sample_table = pd.read_csv(config["samples"], sep="\t").set_index("LabID")
sample_table["path"] = sample_table.apply(lambda x: find_path(x.Experiment, x.SampleID, x.Barcode, data_dir), axis = 1)

SAMPLES = sample_table.index.values

TIME = config["time"]# move to config file

print(SAMPLES, TIME)

include: "rules/move_files.smk"
include: "rules/host_remove.smk"
include: "rules/centrifuge.smk"
include: "rules/amr.smk"
include: "rules/qc.smk"
include: "rules/report.smk"
include: "rules/mlst.smk"
include: "rules/viral.smk"


rule all:
    input:
        expand("reports/{sample}/{sample}_{time}_hours_report.pdf", sample = SAMPLES, time = TIME),
        # expand("results/{sample}/{time}_hours/transfer/transferred.txt", sample = SAMPLES, time = TIME),
        expand("results/{sample}/{time}_hours/viral/centrifuge_viral_report.tsv", sample = SAMPLES, time = TIME),
        expand("results/{sample}/{time}_hours/mlst/", sample= SAMPLES, time= TIME),
        expand("results/{sample}/{time}_hours/amr/virulence_factor_summary.tsv", sample=SAMPLES, time=TIME)
