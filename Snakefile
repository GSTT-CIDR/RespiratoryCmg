import pandas as pd
import glob
import time
import os
import sys
configfile: "config.yaml"


def find_path(exp, sample, barcode):
    barcode_dict = {"1": "01", "2":"02", "3":"03", "4":"04", "5":"05", "6":"06", "7":"07", "8":"08", "9":"09", "10":"10", "11":"11", "12":"12a"}
    path_str = "/data/{}/{}/**/fastq_pass/barcode{}/".format(exp, sample, barcode_dict[str(barcode)])
    path = glob.glob(path_str, recursive=True)
    if len(path) == 0:
        sys.exit("Error: path error with sample {}. Check variables for experiment, sample or barcode".format(sample))
    else:
        return path[0]   
	
#print("Waiting 1 minute before running")
#time.sleep(60)

sample_table = pd.read_csv(config["samples"], sep="\t").set_index("LabID")
sample_table["path"] = sample_table.apply(lambda x: find_path(x.Experiment, x.SampleID, x.Barcode), axis = 1)

SAMPLES = sample_table.index.values

TIME = config["time"]# move to config file

#print(SAMPLES, TIME)

include: "rules/move_files.smk"
include: "rules/host_remove.smk"
include: "rules/centrifuge.smk"
include: "rules/amr.smk"
include: "rules/qc.smk"
include: "rules/report.smk"
include: "rules/mlst.smk"


rule all:
    input:
        expand("reports/{sample}/{sample}_{time}_hours_report.pdf", sample = SAMPLES, time = TIME),
	expand("results/{sample}/{time}_hours/transfer/transferred.txt", sample = SAMPLES, time = TIME)
        #expand("reports/{sample}/{time}_hours/mlst/", sample= SAMPLES, time= TIME)
