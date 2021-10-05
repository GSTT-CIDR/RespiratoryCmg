import pandas as pd
import glob
import time

configfile: "config.yaml"

def find_path(exp, sample,git  barcode):
    if int(barcode) < 10:
        barcode = "0" + str(barcode)
    # Work on new method for this bit    
    path = glob.glob("/data/{}/{}/**/fastq_pass/barcode*{}/".format(exp, sample, barcode), recursive=True)
    #path = glob.glob("{}/{}".format(dir,barcode))
    return str(path[0])

#print("waiting 2 minute before running")
#time.sleep(120)

sample_table = pd.read_csv(config["samples"], sep="\t").set_index("Sample")
#sample_table["path"] = sample_table.apply(lambda x: find_path(x.Directory, x.Barcode), axis = 1)

SAMPLES = sample_table.index.values

TIME = [120] # move to config file

print(SAMPLES, TIME)

#include: "rules/move_files.smk"
include: "rules/host_remove.smk"
include: "rules/centrifuge.smk"
include: "rules/amr.smk"
include: "rules/qc.smk"
include: "rules/report.smk"


rule all:
    input:
        expand("reports/{sample}/{sample}_{time}_minutes_report.pdf", sample = SAMPLES, time = TIME)
