import pandas as pd
import os
from jinja2 import Environment, FileSystemLoader
import re
from weasyprint import HTML, CSS
from datetime import datetime

#RESULTS_DIR = "../../validation/results/160921_validation_run4_sample15/120_minutes/"
#SAMPLE = "John Hancock"
#INTERVAL = 2
#CFG_PATH = RESULTS_DIR + "centrifuge/centrifuge_report.tsv"
#AMR_SUMMARY = RESULTS_DIR + "amr/scagaire_gene_summary.tsv"
#AMR_REPORT = RESULTS_DIR + "amr/scagaire_report.tsv"
#QC_PATH = RESULTS_DIR + "qc/nanostat_summary.txt"
#OUTPUT = "../Tester.pdf"
#REPORT_HTML = "ref/Template/report_template.html"
#REPORT_CSS = "ref/Template/report.css"
#BOOTSTRAP_CSS = "ref/Template/bootstrap.css"

SAMPLE = snakemake.wildcards.sample
INTERVAL = 2
CFG_PATH = snakemake.input.centrifuge
AMR_SUMMARY = snakemake.input.amr_summary
AMR_REPORT = snakemake.input.amr_report
QC_PATH = snakemake.input.qc
OUTPUT = str(snakemake.output)
REPORT_HTML = snakemake.config["pdf"]["html"]
REPORT_CSS = snakemake.config["pdf"]["css"]
BOOTSTRAP_CSS = snakemake.config["pdf"]["bootstrap"]


def convert_bp(size):
    size = float(str(size).replace(",", ""))
    for x in ["bp", "Kb", "Mb", "Gb", "Tb"]:
        if size < 1000.0:
            return "{:.2f} {}".format(size, x)
        size /= 1000.0
    return size


def summary_qc(path):
    qc_dict = dict()
    with open(path, "r") as f:
        data = f.read().splitlines()
    for line in data:
        if line.startswith("Mean read length"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["mean_read_len"] = dat
        elif line.startswith("Mean read quality"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["mean_read_qual"] = dat
        elif line.startswith("Number of reads"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["num_reads"] = dat
        elif line.startswith("Total bases"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["total_bp"] = convert_bp(dat)
    return qc_dict


def cfg_to_html(path, threshold = 1):
    cfg_dict = dict()
    df = pd.read_csv(path, sep="\t")
    df = df.round(2)
    cfg_dict["cfg_full"] = df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    ic = df[df["Organism"] == "Jonesia denitrificans"]
    if ic.empty:
        cfg_dict["ic"] = "NA/NA"
    else:
        cfg_dict["ic"] = "{}/{}%".format(int(ic["Counts"]), float(ic["Percentage"]))
    above_df = df.copy()
    above_df = above_df.drop(columns=["Tax_ID"])
    above_df = above_df[above_df["Percentage"] > threshold]
    cfg_dict["cfg_top"] = above_df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    return cfg_dict


def amr_summary(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["Species", "Gene", "Hits"])
    res = df[df["Hits"] > 1]
    return res.to_html(classes="table table-striped", border=0, justify="left", index=False)


def amr_report(path):
    with open(path, "r") as f:
        dat = f.readline().strip()
        if dat.startswith("#FILE"):
            df = pd.read_csv(path, sep="\t", usecols=["SEQUENCE", "START", "END", "GENE", "%COVERAGE"])
            return df.to_html(classes="table table-striped", border=0, justify="left", index=False)
        else:
            return dat

print(os.getcwd())
report_dict = {"name": SAMPLE,
               "time": str(INTERVAL) + " hrs",
               "title": "Clinical metagenomics report",
               "date": datetime.now()}
report_dict.update(cfg_to_html(CFG_PATH))
report_dict.update(summary_qc(QC_PATH))
report_dict["amr_summary"] = amr_summary(AMR_SUMMARY)
report_dict["amr_report"] = amr_report(AMR_REPORT)


env = Environment(loader=FileSystemLoader(".")) # Change to "." for grid
# declare our jinja template
template = env.get_template(REPORT_HTML)

html_out = template.render({"report": report_dict})

pdf_name = OUTPUT
HTML(string=html_out).write_pdf(pdf_name,
                                stylesheets=[CSS(REPORT_CSS),
                                             CSS(BOOTSTRAP_CSS)])
