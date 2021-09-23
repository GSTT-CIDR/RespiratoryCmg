import pandas as pd
import numpy as np
from jinja2 import Environment, FileSystemLoader
import re
from weasyprint import HTML, CSS
from datetime import datetime

RESULTS_DIR = "../../Validation/August_september_2021/results/160921_validation_run4_sample15/120_minutes/"
SAMPLE = "John Hancock"
INTERVAL = 2
CFG_PATH = RESULTS_DIR + "centrifuge/centrifuge_report.tsv"
AMR_SUMMARY = RESULTS_DIR + "amr/scagaire_gene_summary.tsv"
AMR_REPORT = RESULTS_DIR + "amr/scagaire_report.tsv"
QC_PATH = RESULTS_DIR + "qc/nanostat_summary.txt"


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
    cfg_dict["cfg_full"] = df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    above_df = df.copy()
    above_df = above_df.rename(columns={"Tax_ID" : "TaxID"})
    above_df = above_df.round(2)
    above_df = above_df[above_df["Percentage"] > threshold]
    cfg_dict["cfg_top"] = above_df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    return cfg_dict


def amr_summary(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["Species", "Gene", "Hits"])
    res = df[df["Hits"] > 1]
    return res.to_html(classes="table table-striped", border=0, justify="left", index=False)


def amr_report(path):
    df = pd.read_csv(path, sep="\t", usecols=["SEQUENCE", "START", "END", "GENE", "%COVERAGE"])
    return df.to_html(classes="table table-striped", border=0, justify="left", index=False)


report_dict = {"name": SAMPLE,
               "time": str(INTERVAL) + " hrs",
               "title": "Clinical metagenomics report",
               "date": datetime.now()}
report_dict.update(cfg_to_html(CFG_PATH))
report_dict.update(summary_qc(QC_PATH))
report_dict["amr_summary"] = amr_summary(AMR_SUMMARY)
report_dict["amr_report"] = amr_report(AMR_REPORT)


env = Environment(loader=FileSystemLoader("../"))
# declare our jinja template
template = env.get_template("ref/Template/report_template.html")

html_out = template.render({"report": report_dict})

pdf_name = "../results/report_testing.pdf"
HTML(string=html_out).write_pdf(pdf_name,
                                stylesheets=[CSS("../ref/Template/report.css"),
                                             CSS("../ref/Template/bootstrap.css")])
