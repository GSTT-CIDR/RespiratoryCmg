import pandas as pd
import csv
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
INTERVAL = float(snakemake.wildcards.time)
CFG_PATH = snakemake.input.centrifuge
AMR_SUMMARY = snakemake.input.amr_summary
AMR_REPORT = snakemake.input.amr_report
QC_PATH = snakemake.input.qc
OUTPUT = str(snakemake.output)
REPORT_HTML = snakemake.config["pdf"]["html"]
REPORT_CSS = snakemake.config["pdf"]["css"]
BOOTSTRAP_CSS = snakemake.config["pdf"]["bootstrap"]
SAMPLE_TABLE = snakemake.config["samples"]
SAMTOOLS_STAT = snakemake.input.stats
THRESHOLD = snakemake.config["abundance_threshold"]
TARGETS = snakemake.config["targets"]



def convert_bp(size):
    size = float(str(size).replace(",", ""))
    for x in ["bp", "Kb", "Mb", "Gb", "Tb"]:
        if size < 1000.0:
            if x == "bp":
                return "{:.0f} {}".format(size, x)
            else:
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
            qc_dict["nano_reads"] = dat
        elif line.startswith("Total bases"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["total_bp"] = convert_bp(dat)
    return qc_dict


def samtools_stats(path):
    with open(path, "r") as f:
        data = f.read().splitlines()
    total_reads = data[0].split()[0]
    human_reads = data[4].split()[0]
    stats = {"total_reads" : total_reads,
             "human_reads": human_reads}
    return stats


def patient_info(path, id):
    df = pd.read_csv(path, sep = "\t")
    df = df[df["Lab_ID"] == str(id)]
    sample_dict = {"Lab_ID" : df["Lab_ID"].values[0],
                   "Experiment": df["Experiment"].values[0],
                   "Sample_ID": df["Sample_ID"].values[0],
                   "Barcode": df["Barcode"].values[0],
                   "Sample_Type": df["Sample_Type"].values[0],
                   "Patient_ID": df["Patient_ID"].values[0]}
    return sample_dict

def is_target(s, target_file = TARGETS):
    with open(TARGETS, "r") as f:
        targets = [t.strip() for t in f]
    if s["Organism"] in targets:
        return ["color: red"] * 3
    else:
        return ["color: black"] * 3

# Change this from hardcoding threshold
def cfg_to_html(path, threshold = THRESHOLD, target_file = TARGETS):
    exceptions = "Aspergillus|Candida"
    cfg_dict = dict()
    df = pd.read_csv(path, sep="\t")
    df = df.round(2)
    df = df.drop(columns=["Tax_ID"])
    full = df.copy()
    s = full.style.apply(is_target, axis=1)
    s = s.format(precision=2)
    cfg_dict["cfg_full"] = df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    # cfg_dict["cfg_full"] = s.set_table_attributes('class="table table-striped"').hide_index().render()
    cfg_dict["micro_reads"] = sum(df["Counts"])
    ic = df[df["Organism"] == "Jonesia denitrificans"]
    if ic.empty:
        cfg_dict["ic"] = "NA/NA"
    else:
        cfg_dict["ic"] = "{}/{}%".format(int(ic["Counts"]), float(ic["Percentage"]))
    above_df = df.copy()
    above_df = above_df[(above_df["Percentage"] > threshold) |(above_df["Organism"].str.contains(exceptions) ]
    cfg_dict["cfg_top"] = above_df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    return cfg_dict


def amr_summary(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["Species", "Gene", "Hits"])
    res = df[df["Hits"] > 2]
    if res.empty:
        res = res.append({"Species": "No genes above threshold", "Gene": "NA", "Hits": "NA"}, ignore_index=True)
    return res.to_html(classes="table table-striped", border=0, justify="left", index=False)


def amr_report(path, summary_path):
    with open(path, "r") as f:
        dat = list(filter(None, [i for i in csv.reader(f, delimiter="\t")]))
    column_names = ["SPECIES", "SEQUENCE", "START", "END", "GENE", "%COVERAGE"]
    master_df = pd.DataFrame(columns = column_names)
    if dat[0][0].startswith("#FILE"):
        with open(summary_path, "r") as f:
            species = [i for i in csv.reader(f, delimiter="\t")][0][0]
        df = pd.read_csv(path, sep="\t", usecols=["SEQUENCE", "START", "END", "GENE", "%COVERAGE"])
        df["SPECIES"] = species
        return df[column_names].to_html(classes="table table-striped", border=0, justify="left", index=False)
    elif dat[0][0].startswith("Results for species"):
        coord = []
        species = []
        pos = []
        no_res = []
        for num, line in enumerate(dat):
            if line[0].startswith("Results for species"):
                species.append(line[1])
                if len(pos) == 1:
                    pos.append(num)
                    coord.append(pos)
                    pos = []
            if line[0].startswith("#FILE"):
                pos.append(num)
            if line[0].startswith("No results"):
                no_res.append(line[0])
        if len(pos) == 1:
            pos.append(len(dat))
            coord.append(pos)
        if len(coord) == 0:
            df = pd.DataFrame(columns=["SPECIES", "RESULTS"])
            for s, n in zip(species, no_res):
                df = df.append({"SPECIES": s,
                                "RESULTS": n}, ignore_index=True)
            return df.to_html(classes="table table-striped", border=0, justify="left", index=False)
        for s, c in zip(species,coord):
            tmp = dat[c[0]:c[1]]
            df = pd.DataFrame(tmp[1:], columns=tmp[0])
            df["SPECIES"] = s
            master_df = master_df.append(df[column_names])
        return master_df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    else:
        return "No Results"


report_dict = {"time": str(INTERVAL) + " hrs",
               "title": "Clinical metagenomics report",
               "date": datetime.now()}

report_dict.update(samtools_stats(SAMTOOLS_STAT))
report_dict.update(patient_info(SAMPLE_TABLE, SAMPLE))
report_dict.update(cfg_to_html(CFG_PATH))
report_dict.update(summary_qc(QC_PATH))
report_dict["amr_summary"] = amr_summary(AMR_SUMMARY)
report_dict["amr_report"] = amr_report(AMR_REPORT, AMR_SUMMARY)


env = Environment(loader=FileSystemLoader(".")) # Change to "." for grid
# declare our jinja template
template = env.get_template(REPORT_HTML)

html_out = template.render({"report": report_dict})

pdf_name = OUTPUT
HTML(string=html_out).write_pdf(pdf_name,
                                stylesheets=[CSS(REPORT_CSS),
                                             CSS(BOOTSTRAP_CSS)])
