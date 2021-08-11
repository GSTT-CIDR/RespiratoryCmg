import pandas as pd
import mappy as mp
import os
import time
import taxonomy
from collections import Counter
import json

# SNAKEMAKE ARGUMENTS
NODES = snakemake.config["taxonomy"]["nodes"]
NAMES = snakemake.config["taxonomy"]["names"]
CENTRIFUGE_FILE = snakemake.input.file
TAXMETA = snakemake.config["taxonomy"]["speciesTaxMeta"]
FILENAMES = snakemake.config["taxonomy"]["speciesFileNames"]
FASTA = snakemake.input.fasta
GENOME_DIR = snakemake.config["taxonomy"]["refseqDir"]
MULTI_OUTPUT = snakemake.output.multi
FAILED_OUTPUT = snakemake.output.failed
READ_OUTPUT = snakemake.output.read
REPORT_OUTPUT = snakemake.output.report
THRESHOLD = snakemake.wildcards.thresholds




def remove_low_qual_reads(df, threshold = 300):
    """
    Remove hits below the threhold and collects readIDs of failed reads
    """
    above_threshold = df[df["score"] > threshold]
    rejected_reads = df[df["score"] < threshold]["readID"].unique()
    failed = {k:{"alignment": False, "reason": "score below threshold"} for k in rejected_reads}
    return above_threshold, failed


def species_id(taxID, tax, rank = "species"):
    """checks if taxID is a below species rank and returns its species ID"""
    taxobj = tax.parent(str(taxID), at_rank = rank)
    if taxobj is None:
        return tax.node(str(taxID))
    if taxobj.rank == rank:
        return taxobj
    return None


def get_name(taxID, tax):
    name = tax.node(str(taxID))
    if name is None:
        name = "NA"
        return name
    else:
        return name.name


def split_hits(df, tax):
    """
    Split centrifuge output into two categories: Unique hits and multi-match hits.
    Multi-match hits are further checked to see if hits are to the same taxonomy id/rank
    (ranks below species bought up to species)
    Also filters out reads below threshold
    Parameters
    ----------
    dataFrame
        centrifuge output as pandas DataFrame object
    Returns
    -------
    unique_dict
        dictionary in format;  {readID: taxID}
    multi_dicy
        dictionary in format; {readID:[taxIDs]}
    """

    # Split hits
    unique_df = df[df["numMatches"] == 1]
    multi_df = df[df["numMatches"] != 1]

    unique_dict = dict()
    multi_dict = dict()

    zipped = zip(unique_df["readID"], unique_df["taxID"])
    for r, t in zipped:
        unique_dict[r] = species_id(str(t), tax)

    # Check if multiple hits for read belong to same taxID i.e. Are actually unique hit
    multi_group = multi_df.groupby("readID")
    for read, dat in multi_group:
        tax_nodes = [species_id(str(t),tax) for t in dat["taxID"].unique() if species_id(str(t),tax) is not None]
        tax_ids = set(t.id for t in tax_nodes)
        if len(tax_ids) == 1:
            unique_dict[read] = tax_nodes[0]
        else:
            multi_dict[read] = tax_nodes

    # Bring taxIDS with ranks below species up to level

    return unique_dict, multi_dict


def create_fasta_dict(fasta = FASTA):
    fasta_dict = {}
    for name, seq, qual in mp.fastx_read(fasta, read_comment=False):
        fasta_dict[name] = seq
    return fasta_dict


def create_tax_dict(taxmeta, filenames):
    tax_dict = dict()
    meta_df = pd.read_csv(taxmeta, sep="\t")
    file_df = pd.read_csv(filenames, sep="\t")
    for key, value in zip(meta_df["taxid"], meta_df["species_taxid"]):
        try:
            fileName = file_df.loc[file_df["species_taxid"].isin([value]), "file_name"].values[0]
        except IndexError:
            fileName = "NA"
        tax_dict[key] = {"species_tax": value,
                        "fileName": fileName}
    for key, value in zip(file_df["species_taxid"], file_df["file_name"]):
        tax_dict[key] = {"species_tax": key,
                        "fileName": value}
    return tax_dict


def make_aligner(tax_dict, genome_dir, tax_id):
    if tax_id not in tax_dict.keys():
        return None
    path = os.path.join(genome_dir, tax_dict[tax_id]["fileName"])
    if os.path.isfile(path):
        aligner = mp.Aligner(path, preset="map_ont", best_n=1)
        return aligner
    else:
        return None


def get_fasta(read_id, fasta_dict):
    fasta = fasta_dict[read_id]
    return fasta


def alignment(fasta, aligner):
    mapping = aligner.map(fasta)
    try:
        hit = next(mapping)
        res = {"hit": hit.ctg,
               "alignment_length": hit.blen,
               "total_matches": hit.mlen,
               "mapping_quality": hit.mapq,
               "identity": round(hit.mlen / hit.blen, 3),
               "alignment": True}
    except StopIteration:
        res = {"alignment": False,
               "reason": "No alignment found"}
    return res


def get_read_alignments(tax_ids, genome_dir, tax_dict, fasta, tax):
    hit_dict = dict()
    for tax_id in tax_ids:
        node = tax.node(str(tax_id))
        if node.rank == "species":
            tax_id = int(tax_id)
            aligner = make_aligner(tax_dict, genome_dir, tax_id)
            if aligner is None:
                hit_dict[tax_id] = {"alignment": False,
                                    "reason": "No reference genome"}
            else:
                hit_dict[tax_id] = alignment(fasta, aligner)
        else:
            hit_dict[tax_id] = {"alignment": False,
                                "reason": "tax above species level: {}".format(node.rank)}
    return hit_dict


def pick_best_alignment(align_dict):
    top = None
    best_identity = 0
    for k, v in align_dict.items():
        if v["alignment"] is True:
            if v["identity"] > best_identity:
                top = k
                best_identity = v["identity"]
    if top is None:
        return False
    else:
        return top


def main():
    start = time.time()
    print("loading data")
    tax = taxonomy.Taxonomy.from_ncbi(snakemake.config["taxonomy"]["nodes"], snakemake.config["taxonomy"]["names"])
    df = pd.read_csv(snakemake.input.file, sep="\t")
    taxDict = create_tax_dict(snakemake.config["taxonomy"]["speciesTaxMeta"], snakemake.config["taxonomy"]["speciesFileNames"])
    fasta_dict = create_fasta_dict()

    print("data loaded.... starting counts")
    print("Threshold is set to: {}".format(THRESHOLD))

    above_threshold, failed = remove_low_qual_reads(df, int(THRESHOLD))
    report_dict, multi_hit_dict = split_hits(above_threshold, tax)

    print("Starting multi-match resolution")
    # Multi match sort
    multi_res = dict()
    for key, value in multi_hit_dict.items():
        fasta = get_fasta(key, fasta_dict)
        taxIDs = [int(t.id) for t in value]
        try:
            multi_res[key] = get_read_alignments(taxIDs, GENOME_DIR, taxDict, fasta, tax)
        except KeyError:
            failed[key] = taxIDs

    for k, v in multi_res.items():
        res = pick_best_alignment(v)
        if res:
            report_dict[k] = tax.node(str(res))
        else:
            failed[k] = v

    print("writing files")

    # Mutli match read assignments
    with open(MULTI_OUTPUT, "w") as f:
        json.dump(multi_res, f)

    # Failed reads
    with open(FAILED_OUTPUT, "w") as f:
        json.dump(failed, f)

    # Read assignments
    read_df = pd.DataFrame.from_dict(report_dict, orient="index").reset_index()
    read_df.columns = ["readID", "TaxID"]
    read_df["TaxID"] = read_df["TaxID"].apply(lambda x: x.id if (x is not None) else "NA")
    read_df["Organism"] = read_df["TaxID"].apply(lambda x: get_name(x, tax))  ## swap around
    read_df.to_csv(READ_OUTPUT, index=False, sep="\t")

    report_values = Counter([int(i.id) for i in report_dict.values()])
    total_counts = sum(report_values.values())

    # Report output
    report_df = pd.DataFrame.from_dict(report_values, orient="index").reset_index()
    report_df.columns = ["Tax_ID", "Counts"]
    report_df["Organism"] = report_df["Tax_ID"].apply(lambda x: tax.node(str(x)).name)
    report_df["Percentage"] = report_df["Counts"] / total_counts * 100
    report_df = report_df.sort_values(by="Percentage", ascending=False)
    report_df = report_df[["Organism", "Tax_ID", "Counts", "Percentage"]]
    report_df.to_csv(REPORT_OUTPUT, index=False, sep="\t")

    end = time.time()
    print("Time to run is {} seconds".format(end - start))


if __name__ == "__main__":
    main()

