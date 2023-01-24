import pandas as pd
import taxonomy
from collections import Counter


# SNAKEMAKE ARGUMENTS
NODES = snakemake.config["taxonomy"]["nodes"]
NAMES = snakemake.config["taxonomy"]["names"]
TARGETS = snakemake.config["viral"]["targets"]
CENTRIFUGE_FILE = snakemake.input.raw
READ_OUTPUT = snakemake.output.read
REPORT_OUTPUT = snakemake.output.report
TARGETS_OUTPUT = snakemake.output.targets

# NODES = "ref/refseq/taxonomy/nodes.dmp"
# NAMES = "ref/refseq/taxonomy/names.dmp"


def species_id(taxID, tax, rank = "species"):
    """
    checks if taxID is a below taxonomy rank specified and returns that rank.

    Parameters
    ----------
    taxID
    tax
    rank

    Returns
    -------

    """
    taxobj = tax.parent(str(taxID), at_rank = rank)
    if taxobj is None:
        return tax.node(str(taxID))
    if taxobj.rank == rank:
        return taxobj
    return None

def split_hits(df, tax):
    """
    Split centrifuge output into two categories: Unique hits and multi-match hits.
    Multi-match hits are further checked to see if hits are to the same taxonomy id/rank
    (ranks below species bought up to species)
    Also filters out reads below threshold
    Parameters
    ----------
    df
        centrifuge output as pandas DataFrame object
    tax
        taxonomy object from module taxonomy
    Returns
    -------
    unique_dict
        dictionary in format;  {readID: taxID}
    multi_dict
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


def lca(multi_dict, tax):
    lca_dict = dict()
    for k,v in multi_dict.items():
        lca = v[0]
        for i in range(1,len(v)):
            lca = tax.lca(lca.id, v[i].id)
        lca_dict[k] = lca
    return lca_dict


def get_name(taxID, tax):
    name = tax.node(str(taxID))
    if name is None:
        name = "NA"
        return name
    else:
        return name.name

def get_genus(taxID,tax):
    try:
        genus = tax.parent(str(taxID), at_rank="genus")
        return int(genus.id)
    except:
        return "None"



     

def main():
    tax = taxonomy.Taxonomy.from_ncbi(NODES, NAMES)
    df = pd.read_csv(CENTRIFUGE_FILE, sep="\t")
    df = df[df["taxID"] > 0]
    if df.shape[0] != 0:
        report_dict, multi_hit_dict = split_hits(df, tax)
        lca_dict = lca(multi_hit_dict, tax)
        report_dict.update(lca_dict)
        read_df = pd.DataFrame.from_dict(report_dict, orient="index").reset_index()
        read_df.columns = ["readID", "TaxID"]
        read_df["TaxID"] = read_df["TaxID"].apply(lambda x: x.id if (x is not None) else "NA")
        read_df["Organism"] = read_df["TaxID"].apply(lambda x: get_name(x,tax)) ## swap around
        read_df.to_csv(READ_OUTPUT, index=False, sep="\t")

        report_values = Counter([int(i.id) if i is not None else "No ID" for i in report_dict.values()])
        del report_values[9606]
        total_counts = sum(report_values.values())

        #Report output
        report_df = pd.DataFrame.from_dict(report_values, orient="index").reset_index()
        report_df.columns = ["TaxID", "Counts"]
        report_df["Organism"] = report_df["TaxID"].apply(lambda x: tax.node(str(x)).name if (tax.node(str(x)) is not None) else "ARGOS ISOLATE")
        report_df["Percentage"] = round(report_df["Counts"] / total_counts * 100, 3)
        report_df = report_df.sort_values(by="Percentage", ascending=False)
        report_df = report_df[["Organism", "TaxID", "Counts", "Percentage"]]
        report_df.to_csv(REPORT_OUTPUT, index=False, sep="\t")

        #Target output
        target = pd.read_csv(TARGETS, sep="\t")
        target_ids = target["TaxID"].values.tolist()
        report_df["Genus"] = report_df["TaxID"].apply(lambda x: get_genus(x, tax))
        target_df = report_df[report_df["Genus"].isin(target_ids)]
        target_df.to_csv(TARGETS_OUTPUT, index=False, sep="\t")
        


    else:
        fail_df = pd.DataFrame({"Organisms": "N/A", "TaxID": "N/A", "Counts": "N/A": "Percentage": "N/A"})
        fail_df.to_csv(READ_OUTPUT, sep="\t")
        fail_df.to_csv(REPORT_OUTPUT, sep="\t")
        fail_df.to_csv(TARGETS_OUTPUT, sep="\t")

if __name__ == "__main__":
    main()