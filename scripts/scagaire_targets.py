import pandas as pd

df = pd.read_csv(snakemake.input[0], sep="\t")
top_hits = set(df[df["Percentage"] > 1]["Organism"])
scagaire_file = snakemake.config["amr"]["scagaire"]
with open(scagaire_file, "r") as f:
    scagaire_genes = {line.rstrip() for line in f}

top_parse = top_hits.intersection(scagaire_genes)

with open(snakemake.output[0], "w") as f:
    for hit in top_parse:
        f.write("{}\n".format(hit))


