configfile: "config_grid.yaml"

include: "rules/centrifuge_wildcard.smk"
include: "rules/amr.smk"
include: "rules/qc.smk"

rule all:
    input:
        expand("results/{sample}/amr/scagaire_gene_summary.tsv", sample = SAMPLES),
        expand("results/{sample}/qc/nanostat_summary.txt", sample = SAMPLES)    

	
