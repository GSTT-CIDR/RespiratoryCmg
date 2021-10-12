rule compile_report:
    input:
        centrifuge = "results/{sample}/{time}_hours/centrifuge/centrifuge_report.tsv",
        amr_summary= "results/{sample}/{time}_hours/amr/scagaire_gene_summary.tsv",
        amr_report="results/{sample}/{time}_hours/amr/scagaire_report.tsv",
        qc="results/{sample}/{time}_hours/qc/nanostat_summary.txt",
        stats = "results/{sample}/{time}_hours/host/{sample}_{time}_hours_map_stats.txt"
    output:
        "reports/{sample}/{sample}_{time}_hours_report.pdf"

    script:
        "../scripts/generate_report.py"
