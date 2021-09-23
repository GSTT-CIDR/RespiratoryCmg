rule compile_report:
    input:
        centrifuge = "results/{sample}/{time}_minutes/centrifuge/centrifge_report.tsv",
        amr_summary= "results/{sample}/{time}_minutes/amr/scagaire_gene_summary.tsv",
        amr_report="results/{sample}/{time}_minutes/amr/scagaire_report.tsv",
        qc="results/{sample}/{time}_minutes/qc/nanostat_summary.txt"
    output:
        "results/{sample}/{time}_minutes/report/{sample}_report.pdf"

    script:
        "../scripts/generate_report.py"
