rule compile_report:
    input:
        centrifuge = "results/{sample}/{time}_hours/centrifuge/centrifuge_report.tsv",
        centrifuge_raw = "results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv",
        viral = "results/{sample}/{time}_hours/viral/viral_target_report.tsv",
        amr_summary= "results/{sample}/{time}_hours/amr/scagaire_gene_summary.tsv",
        amr_report="results/{sample}/{time}_hours/amr/scagaire_report.tsv",
        qc="results/{sample}/{time}_hours/qc/nanostat_summary.txt",
        stats = "results/{sample}/{time}_hours/host/{sample}_{time}_hours_map_stats.txt"
    output:
        "reports/{sample}/{sample}_{time}_hours_report.pdf"

    script:
        "../scripts/generate_report.py"

rule transfer_qnap:
    input:
        "reports/{sample}/{sample}_{time}_hours_report.pdf"
    output:
        encrypt = "reports/encrypted/{sample}/{sample}_{time}_hours_report_encrypt.pdf",
        transfer = "results/{sample}/{time}_hours/transfer/transferred.txt"
    shell:
        """
        pdftk {input} output {output.encrypt} user_pw cidr22
        rsync -r {output.encrypt} qnap://mnt/flavia/metagenomics/pilot/reports/ > {output.transfer}
        """
