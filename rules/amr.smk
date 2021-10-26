rule abricate:
    input:
             "results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
    output:
        amr = "results/{sample}/{time}_hours/amr/amr_results.tsv"
    shell:
        "abricate --mincov 90 --db resfinder {input} > {output}"


rule top_centrifuge:
    input:
        "results/{sample}/{time}_hours/centrifuge/centrifuge_report.tsv"
    output:
        "results/{sample}/{time}_hours/amr/centrifuge_top_hits.tsv"
    script:
        "../scripts/scagaire_targets.py"


rule scagaire:
    input:
        amr_res = "results/{sample}/{time}_hours/amr/amr_results.tsv",
        top_hit = "results/{sample}/{time}_hours/amr/centrifuge_top_hits.tsv"
    output:
        summary = "results/{sample}/{time}_hours/amr/scagaire_gene_summary.tsv",
        report = "results/{sample}/{time}_hours/amr/scagaire_report.tsv"
    shell:
        """
        t=$(cat {input.top_hit} | paste -sd "," -)
        scagaire "$t" {input.amr_res} -s {output.summary} -o {output.report}
        if [ ! -f "{output.summary}" ]; then
            echo "No species in database" > {output.summary}
            echo "No Report" > {output.report}
        fi
        """
