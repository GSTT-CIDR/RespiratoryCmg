rule fastq2fasta:
    input:
        "results/{sample}/{time}_minutes/microbial/hg38_unmapped.fastq"
    output:
        temp("results/{sample}/{time}_minutes/{sample}.fasta")
    shell:
        "seqtk seq -a {input} > {output}"


rule abricate:
    input:
        "results/{sample}/{time}_minutes/{sample}.fasta"
    output:
        amr = "results/{sample}/{time}_minutes/amr/amr_results.tsv"
    shell:
        "abricate --mincov 90 {input} > {output}"


rule top_centrifuge:
    input:
        "results/{sample}/{time}_minutes/centrifuge/centrifuge_report.tsv"
    output:
        "results/{sample}/{time}_minutes/amr/centrifuge_top_hits.tsv"
    script:
        "../scripts/scagaire_targets.py"


rule scagaire:
    input:
        amr_res = "results/{sample}/{time}_minutes/amr/amr_results.tsv",
        top_hit = "results/{sample}/{time}_minutes/amr/centrifuge_top_hits.tsv"
    output:
        summary = "results/{sample}/{time}_minutes/amr/scagaire_gene_summary.tsv",
        report = "results/{sample}/{time}_minutes/amr/scagaire_report.tsv"
    shell:
        """
        t=$(cat {input.top_hit} | paste -sd "," -)
        scagaire "$t" {input.amr_res} -s {output.summary} -o {output.report}
        """
