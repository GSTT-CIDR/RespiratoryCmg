rule fastq2fasta:
    input:
        "files/{sample}.fastq"
    output:
        temp("results/{sample}/{sample}.fasta")
    shell:
        "seqtk seq -a {input} > {output}"


rule abricate:
    input:
        "results/{sample}/{sample}.fasta"
    output:
        amr = "results/{sample}/amr/amr_results.tsv"
    shell:
        "abricate {input} > {output}"


rule top_centrifuge:
    input:
        "results/{sample}/centrifuge/centrifuge_report.tsv"
    output:
        "results/{sample}/amr/centrifuge_top_hits.tsv"
    script:
        "../scripts/scagaire_targets.py"


rule scagaire:
    input:
        amr_res = "results/{sample}/amr/amr_results.tsv",
        top_hit = "results/{sample}/amr/centrifuge_top_hits.tsv"
    output:
        "results/{sample}/amr/scagaire_gene_summary.tsv"
    shell:
        """
        t=$(cat {input.top_hit} | paste -sd "," -)
        scagaire "$t" {input.amr_res} -s {output}
        if [ ! -f {output} ]; then
            echo "No species in database" > {output}
        fi
        """


