rule findViral:
    input:
        fq ="results/{sample}/{time}_hours/unclassified/{sample}_{time}_hours_unclassified.fastq"
    output:
        raw = "results/{sample}/{time}_hours/viral/centrifuge_viral_raw.tsv",
        report = "results/{sample}/{time}_hours/viral/centrifuge_viral_abundance.tsv"
    shell:
        """
        centrifuge -p 4 --mm -x {config[parameters][centrifuge][index][viral]} -q {input.fq} -S {output.raw} \
        --report-file {output.report}
        """

rule parseViral:
    input:
        raw = "results/{sample}/{time}_hours/viral/centrifuge_viral_raw.tsv"
    output:
        report = "results/{sample}/{time}_hours/viral/centrifuge_viral_report.tsv",
        read = "results/{sample}/{time}_hours/viral/read_assignments.tsv"
    script:
        "../scripts/parse_viral.py"


