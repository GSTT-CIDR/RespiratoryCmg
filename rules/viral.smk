rule findViral:
    input:
        fq ="results/{sample}/{time}_hours/unclassified/{sample}_{time}_hours_unclassified.fastq"
    output:
        raw = "results/{sample}/{time}_hours/viral/centrifuge_raw.tsv",
        report = "results/{sample}/{time}_hours/viral/centrifuge_report_raw.tsv"

    shell:
        """centrifuge -p 4 --mm -x {config[parameters][centrifuge][index][viral]} -q {input.fq} -S {output.raw} \
        --report-file {output.report}
        """