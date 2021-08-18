#SAMPLES, = glob_wildcards("files/{sample}.fastq")


rule run_centrifuge:
    input:
        "files/{sample}.fastq"
    output:
        raw = "results/{sample}/centrifuge/centrifuge_raw.tsv",
        report = temp("results/{sample}/centrifuge/centrifuge_report_raw.tsv")
    shell:
        "centrifuge -p 4 --min-hitlen 25 --mm -x {config[parameters][centrifuge][index]} -q {input} -S {output.raw} \
        --report-file {output.report}"




rule parse_centrifuge:
    input:
        file = "results/{sample}/centrifuge/centrifuge_raw.tsv",
        fasta = "files/{sample}.fastq"
    output:
        report = "results/{sample}/centrifuge/centrifuge_report.tsv",
        read = "results/{sample}/centrifuge/read_assignments.tsv",
        failed = "results/{sample}/centrifuge/failed_reads.json",
        multi = "results/{sample}/centrifuge/multi_read.json"
    # conda:
    #     "../env/centrifuge.yml"
    script:
        "../scripts/centrifuge_multi_match.py"



