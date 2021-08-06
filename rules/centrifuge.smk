

rule run_centrifuge:
    input:
        "results/microbial/hg38_unmapped.fastq"
    output:
        raw = "results/centrifuge/centrifuge_raw.tsv",
        report = temp("results/centrifuge/centrifuge_report_raw.tsv")
    # conda:
    #     "../env/centrifuge.yml"
    shell:
        "centrifuge -p 4 --min-hitlen 25 --mm -x {config[parameters][centrifuge][index]} -q {input} -S {output.raw} \
        --report-file {output.report}"


rule parse_centrifuge:
    input:
        file = "results/centrifuge/centrifuge_raw.tsv",
        fasta = "results/microbial/hg38_unmapped.fastq"
    output:
        report = "results/centrifuge/centrifuge_report.tsv",
        failed = "results/centrifuge/failed_reads.json"
    # conda:
    #     "../env/centrifuge.yml"
    script:
        "../scripts/centrifuge_multi_match.py"


