rule run_centrifuge:
    input:
        "results/{sample}/{time}_minutes/microbial/hg38_unmapped.fastq"
    output:
        raw = "results/{sample}/{time}_minutes/centrifuge/centrifuge_raw.tsv",
        report = temp("results/{sample}/{time}_minutes/centrifuge/centrifuge_report_raw.tsv")
    shell:
        "centrifuge -p 4 --min-hitlen 25 --mm -x {config[parameters][centrifuge][index]} -q {input} -S {output.raw} \
        --report-file {output.report}"




rule parse_centrifuge:
    input:
        file = "results/{sample}/{time}_minutes/centrifuge/centrifuge_raw.tsv",
        fasta = "results/{sample}/{time}_minutes/microbial/hg38_unmapped.fastq"
    output:
        report = "results/{sample}/{time}_minutes/centrifuge/centrifuge_report.tsv",
        read = "results/{sample}/{time}_minutes/centrifuge/read_assignments.tsv",
        failed = "results/{sample}/{time}_minutes/centrifuge/failed_reads.json",
        multi = "results/{sample}/{time}_minutes/centrifuge/multi_read.json"
    # conda:
    #     "../env/centrifuge.yml"
    script:
        "../scripts/centrifuge_multi_match.py"



