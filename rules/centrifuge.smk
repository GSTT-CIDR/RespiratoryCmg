rule run_centrifuge:
    input:
        "results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
        # "results/{sample}/{time}_hours/microbial/hg38_unmapped.fastq"
    output:
        raw = "results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv",
        report = temp("results/{sample}/{time}_hours/centrifuge/centrifuge_report_raw.tsv")
    shell:
        "centrifuge -p 4 --min-hitlen 25 --mm -x {config[parameters][centrifuge][index]} -q {input} -S {output.raw} \
        --report-file {output.report}"




rule parse_centrifuge:
    input:
        file = "results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv",
        fasta = "results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
    output:
        report = "results/{sample}/{time}_hours/centrifuge/centrifuge_report.tsv",
        read = "results/{sample}/{time}_hours/centrifuge/read_assignments.tsv",
        failed = "results/{sample}/{time}_hours/centrifuge/failed_reads.json",
        multi = "results/{sample}/{time}_hours/centrifuge/multi_read.json"
    script:
        "../scripts/centrifuge_multi_match.py"



