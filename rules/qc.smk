rule nanostat:
    input:
        "results/{sample}/{time}_minutes/microbial/hg38_unmapped.fastq"
    output:
        "results/{sample}/{time}_minutes/qc/nanostat_summary.txt",
 
    shell:
        "NanoStat --fastq {input} -n {output}"


