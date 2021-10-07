rule nanostat:
    input:
        "results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
    output:
        "results/{sample}/{time}_hours/qc/nanostat_summary.txt",
 
    shell:
        "NanoStat --fastq {input} -n {output}"


