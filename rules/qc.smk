rule nanostat:
    input:
        "files/{sample}.fastq"
    output:
        "results/{sample}/qc/nanostat_summary.txt",
 
    shell:
        "NanoStat --fastq {input} -n {output}"


