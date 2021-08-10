rule human_reads:
    input:
        "results/{sample}/{time}_minutes/files/"
    output:
        "results/{sample}/{time}_minutes/host/host_mapped.sam"
    shell:
        "minimap2 -x map-ont --secondary no -a {config[parameters][hg38][index]} {input}*.fastq > {output}"

rule micro_fastq:
    input:
        "results/{sample}/{time}_minutes/host/host_mapped.sam"
    output:
        micro_reads = "results/{sample}/{time}_minutes/microbial/hg38_unmapped.fastq"
    shell:
        "samtools fastq -f 4 {input} > {output}"


