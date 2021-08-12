rule human_reads:
    input:
        "results/{sample}/{time}_minutes/files/{sample}_{time}_minutes_merged.fastq"
    output:
        "results/{sample}/{time}_minutes/host/host_mapped.sam"
    shell:
        "minimap2 -x map-ont --secondary no -a {config[parameters][hg38][index]} {input} > {output}"

rule micro_fastq:
    input:
        "results/{sample}/{time}_minutes/host/host_mapped.sam"
    output:
        "results/{sample}/{time}_minutes/microbial/hg38_unmapped.fastq"
    shell:
        "samtools fastq -f 4 {input} > {output}"


