
SAMPLES, = glob_wildcards("results/files/{samples}.fastq")

rule human_map:
    input:
        rules.move.output.analysis
    output:
        mapped = "results/hg38/human_mapped.sam"
    shell:
        "minimap2 -x map-ont --secondary no -a {config[parameters][hg38][index]} {input}/*.fastq > {output}"

rule micro_fastq:
    input:
        rules.human_map.output.mapped
    output:
        micro_reads = "results/microbial/hg38_unmapped.fastq"
    shell:
        "samtools fastq -f 4 {input} > {output}"


