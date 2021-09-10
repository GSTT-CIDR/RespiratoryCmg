rule move:
    input:
        seq_dir = lambda wildcards: sample_table.path[wildcards.sample]
    output:
        analysis = "results/{sample}/{time}_minutes/files/{sample}_{time}_minutes_merged.fastq"
    log:
        "results/{sample}/{time}_minutes/log/file_move.log"
    script:
        "../scripts/transfer_reads.py"
