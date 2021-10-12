rule move:
    input:
        seq_dir = lambda wildcards: sample_table.path[wildcards.sample]
    output:
        analysis = "results/{sample}/{time}_hours/files/{sample}_{time}_hours_merged.fastq"
    log:
        "results/{sample}/{time}_hours/log/file_move.log"
    script:
        "../scripts/transfer_reads.py"
