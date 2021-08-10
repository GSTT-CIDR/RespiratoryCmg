rule move:
    input:
        seq_dir = lambda wildcards: sample_table.path[wildcards.sample]
    output:
        analysis = directory("results/{sample}/{time}_minutes/files/")
    script:
        "../scripts/move_fastq.py"
