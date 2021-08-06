rule move:
    input:
        seq_dir = "test/"
    output:
        analysis = directory("results/files"),
    log:
        files = "log/move.log"
    script:
        "../scripts/move_fastq.py"
