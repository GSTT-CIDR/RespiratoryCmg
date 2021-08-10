"""Given a time (default 2hrs), this script will take the fastq files from the beginning of the run
 until the end time parameter and move them into a separate folder for analysis"""

from datetime import datetime, timedelta
import os
import shutil

TIME= int(snakemake.wildcards.time)
FQ_DIR = snakemake.input.seq_dir
OUT_DIR = snakemake.output.analysis


def move_files(fq_dir=".", end_dir="results/files/"):
    """
    Moves fastq files from sequencing results folder to analysis folders
    Parameters
    ----------
    fq_dir: str

    end_dir: str

    Returns: None
    -------

    """
    file_names = os.listdir(fq_dir)
    for f in file_names:
        if f.endswith(".fastq"):
            shutil.copy(os.path.join(fq_dir, f), end_dir)
    return file_names


def main():
    time = TIME
    print("Time cut off set to {} minutes".format(time))
    end_time = datetime.now() + timedelta(seconds=time)

    still_run = True
    while still_run:
        if datetime.now() > end_time:
            # move_files(config["fastq_dir"])
            print("Moving Files")
            log = move_files(fq_dir=FQ_DIR, end_dir = OUT_DIR)
            still_run = False
    # with open(snakemake.log[0], "w") as f:
    #     for item in log:
    #         f.write("%s\n" % item)


if __name__ == "__main__":
    main()


