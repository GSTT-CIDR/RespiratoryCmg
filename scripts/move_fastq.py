"""Given a time (default 2hrs), this script will take the fastq files from the beginning of the run
 until the end time parameter and move them into a separate folder for analysis"""

from datetime import datetime, timedelta
import os
import shutil




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
    if not os.path.exists(end_dir):
        print("path doesn't exist. trying to make")
        os.makedirs(end_dir)
    for f in file_names:
        shutil.move(os.path.join(fq_dir, f), end_dir)
    return file_names


def main():
    time = snakemake.config["parameters"]["time"]
    # time = 1
    print("Time cut off set to {} minutes".format(time))
    end_time = datetime.now() + timedelta(minutes=time)
    # end_time = datetime.now() + timedelta(hours=snakemake.input[0])
    still_run = True
    while still_run:
        if datetime.now() > end_time:
            # move_files(config["fastq_dir"])
            print("Moving Files")
            log = move_files(fq_dir=snakemake.input.seq_dir, end_dir = snakemake.output[0])
            still_run = False
    with open(snakemake.log[0], "w") as f:
        for item in log:
            f.write("%s\n" % item)



if __name__ == "__main__":
    main()


