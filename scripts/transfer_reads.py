import glob
import datetime
from Bio import SeqIO
from dateutil.parser import parse as dparse
import pytz
import time
utc=pytz.UTC

THRESHOLD = int(snakemake.wildcards.time)
FQ_DIR = snakemake.input.seq_dir
OUTFILE = snakemake.output.analysis
LOGFILE = snakemake.log[0]

def main():
    read_files = []
    start_time= utc.localize(datetime.datetime.now())

    # Adapted from WouterDeCoster answer from Biostars
    sleep_interval = 10 # minutes
    fastq_list = []

    print("Wait interval set to {} minutes".format(sleep_interval))
    KEEP_GOING = True
    while KEEP_GOING:
        current_time = utc.localize(datetime.datetime.now())
        time.sleep(sleep_interval * 60)
        file_list = glob.glob("{}/*".format(file_path))
        to_read = [i for i in file_list if i not in read_files]
        print("Processing files {}".format(to_read))
        for file in to_read:
            for record in SeqIO.parse(open(file, "rt"), "fastq"):
                read_time = dparse(
                    [i for i in record.description.split() if i.startswith("start_time")][0].split("=")[1])
                fastq_list.append([read_time, record.format("fastq")])
                if read_time < start_time:
                    start_time = read_time
                    cutoff_time = start_time + datetime.timedelta(hours=THRESHOLD)
                if read_time > max_time:
                    max_time = read_time
        read_files.extend(to_read)
        buffer_time = cutoff_time + datetime.timedelta(minutes=15)
        if max_time > cutoff_time:
            print("Past time threshold: writing relevant reads to file")
            KEEP_GOING = False
        elif current_time > buffer_time:
            print("Time exceeded: Writing files")
            KEEP_GOING = False
        else:
            print("Threshold of {} hours not met, still running".format(THRESHOLD))

    with open(OUTFILE, "w") as of:
        for r_time, fq in fastq_list:
            if r_time < cutoff_time:
                of.write(fq)

    with open(LOGFILE, "w") as log:
        for f in read_files:
            log.write("{}\n".format(f))
