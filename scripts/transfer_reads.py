import glob
import datetime
from Bio import SeqIO
from dateutil.parser import parse as dparse
import pytz
import time
import pyfastx
import shutil
utc=pytz.UTC	
gb = pytz.timezone("GB")


THRESHOLD = float(snakemake.wildcards.time)
FQ_DIR = snakemake.input.seq_dir
OUTFILE = snakemake.output.analysis
LOGFILE = snakemake.log[0]
MOVE = snakemake.config["move"]

def transfer():
    read_files = []
    start_time= gb.localize(datetime.datetime.now())
    cutoff_time = start_time + datetime.timedelta(hours=THRESHOLD)
    max_time = utc.localize(datetime.datetime.min)
    # Adapted from WouterDeCoster answer from Biostars
    sleep_interval = 1 # minutes
    fastq_dict = dict()

    print("Wait interval set to {} minutes".format(sleep_interval))
    KEEP_GOING = True
    while KEEP_GOING:
        current_time = gb.localize(datetime.datetime.now())
        time.sleep(sleep_interval * 60)
        file_list = glob.glob(f"{FQ_DIR}/*")
        to_read = [i for i in file_list if i not in read_files]
        print("Processing files {}".format(to_read))
        for file in to_read:
            # New module pyfastx 
            for name,seq,qual,comment in pyfastx.Fastx(file):
                read_time = dparse(
                    [i for i in comment.split() if i.startswith("start_time")][0].split("=")[1])
                raw = f"@{name} {comment}\n{seq}\n+\n{qual}\n"
                fastq_dict[name] = [read_time, raw]
            # for record in SeqIO.parse(open(file, "rt"), "fastq"):
            #     read_time = dparse(
            #         [i for i in record.description.split() if i.startswith("start_time")][0].split("=")[1])
                # fastq_list.append([read_time, record.format("fastq")])
                if read_time < start_time:
                    start_time = read_time
                    cutoff_time = start_time + datetime.timedelta(hours=THRESHOLD, minutes=sleep_interval) # Added the 5 minute sleep to allow read queue to write to files
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
            print("Cut-off time set to {} (including 5 minute buffer time) with Threshold of {} hours, still running".format(cutoff_time, THRESHOLD))
    print("Reading {} reads".format(len(fastq_dict)))
    passed = 0
    with open(OUTFILE, "w") as of:
        for r_time, fq in fastq_dict.values():
            if r_time < cutoff_time:
                passed += 1
                of.write(fq)
    print("{} reads below threshold".format(passed))

    with open(LOGFILE, "w") as log:
        for f in read_files:
            log.write("{}\n".format(f))

def copy():
    shutil.copyfile(FQ_DIR, OUTFILE)
    with open(LOGFILE, "w") as log:
        log.write(f"{FQ_DIR} - transferred")

def main():
    if MOVE is True:
        transfer()
    else:
        copy()

if __name__ == "__main__":
    main()


