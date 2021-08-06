"""
Runs centrifuge and creates report of all Species which account for at least a threshold percentage of the microbial
reads (default 1%)


REDUNDANT WITH SNAKEMAKE/ NEXTFLOW
"""
import subprocess
import logging
import os
import pandas as pd
import argparse
from collections import defaultdict
import io
import time

DEFAULT_REPORT = "centrifuge_report.tsv"
DEFAULT_PERCENTAGE = 1
DEFAULT_CSCORE = 2503
DEFAULT_PATH = ""

parser = argparse.ArgumentParser(description="")
parser.add_argument("path", type = str, default= DEFAULT_PATH,
                    help = "Absolute path to sequencing run with de-multiplexed FastQ files")

## Logging testing
# logging.basicConfig(filename='example.log', level=logging.DEBUG)
# log = logging.getLogger(__name__)
# log.debug('This message should go to the log file')
# log.info('So should this')
# log.warning('And this, too')
# log.error('And non-ASCII stuff, too, like Øresund and Malmö')
#
# files = logging.getLogger("files")
# files.debug("hello")
# files.info("Let's see if this works")
# files.warning('I hope it does')

def get_barcode_path(path):
    """

    Parameters
    ----------
    path

    Returns
    -------

    """
    path_dict = defaultdict(list)
    if os.path.isdir(path):
        for path, dirs, file, in os.walk(path):
            for d in dirs:
                for f in file:
                    path_dict[d].append(f)
    return path_dict

#
# files = os.walk("/Users/alderc/1-projects/12-GSST/2-Data/Test/fastq_pass")
# for t in files:
#     print(t)


def get_fastq(path, prefix):
    """
    Extracts all read files within directory that have not been mapped to human genome

    This is currently a workaround whilst we create a workaround to automate the pipeline

    Parameters
    ----------
    path : str
        path to directory which contains fastq files

    prefix : str
        prefix of reads
    Returns
    -------
    fastqFileList: list
        A list of fastq file paths containing prefix
    """
    fastqFileList = list()

    for paths, dirs, files in os.walk(path):
        for f in files:
            if f.endswith(".fastq") and f.startswith(prefix):
                print(f)
                fastqFileList.append(os.path.join(paths,f))

    return fastqFileList


if __name__ == '__main__':
    path = "/Users/alderc/1-projects/12-GSST/2-Data/EPI2ME_WIMP_Test/"
    prefix = "FAL92041"
    centrifuge_index = "/Users/alderc/1-projects/12-GSST/1-Project/2-metagenomics/index_hpvc/hpvc"

    fastq_str = ",".join(get_fastq(path, prefix)) # Main
    # fastq_str = os.path.join(path,"FAL92041_pass_barcode01_82fc1539_0.fastq") # Single file files
    centrifuge_cmd = "centrifuge -p {} --min-hitlen 25 --mm -x {} -q {}".format(4, centrifuge_index, fastq_str)

    start_time = time.time()

    proc = subprocess.Popen(
        centrifuge_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        stdin=subprocess.PIPE,
        shell=True,
        # Aliased by `text=True` in 3.7
        universal_newlines=True,
        )

    out, err = proc.communicate()

    proc.stdout.close()
    end_time = time.time()
    print("centrifuge done in %s seconds" % (end_time - start_time))

    # Read results output of centrifuge - stdout results
    df = pd.read_csv(io.StringIO(out),
                     sep="\t")
    df.to_csv("{}_barcode01_run_centrifuge_phv.tsv".format(time.strftime("%Y%m%d")),
              sep="\t",
              index = False)
