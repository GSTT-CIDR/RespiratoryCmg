#!/usr/env/bin python3

"""
Using centrifuge results, extract reads matching a particular target - Default is 'unclassified' reads in centrifuge,
but is in development to extract reads matching a particular genus
"""

from Bio import SeqIO
import pandas as pd
import argparse
import sys


def get_readIDs(args):
    """
    Extract IDs of reads matching target seqID (contig) of centrifuge results

    Parameters
    ----------
    cfg : 
        path to centrifuge_raw.tsv file
    target : str, optional
        name of target seqID from centrifuge report, by default "unclassified"

    Returns
    -------
    list(str)
        list of readIDs matching target
    """
    cfg = args.cfg
    target = args.target
    df = pd.read_csv(cfg, sep="\t")
    if target == "unclassifed":
        readIDs = df[df["seqID"].str.contains(target)]["readID"].tolist()
    else:
        readIDs = df[df["Organism"].str.contains(target)]["readID"].tolist()
    return readIDs


def extract_reads(args,readIDs):
    """


    Parameters
    ----------
    fastq : str
        path to fastq file
    readIDs : list(str)
        list of readIDs to extract from fastq file
    outfile: str
        name of output fasta file

    Returns
    -------
    file
        fasta file containing target readIDs from fastq
    """
    fastq = args.fastq
    reads = []
    for record in SeqIO.parse(fastq, "fastq"):
        if record.id in readIDs:
            reads.append(record.format("fastq"))
    output = open(args.output, "w") if args.output else sys.stdout
    for r in reads:
        output.write(r)
    output.close()
    return None


def main(args):
    readIDs=get_readIDs(args)
    extract_reads(args,readIDs)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", required=True, help="Fastq file of sample")
    parser.add_argument("-c", "--cfg", required=True, help="Centrifuge data of sample")
    parser.add_argument("-t", "--target", default="unclassified", help="Target organisms. WARNING: Substring style matching (Default: unclassified)")
    parser.add_argument("-o", "--output", nargs="?", default = "", help="name of output fasta (Default: stdout)")
    args = parser.parse_args()
    main(args)


