#!/usr/env/bin python3

import pandas as pd
import argparse 
import os


def get_mlst_scheme(args):
    mlst_list = args.mlst
    target = args.target
    genus = target.split()[0]
    directory = args.dir
    df = pd.read_csv(mlst_list, sep="\t")
    tmp = df[df["Organism"].str.contains(target)]["Folder"].values
    if len(tmp) == 1:
        return os.path.join(directory, tmp[0])
    tmp = df[(df["Organism"].str.contains(genus)) & (df["Rank"].str.contains("Genus"))].values
    if len(tmp) == 1:
        return os.path.join(directory, tmp[0])
    return None

def main(args):
    mlst_path = get_mlst_scheme(args)
    print(mlst_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--target", default="unclassified", help="Target organisms. WARNING: Substring style matching (Default: unclassified)")
    parser.add_argument("-d", "--dir", required=True, help="Directory of mlst")
    parser.add_argument("-m", "--mlst", required=True, help="Path to mlst scheme list")
    args = parser.parse_args()
    main(args)