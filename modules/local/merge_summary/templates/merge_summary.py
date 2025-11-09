#!/usr/bin/env python3

"""Merge pre-mapped results with genome-computed iCount Summary"""

import platform
import argparse
from sys import exit
import pandas as pd
import os


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")
        out_f.write("    pandas: " + pd.__version__ + "\n")


def extract_cdna_from_bed(file_path):
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep="\t", header=None)
    # Sum the values in the 5th column (index 4)
    total_sum = df[4].sum()
    return total_sum


def adjust_summary_file(file_path, number_cdnas_premapped):
    # Read the file into a DataFrame
    df = pd.read_csv(file_path, sep="\t", header=0)
    # Add the new values
    new_row = ["premapped RNA", "NA", number_cdnas_premapped, 0]
    # Append the new row using loc indexer
    df.loc[len(df)] = new_row
    # Correct the percentages
    # Calculate the total cDNA #
    total_cDNA = df["cDNA #"].sum()
    # Update the cDNA % column
    df["cDNA %"] = (df["cDNA #"] / total_cDNA) * 100
    # Create the output file name
    base_name, extension = os.path.splitext(os.path.basename(file_path))
    output_file = base_name + "_premapadjusted" + extension
    print(f"Saving to: {output_file}")  # Debugging
    # Write the updated DataFrame to the new file
    df.to_csv(output_file, sep="\t", index=False)


def main(processname, summaries):
    type, subtype, gene, cdna = summaries.split(" ")
    
    # Dump version file
    dump_versions(processname)

    # Get number of cDNAs
    number_cdnas_premapped = extract_cdna_from_bed(cdna)
    print("Number of cDNAs premapped: " + str(number_cdnas_premapped))

    adjust_summary_file(type, number_cdnas_premapped)
    adjust_summary_file(subtype, number_cdnas_premapped)
    adjust_summary_file(gene, number_cdnas_premapped)


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--processname", default="!{process_name}")
    parser.add_argument("--summaries", default="!{summaries}")
    args = parser.parse_args()

    print("Process name: " + args.processname)
    print("Summaries: " + args.summaries)
    main(args.processname, args.summaries)
    
