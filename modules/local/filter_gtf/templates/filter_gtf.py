#!/usr/bin/env python3

"""Filter GTF file for optimal clipseq execution. Filters GENCODE or ENSEMBL genomic annotation in GTF format."""

import platform
import argparse
from sys import exit
import csv
import pandas as pd


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")
        out_f.write("    pandas: " + pd.__version__ + "\n")


def main(process_name, gtf, output):
    # Dump version file
    dump_versions(process_name)

    # Parse gtf file into pandas dataframe.
    print("Reading annotation file.")
    annotation = pd.read_csv(
        gtf,
        sep="\t",
        names=["chrom", "source", "feature", "start", "end", "name", "strand", "name2", "annotations"],
        header=None,
        comment="#",
        dtype={
            "chrom": str,
            "source": str,
            "feature": str,
            "start": int,
            "end": int,
            "name": str,
            "strand": str,
            "name2": str,
            "annotations": str,
        },
    )

    # Filter transcripts by basic tag
    print("Number of entries in input annotation:", len(annotation))

    # Check if annotation contains tag "basic"
    print("Checking for basic flag...")

    basic = annotation["annotations"].str.contains("basic", regex=True)
    if basic.any():
        print("Basic flag available.")

        nbasic = basic.value_counts()[True]
        print(f"{nbasic} entries flagged as basic.")

        annotation = annotation.loc[
            annotation["annotations"].str.contains('tag "basic"') | (annotation["feature"] == "gene"), :
        ]
        print('Number of entries after filtering for tag "basic":', len(annotation))

        # Filter annotation gene-by-gene by transcript level support (TSL), to keep higher confidence transcripts where possible (TSL1 and 2)
        df_TSL = annotation.loc[
            annotation["annotations"].str.contains(
                'transcript_support_level "1|transcript_support_level "2', regex=True
            ),
            :,
        ]
        gene_ids = df_TSL["annotations"].str.split(";", n=1, expand=True)[0].unique().tolist()
        print("Number of genes that contain TSL1 or TSL2 transcripts:", len(gene_ids))

        # Keeping only TSL1 and TSL2 entries for genes that contain them, discardig other entries (no TSL information or TSL3-5)
        print("Filtering out low-confidence transcripts.")
        df_t = annotation.loc[
            (annotation["feature"] != "gene") & (annotation["annotations"].str.contains("|".join(gene_ids))), :
        ]
        df_t = df_t.loc[
            ~df_t["annotations"].str.contains('transcript_support_level "1"|transcript_support_level "2"', regex=True)
        ]
        annotation.drop(index=df_t.index, inplace=True)

        print("Number of entries in filtered annotation.", len(annotation))
        print("Saving filtered gtf file.")
    else:
        print('No tag "basic". Returning input annotation as output. Exiting.')
    annotation.to_csv(output, header=None, index=None, sep="\t", quoting=csv.QUOTE_NONE)


if __name__ == "__main__":

    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--gtf", default="!{gtf}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    main(args.process_name, args.gtf, args.output)
