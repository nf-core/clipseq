#!/usr/bin/env python3

"""Annotates genome segments that are not annotated by iCount segmentation."""


import platform
import argparse
from sys import exit
import tempfile
import csv
import pandas as pd
import pybedtools as pbt
import plumbum as pb


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")
        out_f.write("    pandas: " + pd.__version__ + "\n")
        out_f.write("    pybedtools: " + pbt.__version__ + "\n")
        out_f.write("    plumbum: " + pb.__version__ + "\n")


def read_gtf(segmentation):
    df_segment = pd.read_csv(
        segmentation,
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
    return df_segment


def fai2bed(fai):
    df_chromosomes = pd.read_csv(fai, sep="\t", header=None, names=["chr", "end", "offset", "linebases", "linewidth"])
    df_chromosomes = df_chromosomes[["chr", "end"]].assign(start=0, name=".", score=0)
    df_chromosomes_p = df_chromosomes.copy()
    df_chromosomes_p["strand"] = "+"
    df_chromosomes_p = df_chromosomes_p[["chr", "start", "end", "name", "score", "strand"]]

    df_chromosomes_m = df_chromosomes.copy()
    df_chromosomes_m["strand"] = "-"
    df_chromosomes_m = df_chromosomes_m[["chr", "start", "end", "name", "score", "strand"]]

    df_chromosomes = pd.concat([df_chromosomes_p, df_chromosomes_m], ignore_index=True)
    bed_chr = pbt.BedTool.from_dataframe(df_chromosomes).sort()
    return bed_chr


def main(process_name, segmentation, filt_segmentation, annotation, fai, output, genic_other):
    # Dump version file
    dump_versions(process_name)

    # Read filtered iCount genomic segment and convert it from GTF to BED format.
    print("Reading genomic segmentation.")
    df_segment = read_gtf(filt_segmentation)
    bed_segment = df_segment.assign(start=df_segment["start"] - 1, score=0)[
        ["chrom", "start", "end", "feature", "score", "strand"]
    ]
    bed_segment = pbt.BedTool.from_dataframe(bed_segment).sort()
    # Read unfiltered iCount genomic segment and convert it from GTF to BED format.
    df_unfiltered = read_gtf(segmentation)
    bed_unfiltered = df_unfiltered.assign(start=df_unfiltered["start"] - 1, score=0)[
        ["chrom", "start", "end", "feature", "score", "strand", "annotations"]
    ]
    bed_unfiltered = pbt.BedTool.from_dataframe(bed_unfiltered).sort()

    # Convert fasta index to BED format - one entry spans one chromosome.
    bed_fai = fai2bed(fai)

    # Read annotation GTF, keep only gene-level entries and convert it to BED format.
    print("Getting gene-level annotation...")
    df_annotation = read_gtf(annotation)
    df_annotation = df_annotation.loc[df_annotation["feature"] == "gene"]
    bed_annotation = df_annotation.assign(start=df_annotation["start"] - 1, score=0)[
        ["chrom", "start", "end", "annotations", "score", "strand"]
    ]
    bed_annotation = pbt.BedTool.from_dataframe(bed_annotation).sort()

    # Find regions that are unannotated in the iCount genome segmentation.
    print("Getting unannotated regions...")
    bed_missing = bed_fai.subtract(bed_segment, s=True, nonamecheck=True).sort()
    print(f"Found {len(bed_missing)} unannotated genomic regions.")
    # Use intersect to split unnanotated regions
    intersect = bed_missing.intersect(bed_unfiltered, s=True, nonamecheck=True).sort()

    print("Annotating regions with gene information...")
    if genic_other == "false":
        # Intersect missing regions with unfiltered segment to get transcript region
        print("Annotating missing regions in iCount segment with transcript regions...")
        # Annotate with annotations (column 7) and feature (column 4)
        missingAnnotated = intersect.map(bed_unfiltered, s=True, c=[7, 4], o="collapse", nonamecheck=True).sort()
        df_unnanotated = pd.read_csv(
            missingAnnotated.fn,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "name", "score", "strand", "annotations", "feature"],
        )
        df_unnanotated = df_unnanotated.assign(start=df_unnanotated["start"] + 1, source=".", name2=".")
    else:
        print('Annotationg missing regions in iCount segment as "genic_other".')
        # Annotate with annotations (column 7), feature is genic_other
        missingAnnotated = intersect.map(bed_unfiltered, s=True, c=7, o="collapse", nonamecheck=True).sort()
        df_unnanotated = pd.read_csv(
            missingAnnotated.fn,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "name", "score", "strand", "annotations"],
        )
        # Feature is genic_other.
        df_unnanotated = df_unnanotated.assign(
            feature="genic_other", start=df_unnanotated["start"] + 1, source=".", name2="."
        )
        # Find regions that are missing from main annotation and annotate them with gene entries from annotations gtf file
        # Get complement of raw segment to find missing regions
        regsMissingFromMain = bed_fai.subtract(bed_unfiltered, s=True, nonamecheck=True).sort()
        # Annotate them from annotation gtf (gene level annotation) and format entries to replace "gene_type" with "biotype"
        regsMissingFromMain = regsMissingFromMain.map(bed_annotation, s=True, c=4, o="collapse", nonamecheck=True)
        dfMissingFromMain = pd.read_csv(
            regsMissingFromMain.fn,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "name", "score", "strand", "annotations"],
        )
        dfMissingFromMain["annotations"] = dfMissingFromMain["annotations"].apply(
            lambda x: str(x).replace("gene_type", "biotype")
        )
        dfMissingFromMain = dfMissingFromMain.assign(
            feature="genic_other", start=dfMissingFromMain["start"] + 1, source=".", name2="."
        )
        dfMissingFromMain = dfMissingFromMain[
            ["chrom", "source", "feature", "start", "end", "name", "strand", "name2", "annotations"]
        ]
        # Combine all missing regions
        df_unnanotated = pd.concat([df_unnanotated, dfMissingFromMain])
    df_unnanotated = df_unnanotated[
        ["chrom", "source", "feature", "start", "end", "name", "strand", "name2", "annotations"]
    ]
    # Add missing regions to original iCount segment.
    print("Adding annotated missisng regions to iCount segment...")
    df_segment = pd.concat([df_segment, df_unnanotated], ignore_index=True)
    print("N segment entries:", len(df_segment))
    # Sort GTF segment and write it to file
    if genic_other == "true":
        identifier = "genic_other"
    else:
        identifier = "annotated"
    with tempfile.NamedTemporaryFile(mode="w") as tmpfile:
        df_segment.to_csv(tmpfile.name, index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE)
        cmd = (pb.cmd.sort["-t\t", "-k1,1", "-k4,4n", tmpfile.name]) > output
        print(cmd())
    print(f"Saved the segment as {output}")
    return 0


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--segmentation", default="!{segmentation}")
    parser.add_argument("--filt_segmentation", default="!{filt_segmentation}")
    parser.add_argument("--annotation", default="!{annotation}")
    parser.add_argument("--fai", default="!{fai}")
    parser.add_argument("--genic_other", default="!{genic_other}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    main(
        args.process_name,
        args.segmentation,
        args.filt_segmentation,
        args.annotation,
        args.fai,
        args.output,
        args.genic_other,
    )
