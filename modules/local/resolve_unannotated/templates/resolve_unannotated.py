#!/usr/bin/env python3

"""
Annotates genome segments that are not annotated by iCount segmentation,
when performing segmentation based on a filtered GTF file.
"""


import platform
import argparse
import tempfile
import csv
import pandas as pd
import pybedtools as pbt
import plumbum as pb
import logging


def dump_versions(process_name):
    """Write the software version and process name to a file."""
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")
        out_f.write("    pandas: " + pd.__version__ + "\n")
        out_f.write("    pybedtools: " + pbt.__version__ + "\n")
        out_f.write("    plumbum: " + pb.__version__ + "\n")


def read_gtf(segmentation):
    """Read GTF file and return a pandas DataFrame."""
    df_regions = pd.read_csv(
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
    # Get a set of chromosomes
    chromosomes = set(df_regions.chrom)
    return chromosomes, df_regions


def fai2bed(fai):
    """Convert fasta index to BED6 format and return as pybedtools object."""
    df_chromosomes = pd.read_csv(
        fai,
        sep="\t",
        header=None,
        names=["chrom", "end", "offset", "linebases", "linewidth"],
        dtype={
            "chrom": str,
            "end": int,
            "offset": int,
            "linebases": int,
            "linewidth": int
        },
        )
    df_chromosomes = df_chromosomes[["chrom", "end"]].assign(start=0, name=".", score=0)
    # Get a set of chromosomes
    chromosomes = set(df_chromosomes.chrom)
    # Assign positive strand
    df_chromosomes_p = df_chromosomes.copy()
    df_chromosomes_p["strand"] = "+"
    df_chromosomes_p = df_chromosomes_p[["chrom", "start", "end", "name", "score", "strand"]]
    # Assign negative strand
    df_chromosomes_m = df_chromosomes.copy()
    df_chromosomes_m["strand"] = "-"
    df_chromosomes_m = df_chromosomes_m[["chrom", "start", "end", "name", "score", "strand"]]
    # Combine both strands, convert to pyranges and sort.
    df_chromosomes = pd.concat([df_chromosomes_p, df_chromosomes_m], ignore_index=True)
    bed_chr = pbt.BedTool.from_dataframe(df_chromosomes).sort()
    return chromosomes, bed_chr


def validate_chromosomes(filt_chromosomes, unfilt_chromosomes, index_chromosomes):
    """
    Validate that the chromosomes in the filtered and unfiltered GTF files match the chromosomes in the fasta index.
    """
    # Validate that the set of chromosome matches between filtered and unfiltered GTFs.
    logging.info("Validating that the same chromosomes are represented in filtered and unfiltered GTF and in the fasta index...")
    # Check that all sets filt_chromosomes, unfilt_chromosomes, index_chromosomes are the same.
    if not (filt_chromosomes == unfilt_chromosomes == index_chromosomes):
        logging.error(
            "Mismatch in chromosomes between filtered GTF regions, unfiltered GTF regions, and fasta index. "
            "Check that the sets of chromosomes are the same across both regions files and the fasta index. "
            f"Filtered GTF chromosomes: {filt_chromosomes} "
            f"Unfiltered GTF chromosomes: {unfilt_chromosomes} "
            f"Fasta index chromosomes: {index_chromosomes} "
        )
        raise ValueError(
            "ERROR: Mismatch in chromosomes between filtered GTF regions, unfiltered GTF regions, and fasta index. "
            "Check that the sets of chromosomes are the same across both regions files and the fasta index."
        )


def validate_unfiltered(bed_fai, bed_unfiltered):
    """
    Validate that there are no unannotated regions in the unfiltered GTF file.
    """
    logging.info("Finding whether unannotated regions exist in the unfiltered GTF...")
    unfilt_missing = bed_fai.subtract(bed_unfiltered, s=True, nonamecheck=True).sort()
    if len(unfilt_missing) == 0:
        logging.info("No unannotated regions found in the unfiltered GTF.")
    else:
        logging.error(f"Found {len(unfilt_missing)} unannotated regions in the unfiltered GTF.")
        raise ValueError(
            "ERROR: Unannotated regions found for the unfiltered GTF. "
            "Something went wrong with the ICOUNTMINI_SEGMENT (ICOUNT_SEG_GTF) process, "
            "or the fasta index provided here did not match the index provided to the ICOUNTMINI_SEGMENT (ICOUNT_SEG_GTF) process."
        )

def validate_resolved(bed_fai, bed_complete):
    """
    Validate that no regions remain unannotated after the resolve process.
    """
    bed_missing_c = bed_fai.subtract(bed_complete, s=True, nonamecheck=True).sort()
    if len(bed_missing_c) > 0:
        logging.error(
            f"Found {len(bed_missing_c)} unannotated regions in the RESOLVED iCount segment. "
            "Writing missing segments to a BED file..."
            )
        bed_missing_c.saveas(f"missing_in_resolved_regions.bed")
        raise ValueError(
            "ERROR: Unannotated regions were found in the resolved iCount segment. "
            "Something went wrong with the RESOLVE_UNANNOTATED process."
        )


def main(process_name, unfilt_regs, filt_regs, fai, output):
    # Dump version file
    dump_versions(process_name)

    # Logging
    log_file = f"{output.replace('resolved.gtf', 'resolved.log')}"
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info(f"Command-line arguments: {vars(args)}")

    # Read filtered iCount genomic segment and convert it from GTF to BED format.
    logging.info(f"Reading iCount genomic segmentation (regions) for filtered GTF in {filt_regs}")
    filt_chromosomes, df_regions = read_gtf(filt_regs)
    bed_regions = df_regions.assign(start=df_regions["start"] - 1, score=0)[
        ["chrom", "start", "end", "feature", "score", "strand"]
    ].copy()
    bed_regions = pbt.BedTool.from_dataframe(bed_regions).sort()

    # Read unfiltered iCount genomic segment and convert it from GTF to BED format.
    logging.info(f"Reading iCount genomic segmentation (regions) for unfiltered GTF in {unfilt_regs}")
    unfilt_chromosomes, df_unfiltered = read_gtf(unfilt_regs)
    bed_unfiltered = df_unfiltered.assign(start=df_unfiltered["start"] - 1, score=0)[
        ["chrom", "start", "end", "feature", "score", "strand", "annotations"]
    ].copy()
    bed_unfiltered = pbt.BedTool.from_dataframe(bed_unfiltered).sort()

    # Save bed filtered and unfiltered for diagnostics
    bed_regions.saveas("filtered_regions.bed")
    bed_unfiltered.saveas("unfiltered_regions.bed")

    # Convert fasta index to BED format - one entry spans one chromosome.
    index_chromosomes, bed_fai = fai2bed(fai)

    # Validate the input files
    validate_chromosomes(filt_chromosomes, unfilt_chromosomes, index_chromosomes)

    # Check whether there are unannotated regions in the "unfiltered" GTF file
    validate_unfiltered(bed_fai, bed_unfiltered)

    # Find regions that are unannotated in the iCount genome segmentation.
    logging.info("Getting unannotated regions...")
    bed_missing = bed_fai.subtract(bed_regions, s=True, nonamecheck=True).sort()

    logging.info(f"Found {len(bed_missing)} unannotated genomic regions.")

    # If missing regions are found proceed with annotating them, otherwise return input as output
    if len(bed_missing) > 0:

        # Intersect missing regions with unfiltered segment to get transcript region
        logging.info("Annotating missing regions in iCount regions file using the unfiltered regions file...")
        # Use intersect to split unannotated regions into chunks matching the unfiltered segment.
        intersect = bed_missing.intersect(bed_unfiltered, s=True, nonamecheck=True).sort()
        # Annotate with annotations (column 7) and feature (column 4), this places the feature annotation in the 8-th column
        missingAnnotated = intersect.map(bed_unfiltered, s=True, c=[7, 4], o="collapse", nonamecheck=True).sort()
        # Convert to dataframe
        df_unannotated = pd.read_csv(
            missingAnnotated.fn,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "name", "score", "strand", "annotations", "feature"],
        )

        # Save bed for diagnostics
        df_unannotated[["chrom", "start", "end", "feature", "score", "strand"]] \
            .to_csv("unannotated_in_filtered_regions_IMPUTED.bed", index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE)

        # Convert to GTF format
        df_unannotated = df_unannotated.assign(start=df_unannotated["start"] + 1, source=".", name2=".")
        df_unannotated = df_unannotated[
            ["chrom", "source", "feature", "start", "end", "name", "strand", "name2", "annotations"]
        ]
        # Add missing regions to original iCount segment.
        logging.info("Adding annotated missing regions to iCount segment...")
        df_regions = pd.concat([df_regions, df_unannotated], ignore_index=True)
        logging.info(f"N segment entries: {len(df_regions)}")

        # Validate that no unannotated regions are left on the genome.
        logging.info("Validating that no regions remain unannotated after resolving...")
        bed_complete = df_regions.assign(start=df_regions["start"] - 1, score=0)[
            ["chrom", "start", "end", "feature", "score", "strand", "annotations"]
        ].copy()
        validate_resolved(bed_fai, pbt.BedTool.from_dataframe(bed_complete).sort())

        # Validate that resolved regions have the same feature types as the unfiltered regions.
        reg_types_out = set(df_regions.feature)
        logging.info(f"Types of regions in the resolved output: {reg_types_out}")
        if reg_types_out.issubset(set(df_unfiltered.feature)):
            logging.info("Types of regions in the resolved output correspond to region types in the unfiltered regions file.")
        else:
            logging.error("Types of regions in the resolved output do not correspond to region types in the unfiltered regions file.")
            raise ValueError(
                f"An unexpected error occurred. Types of regions in the resolved output do not correspond to region types in the unfiltered regions file.\n"
                f"Region types in the resolved output: {reg_types_out}\n"
                f"Region types in the unfiltered regions file: {set(df_unfiltered.feature)}"
            )

        # Save the resolved GTF file and sort it.
        with tempfile.NamedTemporaryFile(mode="w") as tmpfile:
            df_regions.to_csv(tmpfile.name, index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE)
            cmd = (pb.cmd.sort["-t\t", "-k1,1", "-k4,4n", tmpfile.name]) > output
            print(cmd())
        logging.info(f"Saved the resolved segment as {output}")
    else:
        # If no missing regions are found in the filtered segment, return input as output.
        logging.info(f"No missing regions found. Saving inputs segment as {output}.")
        df_regions.to_csv(output, index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE)
    return 0


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--unfilt_regs", default="!{unfilt_regs}")
    parser.add_argument("--filt_regs", default="!{filt_regs}")
    parser.add_argument("--fai", default="!{fai}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    main(
        args.process_name,
        args.unfilt_regs,
        args.filt_regs,
        args.fai,
        args.output
    )
