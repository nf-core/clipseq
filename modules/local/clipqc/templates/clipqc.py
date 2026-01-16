#!/usr/bin/env python3

"""Filter GTF file for optimal clipseq execution. Filters GENCODE or ENSEMBL genomic annotation in GTF format."""

import platform
import argparse
import os
from sys import exit
import re
import pybedtools as pbt
import numpy as np
import pandas as pd


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")
        out_f.write("    numpy: " + np.__version__ + "\n")
        out_f.write("    pandas: " + pd.__version__ + "\n")
        out_f.write("    pybedtools: " + pbt.__version__ + "\n")


def main(process_name):
    # Dump version file
    dump_versions(process_name)

    # ==========
    # Mapping
    # ==========

    # First get Bowtie2 logs
    bowtie_logs = sorted(["premap/" + f for f in os.listdir("premap") if f.endswith(".out")])

    ncrna = dict((key, []) for key in ["exp", "input_reads", "ncrna_reads"])

    for bowtie_log in bowtie_logs:
        with open(bowtie_log, "r") as logfile:
            exp = re.sub(".out", "", os.path.basename(bowtie_log))

            lines = logfile.readlines()

            input_reads = int(re.findall(r"\d+", lines[0])[0])
            output_reads = [i for i in lines if "# reads that failed to align" in i]
            output_reads = int(re.findall(r"\d+", output_reads[0])[0])

            ncrna["exp"].append(exp)
            ncrna["input_reads"].append(input_reads)
            ncrna["ncrna_reads"].append(input_reads - output_reads)

    ncrna_df = pd.DataFrame(ncrna)
    print(ncrna_df)

    # Next get STAR logs
    star_logs = sorted(["mapped/" + f for f in os.listdir("mapped") if f.endswith(".Log.final.out")])

    genome = dict((key, []) for key in ["exp", "genome_reads", "unmapped_reads"])

    for star_log in star_logs:
        with open(star_log, "r") as logfile:
            exp = re.sub(".Log.final.out", "", os.path.basename(star_log))

            lines = logfile.readlines()

            input_reads = [i for i in lines if "Number of input reads" in i]
            input_reads = int(re.findall(r"\d+", input_reads[0])[-1])

            genome_reads = [i for i in lines if "Uniquely mapped reads number" in i]
            genome_reads = int(re.findall(r"\d+", genome_reads[0])[-1])

            unmapped_reads = input_reads - genome_reads

            genome["exp"].append(exp)
            genome["genome_reads"].append(genome_reads)
            genome["unmapped_reads"].append(input_reads - genome_reads)

    genome_df = pd.DataFrame(genome)

    # Combine the two
    mapping_df = pd.merge(ncrna_df, genome_df, on="exp")
    mapping_df.to_csv("mapping_metrics.tsv", sep="\t", index=False)

    # Subset for MultiQC plots
    mapping_df.loc[:, ["exp", "ncrna_reads", "genome_reads", "unmapped_reads"]].to_csv(
        "mapping.tsv", sep="\t", index=False
    )

    print(mapping_df)
    print("\n\n")

    # ==========
    # Deduplication
    # ==========

    dedup_logs = sorted(["collapse/" + f for f in os.listdir("collapse") if f.endswith(".log")])
    dedup = dict((key, []) for key in ["exp", "input_reads", "output_reads", "mean_umis", "ratio"])

    for dedup_log in dedup_logs:
        with open(dedup_log, "r") as logfile:
            exp = re.sub(".genomeUnique.dedup_UMICollapse.log", "", os.path.basename(dedup_log))

            lines = logfile.readlines()

            input_reads = [i for i in lines if "Number of input reads" in i]
            input_reads = int(re.findall(r"\d+", input_reads[0])[-1])

            output_reads = [i for i in lines if "Number of reads after deduplicating" in i]
            output_reads = int(re.findall(r"\d+", output_reads[0])[-1])

            mean_umis = [i for i in lines if "Average number of UMIs per alignment position" in i]
            mean_umis = float(re.findall(r"\d+", mean_umis[0])[-1])

            dedup["exp"].append(exp)
            dedup["input_reads"].append(input_reads)
            dedup["output_reads"].append(output_reads)
            dedup["mean_umis"].append(mean_umis)
            dedup["ratio"].append(round(input_reads / output_reads, 2))

    dedup_df = pd.DataFrame(dedup)
    dedup_df.to_csv("dedup_metrics.tsv", sep="\t", index=False)

    # Subset for MultiQC plots
    dedup_df.loc[:, ["exp", "input_reads", "output_reads"]].to_csv("dedup_reads.tsv", sep="\t", index=False)
    dedup_df.loc[:, ["exp", "mean_umis"]].to_csv("dedup_mean_umis.tsv", sep="\t", index=False)
    dedup_df.loc[:, ["exp", "ratio"]].to_csv("dedup_ratio.tsv", sep="\t", index=False)

    print(dedup_df)
    print("\n\n")

    # ==========
    # Crosslinks
    # ==========

    def read_bed(filename):
        df = pd.read_table(
            filename,
            header=None,
            names=["chr", "start", "end", "name", "score", "strand"],
            dtype={
                "chr": str,
                "start": int,
                "end": int,
                "name": str,
                "score": float,
                "strand": str,
            },
        )
        return df

    # First get xlink bed files
    xlinks_files = sorted(["xlinks/" + f for f in os.listdir("xlinks") if f.endswith(".bed")])

    xlinks = dict((key, []) for key in ["exp", "number_of_cDNAs", "total_positions_of_crosslinks", "ratio"])

    for xlinks_file in xlinks_files:
        xlinks_df = read_bed(xlinks_file)

        exp = re.sub(".genome.xl.bed", "", os.path.basename(xlinks_file))
        number_of_cDNAs = xlinks_df["score"].sum()
        total_positions_of_crosslinks = xlinks_df.shape[0]
        ratio = number_of_cDNAs / total_positions_of_crosslinks

        xlinks["exp"].append(exp)
        xlinks["number_of_cDNAs"].append(number_of_cDNAs)
        xlinks["total_positions_of_crosslinks"].append(total_positions_of_crosslinks)
        xlinks["ratio"].append(round(number_of_cDNAs / total_positions_of_crosslinks, 2))

    xlinks_metrics_df = pd.DataFrame(xlinks)
    xlinks_metrics_df.to_csv("xlinks_metrics.tsv", sep="\t", index=False)

    # Subset for MultiQC plots
    xlinks_metrics_df.loc[:, ["exp", "number_of_cDNAs", "total_positions_of_crosslinks"]].to_csv(
        "xlinks_counts.tsv", sep="\t", index=False
    )
    xlinks_metrics_df.loc[:, ["exp", "ratio"]].to_csv("xlinks_ratio.tsv", sep="\t", index=False)

    print(xlinks_metrics_df)
    print("\n\n")

    # ==========
    # Summary cDNA
    # ==========
    summary_type_files = sorted(["summary_type/" + f for f in os.listdir("summary_type") if f.endswith(".tsv")])
    # summary_data = dict((key, []) for key in ["exp", "Type", "Length", "Total_number_cDNA", "perc_number_cDNA"])
    # summary_data = dict((key, []) for key in ["exp", "Type", "perc_number_cDNA"])
    summary_data = []

    for summary_type_file in summary_type_files:

        exp = re.sub(".summary_type_premapadjusted.tsv", "", os.path.basename(summary_type_file))
        df = pd.read_csv(summary_type_file, sep="\t")

        exp_data = {"exp": exp}

        # Add percentage for each RNA type
        for _, row in df.iterrows():
            rna_type = row["Type"]
            exp_data[rna_type] = row["cDNA %"]

        summary_data.append(exp_data)

    summary_df = pd.DataFrame(summary_data)
    desired_columns = ["exp", "CDS", "ncRNA", "intergenic", "premapped RNA"]
    existing_cols = [col for col in desired_columns if col in summary_df.columns]
    summary_df = summary_df[existing_cols]
    summary_df.to_csv("summary_type_metrics.tsv", sep="\t", index=False)

    # ==========
    # Summary subtype cDNA
    # ==========

    summary_subtype_files = sorted(
        ["summary_subtype/" + f for f in os.listdir("summary_subtype") if f.endswith(".tsv")]
    )
    summary_subtypedata = []

    for summary_subtype_file in summary_subtype_files:

        exp = re.sub(".summary_subtype_premapadjusted.tsv", "", os.path.basename(summary_subtype_file))
        df = pd.read_csv(summary_subtype_file, sep="\t")

        exp_data = {"exp": exp}

        # Add percentage for each RNA type
        for _, row in df.iterrows():
            rna_type = row["Subtype"]
            exp_data[rna_type] = row["cDNA %"]

        summary_subtypedata.append(exp_data)

    summary_subtype_df = pd.DataFrame(summary_subtypedata)
    desired_columns = [
        "exp",
        "CDS mRNA",
        "ncRNA rRNA",
        "ncRNA snRNA",
        "ncRNA snoRNA",
        "ncRNA tRNA",
        "ncRNA ncRNA",
        "intergenic",
        "premapped RNA",
    ]
    existing_cols = [col for col in desired_columns if col in summary_subtype_df.columns]
    summary_subtype_df = summary_subtype_df[existing_cols]
    summary_subtype_df.to_csv("summary_subtype_metrics.tsv", sep="\t", index=False)

    # ==========
    # Peaks
    # ==========

    peakcallers = ["icount", "paraclu", "clippy", "pureclip"]

    def get_peaks_metrics(peakcaller):
        peak_files = sorted([peakcaller + "/" + f for f in os.listdir(peakcaller) if f.endswith(".bed")])

        peaks = dict(
            (key, [])
            for key in [
                "exp",
                "peakcaller",
                "total_xlinks",
                "total_xlinksites",
                "total_peaks",
                "median_peak_width",
                "mean_peak_width",
                "xlinks_in_peaks",
                "xlinks_in_peaks_percent",
                "xlinksites_in_peaks",
                "xlinksites_in_peaks_percent",
                "peaks_xlinksite_coverage_percent",
            ]
        )

        for peak_file in peak_files:
            peaks_df = read_bed(peak_file)

            if peakcaller == "icount":
                exp = re.sub(".peaks.bed", "", os.path.basename(peak_file))
            elif peakcaller == "paraclu":
                exp = re.sub(".genome.peaks.clustered.simplified.bed", "", os.path.basename(peak_file))
            elif peakcaller == "clippy":
                exp = re.sub(
                    ".genome.peaks_rollmean10_minHeightAdjust1.0_minPromAdjust1.0_minGeneCount5_Peaks.bed",
                    "",
                    os.path.basename(peak_file),
                )
            elif peakcaller == "pureclip":
                exp = re.sub("_pureclip.peaks.bed", "", os.path.basename(peak_file))

            xlinks_df = read_bed(xlinks_files[0])
            expanded_xlinks_df = xlinks_df.loc[xlinks_df.index.repeat(xlinks_df.score)].reset_index(drop=True)

            # Get metrics
            total_xlinks = xlinks_df["score"].sum()
            total_xlinksites = xlinks_df.shape[0]
            total_peaks = peaks_df.shape[0]

            total_peak_width = sum((peaks_df["end"] - peaks_df["start"]).tolist())
            mean_peak_width = round(np.mean((peaks_df["end"] - peaks_df["start"]).tolist()), 2)
            median_peak_width = np.median((peaks_df["end"] - peaks_df["start"]).tolist())

            xlinks_bed = pbt.BedTool.from_dataframe(xlinks_df)
            expanded_xlinks_bed = pbt.BedTool.from_dataframe(expanded_xlinks_df)
            peaks_bed = pbt.BedTool.from_dataframe(peaks_df)

            xlinks_in_peaks = peaks_bed.intersect(expanded_xlinks_bed, s=True, c=True)
            xlinks_in_peaks = sum([int(c[-1]) for c in xlinks_in_peaks])

            xlinksites_in_peaks = peaks_bed.intersect(xlinks_bed, s=True, c=True)
            xlinksites_in_peaks = sum([int(c[-1]) for c in xlinksites_in_peaks])

            if total_peak_width > 0:
                peaks_xlinksite_coverage_percent = round((xlinksites_in_peaks / total_peak_width) * 100, 2)
            else:
                peaks_xlinksite_coverage_percent = np.nan

            # Write out
            peaks["exp"].append(exp)
            peaks["peakcaller"].append(peakcaller)
            peaks["total_xlinks"].append(total_xlinks)
            peaks["total_xlinksites"].append(total_xlinksites)
            peaks["total_peaks"].append(total_peaks)
            peaks["median_peak_width"].append(median_peak_width)
            peaks["mean_peak_width"].append(mean_peak_width)
            peaks["xlinks_in_peaks"].append(xlinks_in_peaks)
            peaks["xlinks_in_peaks_percent"].append(round(xlinks_in_peaks / total_xlinks * 100, 2))
            peaks["xlinksites_in_peaks"].append(xlinksites_in_peaks)
            peaks["xlinksites_in_peaks_percent"].append(round(xlinksites_in_peaks / total_xlinksites * 100, 2))
            peaks["peaks_xlinksite_coverage_percent"].append(peaks_xlinksite_coverage_percent)

        peaks_metrics_df = pd.DataFrame(peaks)
        return peaks_metrics_df

    peaks_metrics = [get_peaks_metrics(pc) for pc in peakcallers]
    peaks_metrics_df = pd.concat(peaks_metrics)
    peaks_metrics_df.to_csv("peaks_metrics.tsv", sep="\t", index=False)

    # Subset for MultiQC plots
    total_peaks_df = peaks_metrics_df.pivot_table(index="exp", columns="peakcaller", values="total_peaks")
    total_peaks_df.to_csv("total_peaks.tsv", sep="\t", index=True)

    xlinks_in_peaks_percent_df = peaks_metrics_df.pivot_table(
        index="exp", columns="peakcaller", values="xlinks_in_peaks_percent"
    )
    xlinks_in_peaks_percent_df.to_csv("xlinks_in_peaks.tsv", sep="\t", index=True)

    xlinksites_in_peaks_percent_df = peaks_metrics_df.pivot_table(
        index="exp", columns="peakcaller", values="xlinksites_in_peaks_percent"
    )
    xlinksites_in_peaks_percent_df.to_csv("xlinksites_in_peaks.tsv", sep="\t", index=True)

    peaks_xlinksite_coverage_percent_df = peaks_metrics_df.pivot_table(
        index="exp", columns="peakcaller", values="peaks_xlinksite_coverage_percent"
    )
    peaks_xlinksite_coverage_percent_df.to_csv("peaks_xlinksite_coverage.tsv", sep="\t", index=True)

    print(peaks_xlinksite_coverage_percent_df)
    print("\n\n")


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    args = parser.parse_args()

    main(args.process_name)
