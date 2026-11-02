#!/usr/bin/env python3
"""
Join the output from both runs of multimap_class_binning - run once at the specific
category level (snRNA | TE | tRNA), and then once at the region level (3UTR | CDS | intron).
"""

import pandas as pd
import argparse
import os
import yaml

parser = argparse.ArgumentParser(
    prog='Merge assigned reads from binning',
    description='Creates bam from unassigned reads - reads which were discarded by previous iteration of binning due to '
    'ambiguity.'
)

parser.add_argument('-a',
                    help = 'Dataframe a, original output')

parser.add_argument('-b',
                    help = 'Dataframe b, secondary output')

parser.add_argument('-premap',
                    nargs = "+",
                    help = "List of samtools view -c bams with counts of reads")

args = parser.parse_args()


TE_counts = pd.read_csv(args.a)

region_counts = pd.read_csv(args.b)

## Create df from bowtie logs

linenum = 0
d = []

for countfile in args.premap:
        sample_id = os.path.basename(countfile).replace(".bam", "").replace(".count_premapped", "")
        with open(countfile, "rt") as myfile:
            count = int(myfile.read().strip())
        d.append({"sample":sample_id, "pre-mapped":count})

df = pd.DataFrame(d)

print(df)

merged = pd.merge(TE_counts, region_counts, on="sample", how="left")

print(merged)

merged_withpremap = pd.merge(merged, df, on = "sample", how = "left")

print(merged_withpremap)

merged_withpremap.to_csv("multimap_binning.csv", index=False)

## Create custom yaml with enough colours for categories --------------------------------------------------------------------------

# codes from: DOI: 10.18637/jss.v090.c01
preset_hex = ("#4B4749", "#E2E1DF", "#CC0000", "#F700FB", "#1CFC00", "#1688FB", "#FCB600", "#DF5F9F", "#16FECA",
              "#8FB22E", "#0DCDEE",
              "#8238A2", "#894216", "#9D0DFD", "#FE168E", "#FCBDFC", "#22865C", "#FC7F00", "#FCDE9F", "#FC00C4",
              "#1C568A", "#FFB6BA",
              "#6E6300", "#F3E522", "#ABAAFE", "#0DC14B", "#FE7A80", "#F676EA", "#873D71", "#ADF816", "#B2F1B7",
              "#491CC4", "#89EAE2",
              "#D26DFE", "#97AFD0", "#F5A071")

# Starting from 1 as first col is always "sample"
cats = list(merged_withpremap.columns[1:])
colours = preset_hex[0:len(cats)]

# Output config customised to number of categoties for multiqc
cats_config = {
      "custom_plot_config": {
            "binning_bargraph": {
                cat: {"color": col} for cat, col in zip(cats, colours)
          },
    }
}

with open("colours_mqc_config.yaml", "w") as f:
    yaml.dump(cats_config, f, default_flow_style=False)
