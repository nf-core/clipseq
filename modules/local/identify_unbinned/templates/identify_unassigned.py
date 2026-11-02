#!/usr/bin/env python3
"""
Include only reads discarded due to ambiguous TE/ncRNA identity.

Produce a sorted and indexed .bam file from this to be run through
multimap_class_binning again.
"""
import pysam
import pandas as pd
import argparse
import os

# Argparse

parser = argparse.ArgumentParser(
    prog='Create bam with unassigned reads',
    description='Creates bam from unassigned reads - reads which were discarded by previous iteration of binning due to '
    'ambiguity.'
)

parser.add_argument('-b', '--bam',
                    help = 'List of original bam files, previously run through multimap_class_binning.py',
                    nargs="+")

parser.add_argument('-d', '--discarded',
                    help = 'List of reads discarded by multimap_class_binning.py due to class' \
                    'ambiguity, in .tsv format')

args = parser.parse_args()

bams = args.bam

## Read in discarded read names

discarded_reads =  pd.read_csv(args.discarded, sep='\t')

discarded_reads_grouped = discarded_reads.groupby("sample")

discarded_reads_dict = discarded_reads_grouped["read_name"].apply(set).to_dict()

def extract_sample_name(bam_filename):
    """Sample name - everything before ."""
    return bam_filename.split(".")[0]

for bam in bams:
    # Get ID
    label = os.path.splitext(os.path.basename(bam))[0]

    sample_name = extract_sample_name(label)

    print(sample_name)

    # Open original bam

    orig_bam = pysam.AlignmentFile(bam, "rb")

    # Make sure bam file has discarded reads

    if sample_name not in discarded_reads_dict:
        print(sample_name, "was found to have no discarded reads from initial bucketing")
        continue

    else:
        with pysam.AlignmentFile(f"{label}.discarded_reads.bam", "wb", template=orig_bam) as out:
            for read in orig_bam:
                if read.query_name in discarded_reads_dict[sample_name]:
                    out.write(read)
                else:
                    continue
    orig_bam.close()

    # sort and index

    pysam.sort("-o", f"{label}.discarded_reads_sorted.bam", f"{label}.discarded_reads.bam")
    pysam.index(f"{label}.discarded_reads_sorted.bam")