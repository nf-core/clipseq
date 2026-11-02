# Hannah Jones
# 20/10/2026
# based on gettRNAxlinksites from Charlotte Capitanchik

## Details
# this script takes a bed file intersected with tRNA annotation (and other
# annotations if needs be) and summarises the tRNA anticodon and amino acid level.

# -- Imports ----------------------------------------------------------------------

import argparse
import re
import pandas as pd
import os

# -- Arguments --------------------------------------------------------------------

parser = argparse.ArgumentParser(
    prog = "Summarise tRNA counts",
    description = "Summarises tRNAs by amino acid and anticodon, reporting ambiguity."
)

parser.add_argument('-i', '--input',
                    help = "Bedfile intersected with gtf annotation")

args = parser.parse_args()

# -- Configuration ------------------------------------------------------------------

bedfile = args.input

## ------------------------------BEGIN SCRIPT ---------------------------------------


# -- 1. Create dictionary of tRNAs per read ----------------------------------------
reads_dict = {}

# Adds anticodon, amino acid, isocoder and copy number into dictionary entry for
# each read

with open(bedfile) as bf:
    for line in bf:
        f = line.rstrip("\n").split("\t")
        # check for tRNA
        if f[13] != "tRNA_db":
            continue
        # Make dictionary
        read_name = f[3]
        transcript_summary = re.search(r'transcript_id "([^"]+)"', f[20]).group(1)
        transcript_details = transcript_summary.split("-")
        # take amino acid, anticodon, isocoder, copy number
        mapping = (transcript_details[1], transcript_details[2], transcript_details[3], transcript_details[4])
        mappings = reads_dict.setdefault(read_name, [])
        # multi-mapped reads have multiple entries per read name
            # if statement ensures that a read mapped to the same tRNA > once
            # e.g. if multiple times within one tRNA
            # is only counted as aligning to the one tRNA
        if mapping not in mappings:
            mappings.append(mapping)


## Example output
#'CCAAACACAG:HWI-D00611:179:C7MR1ANXX:5:2215:4632:83486_ATCCCTTCGT', [('Asn', 'GTT', '5'),
# ('Asn', 'GTT', '18'), ('Glu', 'TTC', '4'), ('Gly', 'CCC', '1'), ('Val', 'CAC', '11')]),
# ('TCACCACCCT:HWI-D00611:179:C7MR1ANXX:5:2302:5621:58258_CTTCATCCAT', [('Asn', 'GTT', '5'),
# ('Asn', 'GTT', '18'), ('Glu', 'TTC', '4'), ('Gly', 'CCC', '1'), ('Val', 'CAC', '11')]),

# -------------------------------------------------------------------------------------------

# -- 2. Assign ambiguous and unambiguous reads ----------------------------------------------

entirely_unambiguous = {}
ambiguous_amino_acid = {}
ambiguous_anticodon = {}
ambiguous_isocoder = {}
ambiguous_copy_number = {}

# Sort reads into dictionaries based on number of entries

for reads, mappings in reads_dict.items():
    if len(mappings) == 1:
        entirely_unambiguous[reads] = mappings
    else:
        amino_acids = []
        anticodons = []
        isocoders = []
        copy_numbers = []
        # log all categories
        for i in mappings:
            amino_acids.append(i[0])
            anticodons.append(i[1])
            isocoders.append(i[2])
            copy_numbers.append(i[3])
        # number of different kinds of each category
        amino_acids_set = set(amino_acids)
        anticodons_set = set(anticodons)
        isocoders_set = set(isocoders)
        copy_numbers_set = set(copy_numbers)

        if len(amino_acids_set) > 1:
            ambiguous_amino_acid[reads] = mappings
            # for tRNAs which map to only one amino acid
        else:
            if len(anticodons_set) > 1:
                ambiguous_anticodon[reads] = mappings
                    # for tRNAs which map to 1 amino acid and 1 anticodon
            else:
                if len(isocoders_set) > 1:
                    ambiguous_isocoder[reads] = mappings
                        # for tRNAs which map to 1 amino acid, 1 isocoder, different numbers
                else:
                    if len(copy_numbers_set) > 1:
                        ambiguous_copy_number[reads] = mappings

# -------------------------------------------------------------------------------------------

# -- 3. Format into usable table  -----------------------------------------------------------

# Each read now assigned
unambiguous_resolved = {key: (mappings[0][0], mappings[0][1], mappings[0][2], mappings[0][3]) for key, mappings in entirely_unambiguous.items()}
ambiguous_amino_acid_resolved = {key: ("ambiguous", "ambiguous", "ambiguous", "ambiguous") if len(ambiguous_amino_acid[key]) > 1 else ambiguous_amino_acid[key] for key in ambiguous_amino_acid}
ambiguous_anticodon_resolved = {key: (mappings[0][0], "ambiguous", "ambiguous", "ambiguous") for key, mappings in ambiguous_anticodon.items()}
ambiguous_isocoder_resolved = {key: (mappings[0][0], mappings[0][1], "ambiguous", "ambiguous") for key, mappings in ambiguous_isocoder.items()}
ambiguous_copy_number_resolved = {key: (mappings[0][0], mappings[0][1], mappings[0][2], "ambiguous") for key, mappings in ambiguous_copy_number.items()}

unambiguous_df = pd.DataFrame.from_dict(unambiguous_resolved,
                                        orient="index",
                                        columns=["amino_acid", "anticodon", "isocoder", "copy_number"]).rename_axis("read_name").reset_index()

amino_acid_df = pd.DataFrame.from_dict(ambiguous_amino_acid_resolved,
                                    orient="index",
                                    columns=["amino_acid", "anticodon", "isocoder", "copy_number"]).rename_axis("read_name").reset_index()

anticodon_df = pd.DataFrame.from_dict(ambiguous_anticodon_resolved,
                                    orient="index",
                                    columns=["amino_acid", "anticodon", "isocoder", "copy_number"]).rename_axis("read_name").reset_index()

isocoder_df = pd.DataFrame.from_dict(ambiguous_isocoder_resolved,
                                    orient="index",
                                    columns=["amino_acid", "anticodon", "isocoder", "copy_number"]).rename_axis("read_name").reset_index()

copy_number_df = pd.DataFrame.from_dict(ambiguous_copy_number_resolved,
                                    orient="index",
                                    columns=["amino_acid", "anticodon", "isocoder", "copy_number"]).rename_axis("read_name").reset_index()

all_df = pd.concat([unambiguous_df, amino_acid_df, anticodon_df, isocoder_df, copy_number_df],
ignore_index = True)


# -------------------------------------------------------------------------------------------

# -- 4. Export ------------------------------------------------------------------------------

sample_name = os.path.basename(bedfile).split("_intersect")[0]

all_df.to_csv(f"{sample_name}_tRNA_summary.csv", index=False)
