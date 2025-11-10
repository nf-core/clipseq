#!/usr/bin/env python3

"""
Validates user-provided transcripts or automatically selects the longest
transcript per gene based on CDS and exon length.

Validates and filters genome annotation (GTF) by the list of selected transcript IDs (.txt).

Outputs:
- A list of selected transcript IDs (.txt)
- A transcript length index file (.fai) with transcript ID and exon length
- A transcript GTF
- A filtered genome annotation GTF containing only entries corresponding to representative transcripts
- A log file documenting all settings, validation and filtering steps
"""

import platform
import argparse
from sys import exit
import pyranges as pr
import pandas as pd
import csv
import logging
import os


def dump_versions(process_name):
    """Write the software version and process name to a file."""
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")
        out_f.write("    pandas: " + pd.__version__ + "\n")
        out_f.write("    pyranges: " + pr.__version__ + "\n")


def parse_gtf(gtf_path):
    """
    Reads a GTF file and checks that every gene has at least one transcript.

    Returns:
    - df_gtf: full parsed GTF as a DataFrame
    """
    df_gtf = pr.read_gtf(gtf_path, as_df=True)
    input_genes = set(df_gtf.gene_id)

    # Identify genes that have at least one associated transcript
    transcript_rows = df_gtf[df_gtf['Feature'] == 'transcript']
    genes_with_transcripts = set(transcript_rows.gene_id)

    # Check for genes with no transcripts
    missing_genes = input_genes - genes_with_transcripts
    if missing_genes:
        top = sorted(missing_genes)[:10]
        more = f" and {len(missing_genes) - 10} more" if len(missing_genes) > 10 else ""
        logging.error(f"Some genes have no transcript entries in the GTF: {top}{more}")
        raise ValueError(
            "ERROR: Some genes in the GTF have no associated transcript entries. "
            "Please ensure the GTF includes transcript-level features for each gene. "
            "Or try to run the pipeline with --skip_filter_gtf enabled, but expect many "
            "processes will FAIL because the GTF format is not compatible."
        )
    else:
        logging.info("GTF passed validation: All genes have associated transcript entries.")

    return df_gtf

def compute_gtf_feature_lengths(gtf_path):
    """
    Reads a GTF file and calculates CDS and exon lengths per transcript.

    Returns:
    - df_gtf: full parsed GTF as a DataFrame (with added 'length' column)
    - df_txlengths: DataFrame with gene_id, transcript_id, cds_length, exon_length
    - input_genes: set of all gene_ids present in the GTF
    """
    df_gtf = parse_gtf(gtf_path)
    df_gtf['length'] = df_gtf['End'] - df_gtf['Start']
    input_genes = set(df_gtf.gene_id)

    # Calculate CDS, exon and unspliced total lengths per transcript
    cds_sums = df_gtf[df_gtf['Feature'] == 'CDS'].groupby(['gene_id', 'transcript_id'])['length'].sum()
    exon_sums = df_gtf[df_gtf['Feature'] == 'exon'].groupby(['gene_id', 'transcript_id'])['length'].sum()
    tx_len = df_gtf[df_gtf['Feature'] == 'transcript'].copy().set_index(['gene_id', 'transcript_id'])['length']

    df_txlengths = pd.concat([cds_sums, exon_sums, tx_len], axis=1, keys=['cds_length', 'exon_length', 'unspliced_length']).reset_index()
    # Fill missing values with 0
    df_txlengths[['cds_length', 'exon_length', 'unspliced_length']] = df_txlengths[['cds_length', 'exon_length', 'unspliced_length']].fillna(0)

    # Error if ANY gene has only transcripts with exon length 0
    zero_exon_genes = df_txlengths.groupby('gene_id')['exon_length'].apply(lambda x: (x == 0).all())
    if zero_exon_genes.any():
        logging.error("Some genes have only transcripts with exon length 0. Genes with zero exon length transcripts: %s", zero_exon_genes[zero_exon_genes].index)
        raise ValueError("Some genes have only transcripts withou exons (exon length 0). "
                         "Please ensure the GTF includes at least one exon feature for each gene. "
                         "Or try to run the pipeline with --skip_filter_gtf enabled, but expect many "
                         "processes will FAIL because the GTF format is not compatible.") # Placing responsibility on the user to provide a GTF compatible with the pipeline.
    else:
        logging.info("All genes have a valid exon feature.") # All genes have non-zero total exon lengths.

    return df_gtf, df_txlengths, input_genes


def load_and_validate_transcripts(df_txlengths, transcript):
    """
    Loads the transcript IDs from the provided file and validates them against the GTF.

    Args:
    - df_txlengths: DataFrame containing the GTF transcript data, incl. a transcript_id column
    - transcript: Path to the file containing the user-provided transcript IDs

    Returns:
    - df_txlengths_filtered (pandas.DataFrame): A filtered DataFrame containing the representative transcript per gene.
    - transcript_ids_sorted: Sorted list of valid transcript IDs
    Raises:
    - ValueError: If the transcript list is empty or any IDs are missing from the GTF
    """
    # Log the transcript file being used
    logging.info(f"Using the transcript IDs in: {transcript}")

    # Load the transcript IDs from the file
    with open(transcript, "r") as file:
        transcript_ids = [line.strip() for line in file]
        # Remove potential empty strings caused by trailing newline
        transcript_ids = [x for x in transcript_ids if x]

    # Check if the transcript file was empty or only contained blank lines
    if len(transcript_ids) == 0:
        logging.error("Transcript list file is empty.")
        raise ValueError(
            "ERROR: The transcript list provided is empty. "
            "Please provide a non-empty file with transcript IDs."
        )

    # Transcript list validation: check if all user-provided IDs are in the GTF
    gtf_tx_set = set(df_txlengths.transcript_id)
    tx_set = set(transcript_ids)
    logging.info(f"Number of transcript IDs in {transcript} : {len(tx_set)}")

    # Find matching and missing IDs
    matching_ids = tx_set & gtf_tx_set
    missing_ids = tx_set - gtf_tx_set

    # Check if none or some of the transcript IDs match
    if len(matching_ids) == 0:
        logging.error("None of the provided transcript IDs are found in the GTF.")
        raise ValueError(
            "ERROR: None of the provided transcript IDs are found in the GTF. "
            "Please ensure the transcript list matches the GTF, "
            "or omit --representative_transcript to let the pipeline automatically select representative transcripts."
        )
    elif len(missing_ids) > 0:
        logging.error(f"{len(missing_ids)} transcript IDs not found in the GTF: {sorted(missing_ids)[:10]}{' ...' if len(missing_ids) > 10 else ''}")
        raise ValueError(
            "Some user-provided transcript IDs are missing from the GTF. "
            "Please ensure the transcript list matches the GTF, "
            "or omit --representative_transcript to let the pipeline automatically select representative transcripts."
        )
    else:
        logging.info("All provided transcript IDs are found in the GTF.")

    # Filter features dataframe based on provided transcript IDs
    df_txlengths_filtered = df_txlengths.loc[df_txlengths.transcript_id.isin(transcript_ids)]

    return df_txlengths_filtered, transcript_ids


def select_longest_transcript(df_txlengths, df_gtf):
    """
    Selects the longest transcript per gene based on the sorting hierarchy:
    CDS length > exon length > unspliced length. Ties are broken alphabetically by transcript ID.

    Args:
    - df_txlengths: DataFrame containing transcript lengths (cds_length, exon_length, unspliced_length)
    - df_gtf: Full parsed GTF DataFrame, used to ensure the order of transcript IDs matches the original GTF

    Returns:
    - df_txlengths_filtered (pandas.DataFrame): A filtered DataFrame containing the longest transcript per gene.
    - transcript_ids_sorted (list): A list of transcript IDs sorted in the same order as in the original GTF.
    """
    # Sort by the lengths of CDS, exon, and unspliced regions in descending order,
    # and alphabetically by transcript ID in ascending order
    df_txlengths_sorted = df_txlengths.sort_values(
        by=['cds_length', 'exon_length', 'unspliced_length', 'transcript_id'],
        ascending=[False, False, False, True]
    )

    # For each gene, keep the longest (first) transcript after sorting
    df_txlengths_filtered = df_txlengths_sorted.drop_duplicates(subset='gene_id', keep='first')

    # Extract the unique transcript IDs from the filtered DataFrame
    # Ensure the order of the transcript IDs matches the original GTF order
    transcript_ids_sorted = df_gtf.loc[(df_gtf['Feature'] == 'transcript') & (df_gtf['transcript_id'].isin(df_txlengths_filtered.transcript_id.tolist())), 'transcript_id'] \
        .unique() \
        .tolist()

    return df_txlengths_filtered, transcript_ids_sorted


def main(process_name, gtf, transcript, output, skip_filter_gtf):
    # Dump version file
    dump_versions(process_name)

    # Logging
    log_file = f"{output}.log"
    logging.basicConfig(
        filename=log_file,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info(f"Command-line arguments: {vars(args)}")
    # Load the genome annotation file and calculate CDS and exon lengths per transcript
    df_gtf, df_txlengths, input_genes = compute_gtf_feature_lengths(gtf)
    df_txlengths.to_csv(f"{os.path.basename(gtf)}.lengths.tsv", index=None, sep='\t', quoting=csv.QUOTE_NONE)
    logging.info(f"Total genes in genome GTF: {len(input_genes)}")

    # If user specified transcripts are provided, use these to filter gtf; else find longest CDS / exon transcript.
    if transcript != "":
        df_txlengths_filtered, transcript_ids_sorted = load_and_validate_transcripts(df_txlengths, transcript)

    else:
        df_txlengths_filtered, transcript_ids_sorted = select_longest_transcript(df_txlengths, df_gtf)

    # Check if all filtered genes have exactly one matching transcript ID
    # Count number of transcripts per gene
    transcripts_per_gene = df_txlengths_filtered.groupby("gene_id").transcript_id.nunique()

    # Raise error if any gene has not exactly one transcript and print problematic genes
    multi_assigned_genes = transcripts_per_gene[transcripts_per_gene > 1].index.tolist()
    if not (transcripts_per_gene == 1).all():
        logging.error("Some genes do not have exactly one representative transcript. Offending genes (more than 1 transcript assigned): ", multi_assigned_genes)
        raise ValueError(
            "ERROR: Some genes do not have exactly one representative transcript. "
            "Please make sure you provide a single transcript ID per gene ID, "
            "or omit --representative_transcript to let the pipeline automatically select representative transcripts."
        )

    set_filtgenes = set(df_txlengths_filtered.gene_id)
    logging.info(f"Remaining genes after filtering by transcript IDs: {len(set_filtgenes)}")

    # Check that all genes from GTF have been assigned at least 1 transcript
    missing_genes = input_genes - set_filtgenes
    if len(missing_genes) > 0:
        missing_genes_list = sorted(missing_genes)
        top = missing_genes_list[:10]
        more = f" and {len(missing_genes_list) - 10} more" if len(missing_genes_list) > 10 else ""

        logging.error(f"Not all genes in the GTF have been assigned a representative transcript. Genes without trancripts, after filtering: {top}{more}")
        raise ValueError(
            "ERROR: Not all genes in the GTF have a representative transcript after filtering. "
            "Please make sure you provide a single transcript ID for every gene ID in the GTF, "
            "or omit --representative_transcript to let the pipeline automatically select representative transcripts."
        )

    # Save transcript IDs to file
    logging.info(f"Saving representative transcript IDs to {output}.txt")
    with open(output + ".txt", "w") as f:
        f.write("\n".join(map(str, transcript_ids_sorted)) + "\n")

    # Create a transcript fai file
    # Set transcript_id as index
    df_txlengths_filtered.set_index('transcript_id', inplace=True)
    # Sort in order of transcript_ids
    df_txlengths_filtered = df_txlengths_filtered.loc[transcript_ids_sorted]
    # Reset index
    df_txlengths_filtered.reset_index(inplace=True)

    # Save transcript id and length (sum of all exons) to .fai file, without header
    logging.info(f"Saving representative transcript fai to {output}.fai")
    df_txlengths_filtered[['transcript_id', 'exon_length']].to_csv(f"{output}.fai", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)

    # Create the transcript GTF file
    df_txlengths_filtered['src'] = 'src'
    df_txlengths_filtered['gene'] = 'gene'
    df_txlengths_filtered['start'] = 1
    df_txlengths_filtered['score'] = '.'
    df_txlengths_filtered['strand'] = '+'
    df_txlengths_filtered['frame'] = '.'
    df_txlengths_filtered['attributes'] = df_txlengths_filtered.apply(lambda row: f'gene_id:{row["gene_id"]}; transcript_id:{row["transcript_id"]}', axis=1)

    logging.info(f"Saving representative transcript GTF to {output}.gtf")
    # Order into gtf format
    df_txlengths_filtered[['transcript_id', 'src', 'gene', 'start', 'exon_length', 'score', 'strand', 'frame', 'attributes']].to_csv(f"{output}.gtf", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)

    # Filter the genome annotation if genome GTF filtering enabled
    logging.info(f"skip_filter_gtf: {skip_filter_gtf}")
    if skip_filter_gtf == "true":
        logging.info("Skipping GTF filtering as requested with --skip_filter_gtf.")
    else:
        filt_annot = df_gtf.loc[(df_gtf['transcript_id'].isin(transcript_ids_sorted)) | (df_gtf['Feature'] == 'gene')].copy()
        # Convert the filtered annotation to pyranges
        pr_gtf = pr.PyRanges(filt_annot[[c for c in filt_annot.columns if c!='length']])
        # Save the filtered annotation
        logging.info(f"Saving genome GTF filtered by representative transcripts to {output}_filtered.gtf")
        pr_gtf.to_gtf(f"{output}_filtered.gtf")

    logging.info("Completed.")

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--gtf", default="!{gtf}")
    parser.add_argument("--transcript", default="!{transcript}")
    parser.add_argument("--output", default="!{output}")
    parser.add_argument("--skip_filter_gtf", default="!{skip_filter_gtf}", help="Skip GTF filtering.")
    args = parser.parse_args()

    main(args.process_name, args.gtf, args.transcript, args.output, args.skip_filter_gtf)

