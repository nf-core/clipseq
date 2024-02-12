#!/usr/bin/env python

import os
import sys
import errno
import argparse
import platform
import pandas as pd

print("Starting Python script")

def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    raise AssertionError(error_str)


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")


def check_samplesheet(process_name, file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    group_name,replicate_name,input_name,fastq
    PHO92_rep1,PHO92,INPUT,ERR3988345-yeast-quarter.fastq.gz
    PHO92_rep2,PHO92,INPUT,ERR3988069-yeast-quarter.fastq.gz
    input,INPUT,,input.fastq.gz
    """
    # Dump version file
    dump_versions(process_name)

    # Init
    sample_names_list = []
    sample_run_dict = {}

    with open(file_in, "r") as fin:
        ## Check header
        HEADER = ["sample_name","group_name", "input_name", "fastq"]
        HEADER_LEN = len(HEADER)
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        ACTUAL_HEADER_LEN = len(header)

        # simply check if number of cols in header is correct
        if ACTUAL_HEADER_LEN != HEADER_LEN:
            print_error(
                f"Header formatting is incorrect, ensure you have four columns: {HEADER}"
            )

        # check that columns names are correct and in correct order
        for i in range(HEADER_LEN):
            if HEADER[i] != header[i]:
                print_error(
                    f"Column position '{i+1}' should be '{HEADER[i]}' instead of '{header[i]}', ensure your columns are named correctly and in the correct order"
                )

        ## Check sample entries
        line_no = 1
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            ## Check if its just a blank line so we dont error
            if line.strip() == "":
                continue

            ## Check that sample_name and fastq are populated, this is the minimum required
            if not lspl[0]: #sample_name
                print_error(
                    f"Sample name not provided - line no. {line_no}",
                    "Line",
                    line,
                )
            if not lspl[3]: #fastq
                print_error(
                    f"FastQ not provided - line no. {line_no}",
                    "Line",
                    line,
                )

            ## Check sample name and group name entries for common weird characters
            ## We don't need to check input_name because later on we will check that it matches a sample name
            sample_name, group_name, input_name, fastq = lspl[: len(HEADER)]
            weird_chars = [" ",".",":",";","/","\\"]
            for weirdo in weird_chars:
                if sample_name.find(weirdo) != -1:
                    print_error(
                        f"Sample name entry contains '{weirdo}', remove or replace with '_' - line no. {line_no}"
                    )
                if group_name:
                    if group_name.find(weirdo) != -1:
                        print_error(
                            f"Group name entry contains '{weirdo}', remove or replace with '_' - line no. {line_no}"
                        )

            ## Check FastQ file extension
            if fastq:
                if fastq.find(" ") != -1:
                    print_error("FastQ file name contains spaces!", "Line", line)
                if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                    print_error(
                        "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                        "Line",
                        line,
                        )
            
    #  Do the remaining checks with a pandas dataframe
    df = pd.read_csv(file_in)

    # 1. Check for duplicated rows
    duplicated_rows = df[df.duplicated(subset=['sample_name', 'group_name', 'input_name'])]
    if not duplicated_rows.empty:
        for index, row in duplicated_rows.iterrows():
            print(
                f"Friendly FYI - Row {index + 2}: Sample '{row['sample_name']}' with group '{row['group_name']}' and input '{row['input_name']}' is duplicated. These specific samples will be merged before mapping."
                )

    # 2. Check for rows where sample_name is duplicated but group_name and/or input_name are different
    duplicate_sample_names = df[df.duplicated(subset='sample_name', keep=False)].groupby('sample_name').apply(lambda x: x.nunique())
    for sample, unique_counts in duplicate_sample_names.iterrows():
        if unique_counts['group_name'] > 1 or unique_counts['input_name'] > 1:
            print_error(
                f"You have duplicate sample names '{sample}' where the group and/or input are different. If you want the fastqs to be merged before mapping the group/input columns must also match. If you do not want the samples to be merged they must have different names."
                )
            return

    # 3. Check if samples in the same group have different inputs
    group_inputs = df.groupby('group_name')['input_name'].nunique()
    for group, unique_inputs in group_inputs.items():
        if unique_inputs > 1:
            print_error(
                f"Samples in the same group '{group}', should have the same input."
                )
            return

    # 4. Check if input_name is in set of all sample_names or group_names
    unique_samples = set(df['sample_name'].unique())
    unique_groups = set(df['group_name'].unique())
    for index, row in df.iterrows():
        if row['input_name'] not in unique_samples and row['input_name'] not in unique_groups:
            print_error(
                f"Row {index + 2}: input_name '{row['input_name']}' must refer to either another sample_name or group_name in the samplesheet."
                )
            return


    ## Write validated samplesheet with appropriate columns
    if len(sample_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(HEADER) + "\n")
            for sample in sorted(sample_run_dict.keys()):
                fout.write(",".join([sample_id] + sample_info) + "\n")


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--samplesheet", default="!{samplesheet}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    check_samplesheet(args.process_name, args.samplesheet, args.output)
