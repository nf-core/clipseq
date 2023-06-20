#!/usr/bin/env python

import os
import sys
import errno
import argparse
import platform


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
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def dump_versions(process_name):
    with open("versions.yml", "w") as out_f:
        out_f.write(process_name + ":\n")
        out_f.write("    python: " + platform.python_version() + "\n")


def check_samplesheet(process_name, file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    group,replicate,fastq_1,fastq_2
    WT_PE,1,WT_LIB1_REP1_1.fastq.gz,WT_LIB1_REP1_2.fastq.gz
    WT_PE,1,WT_LIB2_REP1_1.fastq.gz,WT_LIB2_REP1_2.fastq.gz
    WT_SE,2,WT_LIB1_REP2_1.fastq.gz
    """
    # Dump version file
    dump_versions(process_name)

    # Init
    num_fastq_list = []
    sample_names_list = []
    sample_run_dict = {}

    with open(file_in, "r") as fin:
        ## Check header
        MIN_COLS = 3
        HEADER = ["group", "replicate", "fastq_1", "fastq_2"]
        # HEADER_LEN = len(HEADER)
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        ACTUAL_HEADER_LEN = len(header)

        # if header[: len(HEADER)] != HEADER:
        #     print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
        #     sys.exit(1)

        ## Check sample entries
        line_no = 1
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            ## Check if its just a blank line so we dont error
            if line.strip() == "":
                continue

            ## Check valid number of columns per row
            if len(lspl) != ACTUAL_HEADER_LEN:
                print_error(
                    "Invalid number of columns (found {} should be {})! - line no. {}".format(
                        len(lspl), len(HEADER), line_no
                    ),
                    "Line",
                    line,
                )

            ## Check valid number of populated columns per row
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, replicate, fastq_1, fastq_2 = lspl[: len(HEADER)]
            if sample:
                if sample.find(" ") != -1:
                    print_error("Group entry contains spaces!", "Line", line)
            else:
                print_error("Group entry has not been specified!", "Line", line)

            if sample:
                if sample.find(".") != -1:
                    print_error("Group entry contains dots!", "Line", line)

            ## Check replicate entry is integer
            if not replicate.isdigit():
                print_error("Replicate id not an integer", "Line", line)
            replicate = int(replicate)
            if replicate <= 0:
                print_error("Replicate must be > 0", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )
            num_fastq = len([fastq for fastq in [fastq_1, fastq_2] if fastq])
            num_fastq_list.append(num_fastq)

            ## Auto-detect paired-end/single-end
            sample_info = []
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = [sample, str(replicate), "0"]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = [sample, str(replicate), "1"]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Collect additional sample data
            extra_data = lspl[len(HEADER) :]
            sample_info = sample_info + extra_data

            ## Add fastq 1and fastq2
            sample_info = sample_info + [fastq_1, fastq_2]

            ## Create sample mapping dictionary = {sample: {replicate : [ single_end, fastq_1, fastq_2 ]}}
            if sample not in sample_run_dict:
                sample_run_dict[sample] = {}
            if replicate not in sample_run_dict[sample]:
                sample_run_dict[sample][replicate] = [sample_info]
            else:
                if sample_info in sample_run_dict[sample][replicate]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_run_dict[sample][replicate].append(sample_info)

            ## Store unique sample names
            if sample not in sample_names_list:
                sample_names_list.append(sample)

            line_no = line_no + 1

    ## Check data is either paired-end/single-end and not both
    if min(num_fastq_list) != max(num_fastq_list):
        print_error("Mixture of paired-end and single-end reads!")

    ## Calculate output header
    ouput_header = ["id", "group", "replicate", "single_end"]
    extra_header = header[len(HEADER) :]
    ouput_header = ouput_header + extra_header + ["fastq_1", "fastq_2"]

    ## Write validated samplesheet with appropriate columns
    if len(sample_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(ouput_header) + "\n")

            for sample in sorted(sample_run_dict.keys()):
                ## Check that replicate ids are in format 1..<NUM_REPS>
                uniq_rep_ids = set(sample_run_dict[sample].keys())
                if len(uniq_rep_ids) != max(uniq_rep_ids):
                    print_error(
                        "Replicate ids must start with 1!",
                        "Group",
                        sample,
                    )
                for replicate in sorted(sample_run_dict[sample].keys()):
                    ## Write to file
                    for idx, sample_info in enumerate(sample_run_dict[sample][replicate]):
                        sample_id = "{}_R{}_T{}".format(sample, replicate, idx + 1)
                        fout.write(",".join([sample_id] + sample_info) + "\n")


if __name__ == "__main__":
    # Allows switching between nextflow templating and standalone python running using arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--process_name", default="!{process_name}")
    parser.add_argument("--samplesheet", default="!{samplesheet}")
    parser.add_argument("--output", default="!{output}")
    args = parser.parse_args()

    check_samplesheet(args.process_name, args.samplesheet, args.output)
