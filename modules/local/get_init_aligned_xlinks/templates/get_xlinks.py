import pysam
import argparse

# --- Arguments

parser = argparse.ArgumentParser(
    prog='Amend bams to produce xlink sites',
    description='Identifies cross-link sites from bam files without removing'
    'relevant flags.'
)

parser.add_argument('-i', '--input',
                    help = "A comma separated list of bam files to be edited.", \
                    nargs="+")

args = parser.parse_args()

# --- Configuration

BAM_FILES = args.input


# --- Helpers

def extract_sample_name(bam_filename):
    """Sample name - everything before ."""
    return bam_filename.split(".")[0]

# --- Script

for BAM in BAM_FILES:

    file = pysam.AlignmentFile(BAM, "rb")

    print("Editing file: ", BAM)
    sample_name = extract_sample_name(BAM)

    # wb to write to file, template copies header
    xlink_file = pysam.AlignmentFile(f"{sample_name}_init_xl_coord.bam", "wb", template = file)

    for read in file:
        if read.is_reverse:
            xlink_site = read.reference_end
        else:
            # clamp to avoid negative coords for reads starting at position 0
            xlink_site = max(read.reference_start - 1, 0)

        read.reference_start = xlink_site
        read.cigar = [(0,1)]
        read.query_sequence = "N"
        read.query_qualities = [30]
        xlink_file.write(read)

        print("xlink file ", f"{sample_name}_init_xl_coord.bam", " completed.")

    file.close()
    xlink_file.close()
