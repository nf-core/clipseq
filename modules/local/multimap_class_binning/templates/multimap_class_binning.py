#!/usr/bin/env python3
"""
Count unique reads per BAM that overlap features in the Telescope GTF,
summarised by class_id (repeat elements) or gene_type (ncRNA).

Rules:
 - Each read is counted at most once.
 - If a read's alignments overlap features in MORE THAN ONE category
   (class_id / gene_type), the read is discarded and logged.
 - Reads mapping to exactly one category are counted once for that category.

Uses pybedtools for C-speed interval intersections, then pandas.
"""

import pybedtools
import argparse
import pysam
import pandas as pd
import re
import os

# -- Arguments ----------------------------------------------------------------

parser = argparse.ArgumentParser(
    prog='Count repeat reads',
    description='Counts multi-mapped reads which have been assigned multiple times' \
    'to a single class of repeat element. Discards other reads.'
)

parser.add_argument('--gtf',
                    help = "Annotation for re-assigning multi-mapped reads, gtf format")
parser.add_argument('--bamdir',
                    help = "Directory where input coord-sorted bam files are stored")
parser.add_argument('-o', '--outdir',
                    help = "Directory for program to store output",
                    default = ".")
parser.add_argument('-i', '--input',
                    help = "Comma separated list of coord-sorted bam files",
                    nargs="+")
parser.add_argument('--classifiers',
                    help = "Comma separated list of strings which identify classifiers for bucketing from" \
                    "the 9th column of your GTF, e.g. 'gene_type, class'", default = "")
parser.add_argument('--col',
                    help= "Column which you want the features to be read from",
                    default = 9)
parser.add_argument('--name',
                    help = "The suffix for the name of the discarded and aligned reads files",
                    default = "")
parser.add_argument('--strand',
                    help = "Strandedness argument for bedtools intersect, either -s or -S",
                    default = "")

args = parser.parse_args()

# ── Configuration ────────────────────────────────────────────────────────────
GTF_PATH = args.gtf
BAM_DIR  = args.bamdir
BAM_FILES = args.input
OUTPUT_DIR = args.outdir
CLASSIFIERS = args.classifiers
COL = args.col
NAME = args.name
STRAND = args.strand

# ── Helpers ──────────────────────────────────────────────────────────────────
def extract_sample_name(bam_filename):
    """Sample name - everything before ."""
    return bam_filename.split(".")[0]

def parse_attr(attr_str, key):
    m = re.search(rf'{key}\s+"([^"]+)"', attr_str)
    return m.group(1) if m else None


# ── Ensure BAM index exists ─────────────────────────────────────────────────
def ensure_index(bam_path):
    """Create .bai if missing."""
    bai = bam_path + ".bai"
    if not os.path.exists(bai):
        print(f"  Indexing {os.path.basename(bam_path)} …")
        pysam.index(bam_path)

# ── Ensure BAM is coord-sorted ──────────────────────────────────────────────

def check_sorting(bam):
    sort_order = bam.header.get("HD", {}).get("SO", "unknown")
    if sort_order != "coordinate":
        raise ValueError(
            f"'{bam.filename.decode()}' is '{sort_order}'-sorted, "
            "bam files must be coordinate-sorted."
        )

# ── GTF → sorted BED4 with category ─────────────────────────────────────────
def gtf_to_category_bed(gtf_path):
    """Return a sorted pybedtools BedTool of BED4: chrom start end category.

    category = class_id (TEs) or gene_type (ncRNAs).
    Writes to a temp file to avoid holding 4.8 M intervals in memory.
    """
    print("Building feature BED from GTF …")
    bed_path = os.path.join(OUTPUT_DIR, "_features_unsorted.bed")
    n = 0
    skipped = 0
    cats_list = []
    sources_list = []
    with open(gtf_path) as fin, open(bed_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            chrom = cols[0].strip()
            strand = cols[6].strip()
            if not chrom:               # skip lines with empty chromosome
                skipped += 1
                continue
            start = int(cols[3]) - 1    # GTF 1-based → BED 0-based
            end   = int(cols[4])
            #print(cols[8])
            source = cols[1]
            attrs = cols[int(COL)-1]
            cat = ""
            classifiers = CLASSIFIERS.split(", ")

            if COL == 9:
                for i in classifiers:
                    if cat:
                        continue
                    cat = parse_attr(attrs, i)

                # Replace any whitespace in category with underscore
                cat = cat.replace(" ", "_")
                cats_list.append(cat)
                cat = source + ": " + cat
                fout.write(f"{chrom}\t{start}\t{end}\t{cat}\t0\t{strand}\n")
                n += 1
            else:
                cat = cols[int(COL)-1]
                cats_list.append(cat)
                cat = source + ": " + cat
                cat = cat.replace(" ", "_")
                fout.write(f"{chrom}\t{start}\t{end}\t{cat}\t0\t{strand}\n")
                n += 1
            sources_list.append(source)

    # Warning if multiple sources for one category in the GTF
    cats_source_df = pd.DataFrame({"cats": cats_list, "source": sources_list})
    num_sources = cats_source_df.groupby("cats")["source"].transform("nunique")
    print(num_sources.head())
    duplicates = cats_source_df[num_sources > 1]
    print(duplicates.head())
    if (len(duplicates) > 0):
        print("WARNING: Multiple sources found for one category in your GTF, this could "
        "erroneously lead to reads being discarded.")
        print(duplicates.drop_duplicates(subset=["cats", "source"]))

    if skipped:
        print(f"  Skipped {skipped:,} lines with empty chromosome")

    bt = pybedtools.BedTool(bed_path).sort()
    print(f"  {n:,} features loaded")
    return bt


# ── Process one BAM ─────────────────────────────────────────────────────────
def process_bam(bam_path, feature_bed_path):
    """Intersect BAM with feature BED entirely through shell pipeline.

    Phase 1: bedtools bamtobed | sort | bedtools intersect -sorted
             → produces overlap lines in C, fast
    Phase 2: awk collects unique (read, category) pairs per read,
             prints reads with exactly 1 category (counted),
             and dumps multi-category reads to stderr (logged).

    Returns (counts_dict, n_ambiguous).
    """
    import subprocess

    label = os.path.splitext(os.path.basename(bam_path))[0]
    print(f"Processing: {label}")
    ensure_index(bam_path)

    # Build genome file from BAM header for -sorted
    genome_tmp = os.path.join(OUTPUT_DIR, "_genome_tmp.txt")
    import pysam as _pysam
    with _pysam.AlignmentFile(bam_path, "rb") as bam:
        check_sorting(bam)

        chroms = [(sq["SN"], sq["LN"]) for sq in bam.header["SQ"]]
    chroms.sort(key=lambda x: x[0])
    with open(genome_tmp, "w") as f:
        for c, l in chroms:
            f.write(f"{c}\t{l}\n")

    # Entire pipeline in C / awk — no Python iteration
    # intersect output: BED6(read) + BED4(feature) + overlap
    #   $4 = read_name, $10 = category
    #
    # The awk script:
    #   1) For each line, records the read's category set.
    #   2) At END, counts reads with exactly 1 category,
    #      and prints multi-category reads to stderr.
    awk_script = r"""
    {
        if ($10 == ".") { no_feature[$4] = 1; next }
        rn = $4
        cat = $10
        if (!(rn in cats)) {
            cats[rn] = cat
        } else {
            # append only if new category
            n = split(cats[rn], arr, SUBSEP)
            found = 0
            for (i = 1; i <= n; i++) {
                if (arr[i] == cat) { found = 1; break }
            }
            if (!found) cats[rn] = cats[rn] SUBSEP cat
        }
    }
    END {
        for (rn in no_feature) {
            if (!(rn in cats))
            printf "NO_FEATURE\t%s\n", rn > "/dev/stderr"
            }
        for (rn in cats) {
            n = split(cats[rn], arr, SUBSEP)
            if (n == 1) {
                count[arr[1]]++
            } else {
                # ambiguous → stderr
                printf "%s\t%s\n", rn, cats[rn] > "/dev/stderr"
                ambig++
            }
        }
        for (c in count) print count[c] "\t" c
        printf "AMBIG\t%d\n", ambig+0 > "/dev/stderr"
    }
    """

    with open(os.path.join(OUTPUT_DIR, "_awk_dedup.awk"), "w") as f:
        f.write(awk_script)

    # path for intersected bed
    intersect_path = os.path.join(OUTPUT_DIR, f"{label}_intersected.bed")

    cmd = (
        f"bedtools bamtobed -i '{bam_path}' "
        f"| sort -k1,1 -k2,2n "
        f"| bedtools intersect -a stdin -b '{feature_bed_path}' "
        f"  -loj -sorted {STRAND} -g '{genome_tmp}' "
        f"| tee '{intersect_path}' "
        f"| awk -F'\\t' -f '{OUTPUT_DIR}/_awk_dedup.awk'"
    )

    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, executable="/bin/bash",
    )
    if result.returncode != 0 and not result.stdout:
        print(f"  ERROR: {result.stderr.strip()}")
        return {}, 0

    # Parse counts (stdout)
    counts = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t", 1)
        if len(parts) == 2 and parts[0].isdigit():
            counts[parts[1]] = int(parts[0])

    # Parse ambiguous count and reads (stderr)
    ambig_lines = result.stderr.strip().split("\n") if result.stderr.strip() else []
    n_ambiguous = 0
    ambig_reads = []
    for line in ambig_lines:
        if line.startswith("AMBIG\t"):
            n_ambiguous = int(line.split("\t")[1])
        elif line.startswith("NO_FEATURE\t"):
            rn = line.split("\t")[1]
            ambig_reads.append((rn, "no_feature"))
        else:
            parts = line.split("\t", 1)
            if len(parts) == 2:
                ambig_reads.append((parts[0], parts[1]))

    # -------------------------------------------------------------------------------------------

    total = sum(counts.values())
    print(f"  {total:,} reads kept, {n_ambiguous:,} discarded (multi-category)")

    return counts, ambig_reads


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    feature_bt = gtf_to_category_bed(GTF_PATH)
    # Get path of the sorted BED file for shell pipeline
    feature_bed_path = feature_bt.fn

    all_rows = []
    all_ambiguous = {}

    for bam_file in BAM_FILES:
        sample_name  = extract_sample_name(bam_file)
        bam_path = os.path.join(BAM_DIR, bam_file)
        counts, ambig = process_bam(bam_path, feature_bed_path)

        label = os.path.splitext(os.path.basename(bam_path))[0]

        for cat, n in counts.items():
            all_rows.append({"sample": sample_name, "category": cat, "read_count": n})
        all_ambiguous[sample_name] = ambig

    df = pd.DataFrame(all_rows)
    df_wide = df.pivot(index="sample", columns="category", values="read_count").reset_index()

    # Filter so only categories with > 5% reads are on multiQC graph
    df_5percent = df.copy()
    df_5percent['read_percent']= (df_5percent['read_count']/df_5percent.groupby('sample')['read_count'].transform('sum'))*100
    df_5percent = df_5percent[df_5percent['read_percent'] > 5]
    del df_5percent['read_percent']

    df_5perc_wide = df_5percent.pivot(index="sample", columns="category", values="read_count").reset_index()


    # ── Save CSV summary ─────────────────────────────────────────────────

    csv_path = os.path.join(OUTPUT_DIR, f"repeat_ncRNA_read_counts_{NAME}.csv")
    df_wide.to_csv(csv_path, index=False)
    print(f"\nSummary saved to: {csv_path}")

    csv_path_5p = os.path.join(OUTPUT_DIR, f"repeat_ncRNA_read_counts5_perc_{NAME}.csv")
    df_5perc_wide.to_csv(csv_path_5p, index=False)

    pivot = df.pivot_table(
        index="category", columns="sample",
        values="read_count", fill_value=0,
    )
    print(pivot.to_string())

    # ── Log discarded ambiguous reads ────────────────────────────────────
    log_path = os.path.join(OUTPUT_DIR, f"discarded_ambiguous_reads_{NAME}.tsv")
    with open(log_path, "w") as f:
        f.write("sample\tread_name\tcategories\n")
        for sample, ambig_list in all_ambiguous.items():
            for rname, cats_str in ambig_list:
                # cats_str uses SUBSEP (0x1c) as separator from awk
                cats = cats_str.replace("\x1c", ",")
                f.write(f"{sample}\t{rname}\t{cats}\n")
    n_total_ambig = sum(len(v) for v in all_ambiguous.values())
    print(f"\n{n_total_ambig:,} ambiguous reads logged to: {log_path}")


if __name__ == "__main__":
    main()
