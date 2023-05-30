#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# Regexes for new tools in process get_software_versions
regexes = {
    "nf-core/clipseq": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
    "FastQC": ["v_fastqc.txt", r"FastQC v(\S+)"],
    "MultiQC": ["v_multiqc.txt", r"multiqc, version (\S+)"],
    "Cutadapt": ["v_cutadapt.txt", r"(\S+)"],
    "Bowtie2": ["v_bowtie2.txt", r"version (\S+)"],
    "STAR": ["v_star.txt", r"STAR_(\S+)"],
    "Samtools": ["v_samtools.txt", r"samtools (\S+)"],
    "umi_tools": ["v_umi_tools.txt", r"version: (\S+)"],
    "bedtools": ["v_bedtools.txt", r"bedtools v(\S+)"],
    "preseq": ["v_preseq.txt", r"Version: (\S+)"],
    "subread": ["v_subread.txt", r"Subread-align v(\S+)"],
    "rseqc": ["v_rseqc.txt", r"bam2fq.py (\S+)"],
    "iCount": ["v_icount.txt", r"(\S+)"],
    "PureCLIP": ["v_pureclip.txt", r"pureclip version: (\S+)"],
    "Piranha": ["v_piranha.txt", r"Piranha Version (\S+)"],
    "Paraclu": ["v_paraclu.txt", r"(\S+)"],
    "Meme": ["v_meme.txt", r"(\S+)"],
    "Python": ["v_python.txt", r"Python (\S+)"],
    "Pygments": ["v_pygments.txt", r"Pygments version (\S+)"],
    "pigz": ["v_pigz.txt", r"pigz (\S+)"],
    "Perl": ["v_perl.txt", r"\(v(\S+)\)"]
}

results = OrderedDict()
results["nf-core/clipseq"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results["FastQC"] = '<span style="color:#999999;">N/A</span>'
results["MultiQC"] = '<span style="color:#999999;">N/A</span>'
results["Cutadapt"] = '<span style="color:#999999;">N/A</span>'
results["Bowtie2"] = '<span style="color:#999999;">N/A</span>'
results["STAR"] = '<span style="color:#999999;">N/A</span>'
results["Samtools"] = '<span style="color:#999999;">N/A</span>'
results["umi_tools"] = '<span style="color:#999999;">N/A</span>'
results["bedtools"] = '<span style="color:#999999;">N/A</span>'
results["preseq"] = '<span style="color:#999999;">N/A</span>'
results["subread"] = '<span style="color:#999999;">N/A</span>'
results["rseqc"] = '<span style="color:#999999;">N/A</span>'
results["iCount"] = '<span style="color:#999999;">N/A</span>'
results["PureCLIP"] = '<span style="color:#999999;">N/A</span>'
results["Piranha"] = '<span style="color:#999999;">N/A</span>'
results["Paraclu"] = '<span style="color:#999999;">N/A</span>'
results["Meme"] = '<span style="color:#999999;">N/A</span>'
results["Python"] = '<span style="color:#999999;">N/A</span>'
results["Pygments"] = '<span style="color:#999999;">N/A</span>'
results["pigz"] = '<span style="color:#999999;">N/A</span>'
results["Perl"] = '<span style="color:#999999;">N/A</span>'


# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del results[k]

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'nf-core/clipseq Software Versions'
section_href: 'https://github.com/nf-core/clipseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.csv", "w") as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
