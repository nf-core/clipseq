#!/usr/bin/env python

import os
import sys
import gzip
import platform
from Bio import SeqIO
import Bio

input_fq = "!{reads}"
output_fq = "!{prefix}.umi.fastq"

with gzip.open(input_fq, mode = 'rt') as f_in:
    with open(output_fq, mode = 'w') as f_out:
        for record in SeqIO.parse(f_in, 'fastq'):
            header = record.id.split(":")
            if '_' not in header[-1]:
                rearranged = ":".join(header[1:]) + '_rbc:' + header[0]
                record.id = rearranged
                record.name = rearranged
                record.description = rearranged
            SeqIO.write(record, f_out, 'fastq')

os.system('pigz ' + output_fq)

with open("versions.yml", "w") as out_f:
    out_f.write("!{process_name}" + ":\n")
    out_f.write("    python: " + platform.python_version() + "\n")
    out_f.write("    biopython: " + Bio.__version__ + "\n")