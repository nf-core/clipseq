#!/usr/bin/env python

# Script to get tables for custom MultQC plotting
# A. M. Chakrabarti
# 16th December 2020

import os
import re
import pandas as pd

# ==========
# Mapping
# ==========

# First get Bowtie2 logs
bowtie_logs = sorted(['premap/' + f for f in os.listdir('premap') if f.endswith('.log')])

smrna = dict((key, []) for key in ['exp', 'input_reads', 'smrna_reads'])

for bowtie_log in bowtie_logs:

    with open(bowtie_log, 'r') as logfile:

        exp = re.sub('.premap.log', '', os.path.basename(bowtie_log))

        lines = logfile.readlines()
        
        input_reads = int(re.findall(r'\d+', lines[0])[0])
        output_reads = [i for i in lines if 'aligned 0 times' in i]
        output_reads = int(re.findall(r'\d+', output_reads[0])[0])

        smrna['exp'].append(exp)
        smrna['input_reads'].append(input_reads)
        smrna['smrna_reads'].append(input_reads - output_reads)

smrna_df = pd.DataFrame(smrna)

# Next get STAR logs 
star_logs = sorted(['mapped/' + f for f in os.listdir('mapped') if f.endswith('.Log.final.out')])

genome = dict((key, []) for key in ['exp', 'genome_reads', 'unmapped_reads'])

for star_log in star_logs:

    with open(star_log, 'r') as logfile:

        exp = re.sub('.Log.final.out', '', os.path.basename(star_log))

        lines = logfile.readlines()    

        input_reads = [i for i in lines if 'Number of input reads' in i]
        input_reads = int(re.findall(r'\d+', input_reads[0])[-1])  

        genome_reads = [i for i in lines if 'Uniquely mapped reads number' in i]
        genome_reads = int(re.findall(r'\d+', genome_reads[0])[-1])  

        unmapped_reads = input_reads - genome_reads

        genome['exp'].append(exp)
        genome['genome_reads'].append(genome_reads)
        genome['unmapped_reads'].append(input_reads - genome_reads)        

genome_df = pd.DataFrame(genome)

# Combine the two
mapping_df = pd.merge(smrna_df, genome_df, on = 'exp')
mapping_df.to_csv('mapping_metrics.tsv', sep = '\t', index = False)

# Subset for MultiQC plots
mapping_df.loc[:, ['exp', 'smrna_reads', 'genome_reads', 'unmapped_reads']].to_csv('mapping.tsv', sep = '\t', index = False)

# ==========
# Deduplication
# ==========

dedup_logs = sorted(['dedup/' + f for f in os.listdir('dedup') if f.endswith('.log')])
dedup = dict((key, []) for key in ['exp', 'input_reads', 'output_reads', 'mean_umis', 'ratio'])

for dedup_log in dedup_logs:

    with open(dedup_log, 'r') as logfile:

        exp = re.sub('.log', '', os.path.basename(dedup_log))

        lines = logfile.readlines()
        
        input_reads = [i for i in lines if 'INFO Reads: Input Reads:' in i]
        input_reads = int(re.findall(r'\d+', input_reads[0])[-1])

        output_reads = [i for i in lines if 'Number of reads out:' in i]
        output_reads = int(re.findall(r'\d+', output_reads[0])[-1])

        mean_umis = [i for i in lines if 'Mean number of unique UMIs per position:' in i]
        mean_umis = float(re.findall(r'\d+', mean_umis[0])[-1])

        dedup['exp'].append(exp)
        dedup['input_reads'].append(input_reads)
        dedup['output_reads'].append(output_reads)
        dedup['mean_umis'].append(mean_umis)
        dedup['ratio'].append(round(input_reads/output_reads, 2))

dedup_df = pd.DataFrame(dedup)
dedup_df.to_csv('dedup_metrics.csv', sep = '\t', index = False)

# Subset for MultiQC plots
dedup_df.loc[:, ['exp', 'input_reads', 'output_reads']].to_csv('dedup_reads.tsv', sep = '\t', index = False)
dedup_df.loc[:, ['exp', 'mean_umis']].to_csv('dedup_mean_umis.tsv', sep = '\t', index = False)
dedup_df.loc[:, ['exp', 'ratio']].to_csv('dedup_ratio.tsv', sep = '\t', index = False)