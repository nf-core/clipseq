#!/usr/bin/env python

# Script to get tables for custom MultQC plotting
# A. M. Chakrabarti
# 16th December 2020

import os
import sys
import re
import pybedtools
import pandas as pd
import numpy as np

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
dedup_df.to_csv('dedup_metrics.tsv', sep = '\t', index = False)

# Subset for MultiQC plots
dedup_df.loc[:, ['exp', 'input_reads', 'output_reads']].to_csv('dedup_reads.tsv', sep = '\t', index = False)
dedup_df.loc[:, ['exp', 'mean_umis']].to_csv('dedup_mean_umis.tsv', sep = '\t', index = False)
dedup_df.loc[:, ['exp', 'ratio']].to_csv('dedup_ratio.tsv', sep = '\t', index = False)

# ==========
# Crosslinks
# ==========

def read_bed(filename):
    df = pd.read_table(filename, header = None, names = ['chr', 'start', 'end', 'name', 'score', 'strand'], dtype = {'chr':str, 'start':int, 'end':int, 'name':str, 'score':float, 'strand':str})
    return df

# First get xlink bed files
xlinks_files = sorted(['xlinks/' + f for f in os.listdir('xlinks') if f.endswith('.xl.bed.gz')])

xlinks = dict((key, []) for key in ['exp', 'total_xlinks', 'total_xlinksites', 'ratio'])

for xlinks_file in xlinks_files:

    xlinks_df = read_bed(xlinks_file)

    exp = re.sub('.xl.bed.gz', '', os.path.basename(xlinks_file))
    total_xlinks = xlinks_df['score'].sum()
    total_xlinksites = xlinks_df.shape[0]
    ratio = total_xlinks/total_xlinksites

    xlinks['exp'].append(exp)
    xlinks['total_xlinks'].append(total_xlinks)
    xlinks['total_xlinksites'].append(total_xlinksites)
    xlinks['ratio'].append(round(total_xlinks/total_xlinksites, 2))

xlinks_metrics_df = pd.DataFrame(xlinks)
xlinks_metrics_df.to_csv('xlinks_metrics.tsv', sep = '\t', index = False)

# Subset for MultiQC plots
xlinks_metrics_df.loc[:, ['exp', 'total_xlinks', 'total_xlinksites']].to_csv('xlinks_counts.tsv', sep = '\t', index = False)
xlinks_metrics_df.loc[:, ['exp', 'ratio']].to_csv('xlinks_ratio.tsv', sep = '\t', index = False)

# ==========
# Peaks
# ==========

peakcallers = sys.argv[1:]

if not peakcallers:
    sys.exit()

pc_check = all(pc in ['icount', 'paraclu', 'pureclip', 'piranha'] for pc in peakcallers)
if pc_check == False:
    sys.exit('ERROR: peakcallers should be any of: icount, paraclu, pureclip, piranha')


def get_peaks_metrics(peakcaller):

    peak_files = sorted([peakcaller + '/' + f for f in os.listdir(peakcaller) if f.endswith('.peaks.bed.gz')])

    peaks = dict((key, []) for key in ['exp', 'peakcaller', 'total_xlinks', 'total_xlinksites', 'total_peaks', 'median_peak_width', 'mean_peak_width', 'xlinks_in_peaks', 'xlinks_in_peaks_percent', 'xlinksites_in_peaks', 'xlinksites_in_peaks_percent', 'peaks_xlinksite_coverage_percent'])

    for peak_file in peak_files:

        peaks_df = read_bed(peak_file)

        if peakcaller == 'icount':
            exp = re.sub('.[0-9]*nt_[0-9]*nt.peaks.bed.gz', '', os.path.basename(peak_file))
        elif peakcaller == 'paraclu':
            exp = re.sub('.[0-9]*_[0-9]*nt_[0-9]*.peaks.bed.gz', '', os.path.basename(peak_file))
        elif peakcaller == 'pureclip':
            exp = re.sub('.[0-9]*nt.peaks.bed.gz', '', os.path.basename(peak_file))
        elif peakcaller == 'piranha':
            exp = re.sub('.[0-9]*nt_[0-9]*nt.peaks.bed.gz', '', os.path.basename(peak_file))


        xlinks_df = read_bed('xlinks/' + exp + '.xl.bed.gz')
        expanded_xlinks_df = xlinks_df.loc[xlinks_df.index.repeat(xlinks_df.score)].reset_index(drop = True)

        # Get metrics
        total_xlinks = xlinks_df['score'].sum()
        total_xlinksites = xlinks_df.shape[0]
        total_peaks = peaks_df.shape[0]

        total_peak_width = sum((peaks_df['end'] - peaks_df['start']).tolist())
        mean_peak_width = round(np.mean((peaks_df['end'] - peaks_df['start']).tolist()), 2)
        median_peak_width = np.median((peaks_df['end'] - peaks_df['start']).tolist())

        xlinks_bed = pybedtools.BedTool.from_dataframe(xlinks_df)
        expanded_xlinks_bed = pybedtools.BedTool.from_dataframe(expanded_xlinks_df)
        peaks_bed = pybedtools.BedTool.from_dataframe(peaks_df)

        xlinks_in_peaks = peaks_bed.intersect(expanded_xlinks_bed, s = True, c = True)
        xlinks_in_peaks = sum([int(c[-1]) for c in xlinks_in_peaks])

        xlinksites_in_peaks = peaks_bed.intersect(xlinks_bed, s = True, c = True)
        xlinksites_in_peaks = sum([int(c[-1]) for c in xlinksites_in_peaks])

        peaks_xlinksite_coverage_percent = round((xlinksites_in_peaks/total_peak_width) * 100, 2)

        # Write out
        peaks['exp'].append(exp)
        peaks['peakcaller'].append(peakcaller)
        peaks['total_xlinks'].append(total_xlinks)
        peaks['total_xlinksites'].append(total_xlinksites)
        peaks['total_peaks'].append(total_peaks)
        peaks['median_peak_width'].append(median_peak_width)
        peaks['mean_peak_width'].append(mean_peak_width)
        peaks['xlinks_in_peaks'].append(xlinks_in_peaks)
        peaks['xlinks_in_peaks_percent'].append(round(xlinks_in_peaks/total_xlinks * 100, 2))
        peaks['xlinksites_in_peaks'].append(xlinksites_in_peaks)
        peaks['xlinksites_in_peaks_percent'].append(round(xlinksites_in_peaks/total_xlinksites * 100, 2))
        peaks['peaks_xlinksite_coverage_percent'].append(peaks_xlinksite_coverage_percent)
    
    peaks_metrics_df = pd.DataFrame(peaks)
    return(peaks_metrics_df)

peaks_metrics = [get_peaks_metrics(pc) for pc in peakcallers]
peaks_metrics_df = pd.concat(peaks_metrics)
peaks_metrics_df.to_csv('peaks_metrics.tsv', sep = '\t', index = False)

# Subset for MultiQC plots
total_peaks_df = peaks_metrics_df.pivot_table(index = 'exp', columns = 'peakcaller', values = 'total_peaks')
total_peaks_df.to_csv('total_peaks.tsv', sep = '\t', index = True)

xlinks_in_peaks_percent_df = peaks_metrics_df.pivot_table(index = 'exp', columns = 'peakcaller', values = 'xlinks_in_peaks_percent')
xlinks_in_peaks_percent_df.to_csv('xlinks_in_peaks.tsv', sep = '\t', index = True)

xlinksites_in_peaks_percent_df = peaks_metrics_df.pivot_table(index = 'exp', columns = 'peakcaller', values = 'xlinksites_in_peaks_percent')
xlinksites_in_peaks_percent_df.to_csv('xlinksites_in_peaks.tsv', sep = '\t', index = True)

peaks_xlinksite_coverage_percent_df = peaks_metrics_df.pivot_table(index = 'exp', columns = 'peakcaller', values = 'peaks_xlinksite_coverage_percent')
peaks_xlinksite_coverage_percent_df.to_csv('peaks_xlinksite_coverage.tsv', sep = '\t', index = True)