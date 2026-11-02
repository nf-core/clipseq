#!/usr/bin/env python3
"""
Read one or more motif‐enrichment TSVs, z‐score the PEKA‐scores (keeping missing as NaN),
and produce:
  1) Cosine‐distance clustermap across samples (if >1 input)
  2) Top‐10 k‐mer z-score heatmap (clustermap if >1 input,
     simple heatmap if only one)
"""

import argparse
from pathlib import Path

import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_distances
import seaborn as sns
import matplotlib.pyplot as plt


def load_and_z(df_path, score_col='PEKA-score'):
    # load; keep missing scores as NaN
    df = pd.read_csv(df_path, sep='\t', index_col=0)
    print(df)
    # compute mean and std ignoring NaNs
    mean = df[score_col].mean(skipna=True)
    std = df[score_col].std(skipna=True)
    # assign z-score; NaNs propagate
    df['z'] = (df[score_col] - mean) / std
    return df['z']


def plot_cosine_clustermap(zmat, out_png):
    # compute sample×sample cosine‐distance, treating NaN as zero
    filled = zmat.T.fillna(0)
    dist = pd.DataFrame(
        cosine_distances(filled),
        index=zmat.columns,
        columns=zmat.columns
    )
    cg = sns.clustermap(
        dist,
        metric='euclidean',  # using precomputed cosine distance
        row_cluster=True,
        col_cluster=True,
        cmap='Blues',
        vmin=0, vmax=1,
        figsize=(8, 8)
    )
    cg.fig.suptitle('Cosine‐distances based on motif enrichment')
    cg.savefig(out_png)
    plt.close(cg.fig)


def plot_topk_heatmap(zmat, out_png, cluster=True, topk = 5):
    # select top k per sample, unique across samples
    top_kmers = set()
    for col in zmat.columns:
        # drop NaN before selecting
        s = zmat[col].dropna()
        top_kmers |= set(s.nlargest(topk).index)
    # Convert top k-mers to list
    top_kmers = sorted(top_kmers)
    sub = zmat.loc[top_kmers].fillna(0)
    if cluster and sub.shape[1] > 1:
        # For more than one sample plot clustermap
        cg = sns.clustermap(
            sub,
            cmap='Blues',
            figsize=(sub.shape[1]*0.6, sub.shape[0]*0.4),
            dendrogram_ratio=(.1, .1),
        )
        cg.fig.suptitle(f'Top {topk} k-mers per sample')
        cg.savefig(out_png)
        plt.close(cg.fig)

        # SIZE ADJUST, move cbar, perhaps make into a heatmap instead of clustermap - add dendrogram separately.Font adjustment, always position words parallel to x-axis
        # CHANGE INPUT TO PATH - script collects samples itself


    else:
        plt.figure(figsize=(6, max(4, len(top_kmers)*0.4)))
        sns.heatmap(sub, square=True, cmap='Blues')
        plt.title(f'Top {topk} k-mers')
        plt.ylabel('k-mer')
        plt.xlabel('sample')
        plt.tight_layout()
        plt.savefig(out_png)
        plt.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('tsv_string', type=str,
                   help='space-separated list of TSV file paths')
    p.add_argument('-o', '--outdir', type=Path, default=Path('.'),
                   help='output directory for PNGs')
    args = p.parse_args()

    # convert input string to Path objects
    tsv_paths = [Path(x) for x in args.tsv_string.split()]

    args.outdir.mkdir(exist_ok=True, parents=True)

    # load and z-score each file
    zseries = {}
    for f in tsv_paths:
        sample = f.stem.split('.genome_5mer_distribution')[0] # VALIDATE THAT THIS NAMING WILL ALWAYS MATCH
        zseries[sample] = load_and_z(f)

    zmat = pd.DataFrame(zseries)
    # print(zmat.head())

    # Get top n kmers from number of samples
    topn = max(3, round(30/len(zseries)))

    # 1) cosine-distance clustermap if more than 2 samples passed
    if len(zseries) > 1:
        plot_cosine_clustermap(
            zmat,
            args.outdir / 'cosine_distance_clustermap.png'
        )
    else:
        print("Single sample: skipping cosine-distance plot.")

    # 2) top-10 k-mers heatmap
    plot_topk_heatmap(
        zmat,
        args.outdir / 'top10_kmer_heatmap.png',
        cluster=(len(zseries) > 1),
        topk=topn
    )


if __name__ == '__main__':
    main()
