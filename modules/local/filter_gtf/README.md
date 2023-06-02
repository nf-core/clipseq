**Motivation**

Pre-filtering genomic annotation is crucial to ensure that accurate genome-level segmentation is produced by iCount segment. When iCount segmentation is applied to unfiltered annotation, this can result in introns of protein-coding genes to be annotated as ncRNA and in shifting of CDSs and UTRs to include low-confidence transcripts.

**Filtering rules**

1. No filtering is done on gene level.
2. Entries below gene-level tagged as "basic" are retained. The transcripts tagged as "basic" form part of a subset of representative transcripts for each gene. This subset prioritizes full-length protein-coding transcripts over partial or non-protein-coding transcripts within the same gene, and intends to highlight those transcripts that will be useful to the majority of users.
3. If a gene contains entries with transcript_support_level 1 or 2, only those transcripts are retained by filtering and other less-supported trascripts are discarded.
