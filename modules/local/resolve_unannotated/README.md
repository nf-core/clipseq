# RESOLVE_UNANNOTATED

resolve_unannotated.py annotates genome segments that are not annotated when performing the segmentation on the filtered GTF annotation using iCount-Mini.

**Motivation**

When performing iCount-Mini segmentation on GTF annotation, the "regions.gtf" file is produced, in which the whole genome is flattened into the following genomic regions: CDS, UTR3, UTR5, ncRNA, intron, intergenic.
If the segmentation is performed on filtered genomic annotation, some regions remain unannotated.
This happens, because the annotation on a "gene" level covers all transcripts related to this gene, but some transcripts are removed by GTF filtering, potentially resulting in "gene" level annotation extending beyond the boundaries of remaining transcripts. As a result, we get unannotated regions, which are covered by a gene-level entry, but not by transcript-level entries which are required for segmentation.

**How it works**

resolve_unannotated.py finds unannotated regions in the regions GTF resulting from segmenting the filtered annotation (regions.gtf) and annotates them according to the regions obtained with the segmentation of unfiltered GTF. This means that segmentation with iCount-Mini is performed twice whenever GTF filtering is enabled: once on the original unfiltered GTF and once on the filtered GTFs. Regions generated from the unfiltered GTF are used to impute missing regions into the filtered regions GTF.

The process requires that in the unfiltered regions GTF, the entire genome is segmented. If missing regions exist in unfiltered GTF, the process will raise an error.

**Outputs**

Resolved regions GTF file with missing regions imputed from the unfiltered regions GTF.
