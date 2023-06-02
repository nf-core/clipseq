# ResolveUnannotated

ResolveUnannotated.py annotates genome segments that are not annotated by iCount segmentation.

**Motivation**

When performing iCount segmentation on GTF annotation, the whole genome is split into the following genomic regions: CDS, UTR3, UTR5, ncRNA, intron, intergenic.
If the segmentation is performed on filtered genomic annotation, some regions remain unnanotated.
This is due to the fact that the annotation on a "gene" level covers all transcripts related to this gene, including transcripts with low transcript_support_level, which were removed by filtering. As a result, we get unannotated regions, which are covered by a gene, but lack transcripts that are required for segmentation.

**How it works**

ResolveUnannotated.py finds unannotated regions in the iCount genomic segment (regions.gtf), annotates them according to the corresponding gene, and adds them to the original segment under the feature "genic_other".

**Dependencies**:

```
python=3.7
pandas
plumbum
pybedtools

Versions of packages that were used for development and testing:
pandas=0.24.2
plumbum=1.6.8
pybedtools=0.8.0
```

**Usage**

```
python3 ResolveUnnanotated.py [-h] -s SEGMENTATION -a ANNOTATION -fai FASTA_INDEX -o OUTPUTDIR
required arguments:
  -s, --segmentation SEGMENTATION
                        iCount genome level segmentation in GTF format i.e. "regions.gtf.gz".
  -a, --annotation ANNOTATION
                        Annotation file from GENCODE or ENSEMBL in GTF format, which was as input for iCount segmentation.
  -fai, --fasta_index FASTA_INDEX
                       Fasta index file generated with samtools.
  -o, --outputdir OUTPUTDIR
                        Path to output folder.
```

**Outputs**

Gtf file with annotated "genic_other" regions.
