#!/usr/bin/env Rscript

# command line arguments
args <- commandArgs(trailingOnly = TRUE)

print(args)
annoFile = args[1]
out_gns = args[2]
out_regions = args[3]

# for local tests
# annoFile = "/Users/melinaklostermann/Documents/projects/anno/GENCODEv31-p12/gencode.v31.annotation.gtf"
# out_gns = "/Users/melinaklostermann/Documents/projects/nf-core-clipseq/devel_folder/outputs/gns.rds"
# out_regions = "/Users/melinaklostermann/Documents/projects/nf-core-clipseq/devel_folder/outputs/regions.rds"

# libraries
library(GenomicFeatures)




# Make annotation database from gff3 file
annoDb = GenomicFeatures::makeTxDbFromGFF(file = annoFile, format = "gtf")
annoInfo = rtracklayer::import(annoFile, format = "gtf")


# Get genes as GRanges
gns = genes(annoDb)
idx = match(gns$gene_id, annoInfo$gene_id)
meta = cbind(elementMetadata(gns),
             elementMetadata(annoInfo)[idx,])
meta = meta[,!duplicated(colnames(meta))]
elementMetadata(gns) = meta

saveRDS(gns, out_gns)


# Get regions as Granges
cdseq = cds(annoDb)
intrns = unlist(intronsByTranscript(annoDb))
utrs3 = unlist(threeUTRsByTranscript(annoDb))
utrs5 = unlist(fiveUTRsByTranscript(annoDb))
regions = GRangesList(CDS = cdseq, Intron = intrns, UTR3 = utrs3, UTR5 = utrs5)

saveRDS(regions, out_regions)
