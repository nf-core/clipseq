# -----------------------------
# Module make binding sites
# -----------------------------
options(warn = -1)

# libraries
################
 # Suppress messages
suppressMessages( library(BindingSiteFinder))
suppressMessages(  library(GenomicRanges))
suppressMessages(  library(rtracklayer))
suppressMessages( library(tidyverse))
suppressMessages( library(optparse))

# print BSF version
##################
cat("BindingSiteFinder version: ", packageVersion("BindingSiteFinder"), "\n")


# Input
################

# Specify the input arguments
option_list <- list(
  make_option(c("-b", "--bw_files_folder"), type = "character"),
  make_option(c("-p", "--peaks"), type = "character"),
  make_option(c("-g", "--anno_genes"), type = "character"),
  make_option(c("-r", "--anno_regions"), type = "character"),
  # output paths
  make_option(c("-o", "--output_path"), type = "character", default = "."),

  # optional parameters to customize binding sites
  #--------------------------------
  # the default values are the ones specified in the Bioconductor package BindingSiteFinder
  # https://www.bioconductor.org/packages/release/bioc/manuals/BindingSiteFinder/man/BindingSiteFinder.pdf
  # the default values as specified in BindingSiteFinder 2.4.0 are  given in the comment after the parameter name
  # optional parameters for binding site defnition
  make_option(c( "--peak_score_global_cuttoff"), type = "numeric"), # default 0.01
  make_option(c( "--bsWidth"), type = "numeric"), # default automatic estimation
  make_option(c( "--peak_score_genewise_cuttoff"), type = "numeric"), # default automatic estimation
  make_option(c( "--minWidth"), type = "numeric"), # default 2
  make_option(c( "--minCrosslinks"), type = "numeric"), # default 2
  make_option(c( "--minCLSites"), type = "numeric"), # default 1
  make_option(c( "--maxBsWidth"), type = "numeric"), # default 13
  # optional parameters for reproducibility
  # make_option(c( "--reproducibility_cutoff"), type = "numeric"),
  # make_option(c( "--reproducibility_nReps"), type = "numeric"),
  # optional parameters for gene and region assignment
  make_option(c( "--method_gene_overlaps"), type = "character"), # default "frequency"
  make_option(c( "--rule_gene_overlaps"), type = "character"), # default NULL (only needed for --method_gene_overlaps "hierarchy")
  make_option(c( "--method_region_overlaps"), type = "character"), # default "frequency"
  make_option(c( "--rule_region_overlaps"), type = "character"), # default NULL (only needed for --method_region_overlaps "hierarchy")
  # optional parameters to fit non-standard genomes
  make_option(c( "--match_geneID"), type = "character"), # default "gene_id"
  make_option(c( "--match_geneName"), type = "character"), # default "gene_name"
  make_option(c( "--match_geneType"), type = "character"), # default "gene_type"
  # optional parameters to fit peakcaller scores
  make_option(c( "--match_score"), type = "numeric"), # default "score"
  make_option(c( "--match_score_option"), type = "character") # default "max"
)

# Parse arguments
parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Default parameters (NULL means use function defaults)
params.input.output <- list(
  bw_files_folder = args$bw_files_folder,
  peaks = args$peaks,
  anno_genes = args$anno_genes,
  anno_regions = args$anno_regions,
  output_path = args$output_path)

params.pureClipGlobalFilter <- list(
  cuttoff = args$peak_score_global_cuttoff)

params.estimateBsWidth <- list(
  est.minWidth = args$minWidth,
  est.maxBsWidth = args$maxBsWidth)

params.pureClipGeneWiseFilter <- list(
  cutoff = args$peak_score_genewise_cuttoff,
  match.score = args$match_score,
  match.geneID = args$match_geneID,
  overlaps = args$method_gene_overlaps)

params.makeBindingSites <- list(
  minWidth = args$minWidth,
  minCrosslinks = args$minCrosslinks,
  minClSites = args$minCLSites)

params.reproducibilityFilter <- list(
  cutoff = args$reproducibility_cutoff,
  nReps = args$reproducibility_nReps)

params.assignToGenes <- list(
  overlaps = args$method_gene_overlaps,
  overlaps.rule = args$rule_gene_overlaps,
  match.geneID = args$match_geneID,
  match.geneName = args$match_geneName,
  match.geneType = args$match_geneType
  )
params.assignToTranscriptRegions <- list(
  overlaps = args$method_region_overlaps,
  overlaps.rule = args$rule_region_overlaps)

params.annotateWithScore <- list(
  match.score = args$match_score,
  match.option = args$match_score_option
)

# Remove NULL values so function defaults apply
params.pureClipGlobalFilter <- params.pureClipGlobalFilter[!sapply(params.pureClipGlobalFilter, is.null)]
params.estimateBsWidth <- params.estimateBsWidth[!sapply(params.estimateBsWidth, is.null)]
params.pureClipGeneWiseFilter <- params.pureClipGeneWiseFilter[!sapply(params.pureClipGeneWiseFilter, is.null)]
params.makeBindingSites <- params.makeBindingSites[!sapply(params.makeBindingSites, is.null)]
params.reproducibilityFilter <- params.reproducibilityFilter[!sapply(params.reproducibilityFilter, is.null)]
params.assignToGenes <- params.assignToGenes[!sapply(params.assignToGenes, is.null)]
params.assignToTranscriptRegions <- params.assignToTranscriptRegions[!sapply(params.assignToTranscriptRegions, is.null)]
params.annotateWithScore <- params.annotateWithScore[!sapply(params.annotateWithScore, is.null)]

########################
# BindingSiteFinder 
#######################
# crosslinks
bw_files_names <- list.files(params.input.output$bw_files_folder)

clipFilesP <- list.files(params.input.output$bw_files_folder, pattern = "plus.bw$", full.names = TRUE)
clipFilesM <- list.files(params.input.output$bw_files_folder, pattern = "minus.bw$", full.names = TRUE)


# annotation
gns <- readRDS(params.input.output$anno_genes)
regions <- readRDS(params.input.output$anno_regions)


# Peaks from pureclip
peaks  = rtracklayer::import(con = params.input.output$peaks, format = "BED", extraCols=c("additionalScores" = "character"))
peaks$additionalScores = NULL
peaks$name = NULL



# Prepare meta data
meta = data.frame(
  id = c(1:length(clipFilesP)),
  condition = factor(rep("all", length(clipFilesP))), # add option for multiple groups from sample file
  clPlus = clipFilesP, 
  clMinus = clipFilesM)


# run BindingSiteFinder
#######################

cat("############################# \n Running BindingSiteFinder \n############################# \n")
bds = BSFDataSetFromBigWig(ranges = peaks, meta = meta, silent =T)

cat("\nGlobal filter on peak sites \n \n")
bds = do.call(pureClipGlobalFilter, 
              c(list(bds), 
                params.pureClipGlobalFilter)) # param cutoff
cat("\nEstimate binding site width \n \n")
bds = do.call(estimateBsWidth, 
              c(list(bds, anno.genes = gns), 
                params.estimateBsWidth)) # optional param: bsWidth

cat("\nGenewise filter on peak sites \n \n")
bds = do.call(pureClipGeneWiseFilter, 
              c(list(bds, anno.genes = gns),
                params.pureClipGeneWiseFilter)) # param cutoff, overlaps, match score, match geneID

cat("\nMake binding sites \n \n")
bds = do.call(makeBindingSites, 
              c(list(bds),
                params.makeBindingSites)) # params minWidth, minCrosslinks, minCLSites


# bds = do.call(reproducibilityFilter, c(list(bds), params.reproducibilityFilter)) # params cutoff, nReps
cat("\nAssign binding sites to genes \n \n")
bds = do.call(assignToGenes, c(list(bds, anno.genes = gns),
                               params.assignToGenes)) # params overlaps, overlaps.rule, match.geneID, match.geneName, match.geneType

cat("\nAssign binding sites to transcript regions \n \n")
bds = do.call(assignToTranscriptRegions, c(list(bds,anno.transcriptRegionList = regions),
                                           params.assignToTranscriptRegions))# params overlaps, overlaps.rule, 

cat("\nAnotate binding scores \n \n")
bds = do.call(annotateWithScore, c(list(bds, peaks),params.annotateWithScore)) #match.score, match.option

cat("\nSaving bindings sites \n \n")
bs_gr = getRanges(bds)
names(bs_gr) <- 1:NROW(bs_gr)
bs_df <- as.data.frame(bs_gr)


# export outputs
########################

exportToBED(bds, con = paste0(params.input.output$output_path , "/myBindingSites.bed"))
saveRDS(bds, paste0(params.input.output$output_path ,"/bds_object.rds"))
saveRDS(bs_df, paste0(params.input.output$output_path ,"/bindingSites.rds"))
write.csv(bs_df, file = paste0(params.input.output$output_path,"/bindingSites.csv"), row.names = FALSE)


