#!/usr/bin/env Rscript

# Script to generate metrics and plots for nf-core/clipseq
# A. M. Chakrabarti
# 10th December 2020

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rtracklayer))

args <- commandArgs(trailingOnly = TRUE)
args <- c("icount", "paraclu", "pureclip", "piranha")

setwd("/Users/chakraa2/projects/clip_pipeline/results_all")

# ==========
# Functions
# ==========

load_peaks <- function(peakcaller) {
  
  if(peakcaller == "xlinks") {
    
    peaks <- list.files("xlinks", pattern = ".bed.gz$", full.names = TRUE)
    peaks.list <- lapply(peaks, import.bed)
    peaks.list <- lapply(peaks.list, keepStandardChromosomes, pruning.mode = "coarse")
    names(peaks.list) <- gsub(".xl.bed.gz", "", basename(peaks))
    
  } else if(peakcaller == "icount") {
    
    peaks <- list.files("icount", pattern = ".peaks.bed.gz$", full.names = TRUE)
    peaks.list <- lapply(peaks, import.bed)
    peaks.list <- lapply(peaks.list, keepStandardChromosomes, pruning.mode = "coarse")
    names(peaks.list) <- gsub(".[0-9]*nt_[0-9]*nt.peaks.bed.gz", "", basename(peaks))
    
  } else if(peakcaller == "paraclu") {
    
    peaks <- list.files("paraclu", pattern = ".peaks.bed.gz$", full.names = TRUE)
    peaks.list <- lapply(peaks, import.bed)
    peaks.list <- lapply(peaks.list, keepStandardChromosomes, pruning.mode = "coarse")
    names(peaks.list) <- gsub(".[0-9]*_[0-9]*nt_[0-9]*.peaks.bed.gz", "", basename(peaks))
    
  } else if(peakcaller == "pureclip") {
    
    peaks <- list.files("pureclip", pattern = ".peaks.bed.gz$", full.names = TRUE)
    peaks.list <- lapply(peaks, import.bed)
    peaks.list <- lapply(peaks.list, keepStandardChromosomes, pruning.mode = "coarse")
    names(peaks.list) <- gsub(".[0-9]*nt.peaks.bed.gz", "", basename(peaks))
    
  } else if(peakcaller == "piranha") {
    
    peaks <- list.files("piranha", pattern = ".peaks.bed.gz$", full.names = TRUE)
    peaks.list <- lapply(peaks, import.bed)
    peaks.list <- lapply(peaks.list, keepStandardChromosomes, pruning.mode = "coarse")
    names(peaks.list) <- gsub(".[0-9]*nt_[0-9]*nt.peaks.bed.gz", "", basename(peaks))
    
  }
  
  return(peaks.list)
  
}

get_metrics <- function(xlinks.list, peaks.list, peakcaller) {

  # xlinks.list = named list of xlinks GRanges 
  # peaks.list = named list of peaks GRanges (from peak caller output)
    
  metrics.list <- lapply(seq_along(peaks.list), function(i) {
    
    xlinks.gr <- xlinks.list[[i]]
    peaks.gr <- peaks.list[[i]]
    
    stopifnot(names(xlinks.list)[i] == names(peaks.list)[i])
    
    total_xlinks <- sum(xlinks.gr$score)
    total_xlinksites <- length(xlinks.gr)
    total_peaks <- length(peaks.gr)
    
    median_peak_width <- median(width(peaks.gr))
    mean_peak_width <- mean(width(peaks.gr))
    
    expanded.xlinks.gr <- rep(xlinks.gr, times = xlinks.gr$score)
    xlinks_in_peaks <- sum(countOverlaps(peaks.gr, expanded.xlinks.gr))
    xlinksites_in_peaks <- sum(countOverlaps(peaks.gr, xlinks.gr))
    
    peaks_xlinksite_coverage_percent <- (xlinksites_in_peaks/sum(width(peaks.gr))) * 100
    
    metrics.df <- data.frame(exp = names(peaks.list)[i],
                             peakcaller = peakcaller,
                             total_xlinks = total_xlinks,
                             total_xlinksites = total_xlinksites,
                             total_peaks = total_peaks,
                             median_peak_width = median_peak_width,
                             mean_peak_width = mean_peak_width,
                             xlinks_in_peaks = xlinks_in_peaks,
                             xlinks_in_peaks_percent = (xlinks_in_peaks/total_xlinks) * 100,
                             xlinksites_in_peaks = xlinksites_in_peaks,
                             xlinksites_in_peaks_percent = (xlinksites_in_peaks/total_xlinksites) * 100,
                             peaks_xlinksite_coverage_percent = peaks_xlinksite_coverage_percent)
    
    return(metrics.df)
    
  })
  
  # metrics.dt <- rbindlist(metrics.list)
  metrics.df <- do.call(rbind, metrics.list)
  
  return(metrics.df)
  
}

# ==========
# Run analysis
# ==========

# ==========
# Crosslinks
# ==========

xlinks.list <- load_peaks(peakcaller = "xlinks")

xlinks.metrics.df <- data.frame(exp = names(xlinks.list),
                                total_xlinks = unlist(lapply(xlinks.list, function(x) sum(x$score))),
                                total_xlinksites = unlist(lapply(xlinks.list, length)))

write.table(xlinks.metrics.df, "../plotting/xlinks_metrics.tsv", row.names = FALSE, sep = "\t")

xlinks.ratio.df <- with(xlinks.metrics.df, data.frame(exp = exp,
                                                      xlink_ratio = total_xlinks/total_xlinksites))

write.table(xlinks.ratio.df, "../plotting/xlinks_ratio.tsv", row.names = FALSE, sep = "\t")

# ==========
# Peaks
# ==========

peakcallers <- args
if(length(peakcallers) == 0) quit(save = "no") # Exit if no peakcallers run

metrics.list <- lapply(peakcallers, function(x) {
  
  peaks.list <- load_peaks(peakcaller = x)
  peaks.metrics.df <- get_metrics(xlinks.list = xlinks.list,
                                  peaks.list = peaks.list,
                                  peakcaller = x)
  
  return(peaks.metrics.df)
  
})

# metrics.dt <- rbindlist(metrics.list)
metrics.df <- do.call(rbind, metrics.list)

total_peaks.df <- dcast(metrics.df, exp ~ peakcaller, value.var = "total_peaks")
write.table(total_peaks.df, "../plotting/total_peaks.tsv", row.names = FALSE, sep = "\t")

xlinks_in_peaks_percent.df <- dcast(metrics.df, exp ~ peakcaller, value.var = "xlinks_in_peaks_percent")
write.table(xlinks_in_peaks_percent.df, "../plotting/xlinks_in_peaks.tsv", row.names = FALSE, sep = "\t")

xlinksites_in_peaks_percent.df <- dcast(metrics.df, exp ~ peakcaller, value.var = "xlinksites_in_peaks_percent")
write.table(xlinksites_in_peaks_percent.df, "../plotting/xlinksites_in_peaks.tsv", row.names = FALSE, sep = "\t")

peaks_xlinksite_coverage_percent.df <- dcast(metrics.df, exp ~ peakcaller, value.var = "peaks_xlinksite_coverage_percent")
write.table(peaks_xlinksite_coverage_percent.df, "../plotting/peaks_xlinksite_coverage.tsv", row.names = FALSE, sep = "\t")

# # Peak calling metrics
# 
# p.total_peaks <- ggplot(metrics.dt, aes(x = exp, y = total_peaks, fill = peakcaller)) +
#   geom_col(position = "dodge") +
#   scale_fill_tableau() +
#   # facet_grid(. ~ peakcaller) +
#   coord_flip() +
#   scale_y_continuous(labels = comma) +
#   labs(subtitle = "Total peaks",
#        x = "Experiment",
#        y = "Count",
#        fill = "Peak caller") +
#   theme_minimal_grid() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# p.xlinks_in_peaks <- ggplot(metrics.dt, aes(x = exp, y = xlinks_in_peaks_percent, fill = peakcaller)) +
#   geom_col(position = "dodge") +
#   scale_fill_tableau() +
#   # facet_grid(. ~ peakcaller) +
#   coord_flip() +
#   scale_y_continuous(labels = comma) +
#   labs(subtitle = "Crosslinks within peaks",
#        x = "",
#        y = "Fraction",
#        fill = "Peak caller") +
#   theme_minimal_grid() + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
# 
# p.xlinksites_in_peaks <- ggplot(metrics.dt, aes(x = exp, y = xlinksites_in_peaks_percent, fill = peakcaller)) +
#   geom_col(position = "dodge") +
#   scale_fill_tableau() +
#   # facet_grid(. ~ peakcaller) +
#   coord_flip() +
#   labs(subtitle = "Crosslink sites within peaks",
#        x = "",
#        y = "Fraction",
#        fill = "Peak caller") +
#   theme_minimal_grid() + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
# 
# peaks1.patch <- p.total_peaks + p.xlinks_in_peaks + p.xlinksites_in_peaks +
#   plot_layout(guides = "collect") + plot_annotation(title = "CLIP peak metrics")
# 
# ggsave(peaks1.patch, filename = "peak_metrics_1.pdf", width = 297, height = 210, units = "mm")
# 
# p.peaks_xlinksite_coverage <- ggplot(metrics.dt, aes(x = exp, y = peaks_xlinksite_coverage_percent, fill = peakcaller)) +
#   geom_col(position = "dodge") +
#   scale_fill_tableau() +
#   # facet_grid(. ~ peakcaller) +
#   coord_flip() +
#   labs(subtitle = "Coverage of peaks with crosslink sites",
#        x = "",
#        y = "Fraction",
#        fill = "Peak caller") +
#   theme_minimal_grid() + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
# 
# peaks2.patch <- p.total_peaks + p.xlinks_in_peaks + p.xlinksites_in_peaks +
#   plot_layout(guides = "collect") + plot_annotation(title = "CLIP peak metrics")
# 
# ggsave(peaks2.patch, filename = "peak_metrics_2.pdf", width = 297, height = 210, units = "mm")
# 
# # Peak widths
# 
# p.median_peak_width <- ggplot(metrics.dt, aes(x = exp, y = median_peak_width, fill = peakcaller)) +
#   geom_col(position = "dodge") +
#   scale_fill_tableau() +
#   # facet_grid(. ~ peakcaller) +
#   coord_flip() +
#   scale_y_continuous(labels = comma) +
#   labs(title = "",
#        x = "Experiment",
#        y = "Median peak width",
#        fill = "Peak caller") +
#   theme_minimal_grid() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# p.mean_peak_width <- ggplot(metrics.dt, aes(x = exp, y = mean_peak_width, fill = peakcaller)) +
#   geom_col(position = "dodge") +
#   scale_fill_tableau() +
#   # facet_grid(. ~ peakcaller) +
#   coord_flip() +
#   scale_y_continuous(labels = comma) +
#   labs(title = "",
#        x = "",
#        y = "Mean peak width",
#        fill = "Peak caller") +
#   theme_minimal_grid() + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
# 
# peakwidths.patch <- p.median_peak_width + p.mean_peak_width + 
#   plot_layout(guides = "collect") + plot_annotation(title = "CLIP peak widths")
# 
# ggsave(peakwidths.patch, filename = "peak_metrics_3.pdf", width = 297, height = 210, units = "mm")
# 
# # # peak widths
# # icount.peaks.width.dt <- data.table(exp = rep(names(icount.list), elementNROWS(icount.list)),
# #                                     widths = unlist(lapply(icount.list, width)))
# # 
# # ggplot(icount.peaks.width.dt, aes(x = exp, y = widths, colour = exp)) +
# #   geom_boxplot() + 
# #   coord_flip() + 
# #   labs(title = "iCount peak width",
# #        x = "Experiment",
# #        y = "Peak width") +
# #   theme_minimal_grid() + theme(legend.position = "none")