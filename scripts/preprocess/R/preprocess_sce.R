#' Preprocess SingleCellExperiment

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)

#library(scrna.utils)
#library(scrna.sceutils)
#library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Preprocess SingleCellExperiment")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--mito_thres', type='double',
                    help="Mitochondrial threshold")
parser$add_argument('--ribo_thres', type='double',
                    help="Ribosomal threshold")
parser$add_argument('--nmads', type='integer',
                    help="Number of MADs to filter at.", default = 3)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for preprocessed SCE.")

filter_cells <- function (sce, max_mito = Inf, max_ribo = Inf, ...) 
{
    sce_filtered <- scater::filter(sce, !isOutlier(log10_total_features_by_counts, 
        ...) & !isOutlier(log10_total_counts, ...))
    if (max_mito < 100) {
        sce_filtered <- scater::filter(sce_filtered, pct_counts_mitochondrial < 
            max_mito)
    }
    if (max_ribo < 100) {
        sce_filtered <- scater::filter(sce_filtered, pct_counts_ribosomal < 
            max_ribo)
    }
    return(sce_filtered)
}

args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

sce_filtered <- filter_cells(sce, nmads = args$nmads, type = "lower", 
                             log = TRUE, max_mito = args$mito_thres, max_ribo = args$ribo_thres)
#sce_filtered=normalize(sce_filtered)

qclust <- quickCluster(sce_filtered, min.size = 30)
sce_filtered <- computeSumFactors(sce_filtered, clusters = qclust)

print('size factor')
sce_filtered$size_factor <- sizeFactors(sce_filtered)

sce_normalized <- normalize(sce_filtered)
#saveRDS(sce_normalized,'~/Desktop/sce_normalized.rds')

sce_normalized <- runPCA(sce_normalized, ntop = 1000, ncomponents = 50, exprs_values = "logcounts")
sce_normalized <- runTSNE(sce_normalized, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)
sce_normalized <- runUMAP(sce_normalized, use_dimred = "PCA", n_dimred = 50, ncomponents = 2)

# Write outputs
saveRDS(sce_normalized, file = args$outfname)

cat("Completed.\n")


