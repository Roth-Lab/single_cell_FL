# Differential expression accounting for malignant status and timepoint (for B cells)

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(gage)
library(limma)
library(org.Hs.eg.db)
library(ReactomePA)
library(argparse)

parser <- ArgumentParser(description = "Perform DE tests on B cells")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS for B cells")
parser$add_argument('--min_gene_counts', type='integer', default = 500,
                    help="Minimum number of gene counts to filter at")
parser$add_argument('--patients', type='character', nargs ='+',
                    help="Patients to use", default = NULL)
parser$add_argument('--bcell_labels', type='character', nargs ='+',
                    help="B cell labels", default = NULL)
parser$add_argument('--ngene', type='integer',
                    help="Number of top genes to use", default = 50)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for DE results.")

add_symbol <- function (de_result, filter_ribo = TRUE, filter_mito = TRUE, 
    rowdat) 
{
    de_result <- lapply(de_result, function(x) {
        x$ID <- rownames(x)
        x$Symbol <- df_as_map(rowdat, rownames(x), 
            from = "ID", to = "Symbol")
        if (filter_ribo) {
            x <- subset(x, !str_detect(Symbol, "^RP"))
        }
        if (filter_mito) {
            x <- subset(x, !str_detect(Symbol, "^MT"))
        }
        return(x)
    })
    return(de_result)
}

de_analysis <- function (sce, de_method = "voom", formula, cluster, filter_mito = TRUE, 
    filter_ribo = TRUE, block = NULL, coef = NULL, min_gene_counts = 500, 
    full.stats = FALSE) 
{
    sce <- sce[rowData(sce)$total_counts >= min_gene_counts, 
        ]
    if (de_method == "voom") {
        design <- model.matrix(formula, data = colData(sce))
        dge <- edgeR::DGEList(counts(sce))
        dge <- edgeR::calcNormFactors(dge)
        v <- limma::voom(dge, design, plot = TRUE)
        fit <- lmFit(v, design, robust = TRUE)
        fit <- eBayes(fit)
        diffexpr_tab <- limma::topTable(fit, coef = coef, adjust.method = "BH", 
            number = Inf)
        diffexpr_tab_labeled <- add_symbol(list(diffexpr_tab), 
            filter_mito = filter_mito, filter_ribo = filter_ribo, 
            rowData(sce) %>% as.data.frame)[[1]] %>% dplyr::mutate(FDR = adj.P.Val, 
            is_significant = FDR <= 0.05)
    }
    else if (de_method == "scran") {
        diffexpr_tab <- findMarkers(sce, clusters = cluster, 
            block = block, full.stats = full.stats)
        if (full.stats) {
            diffexpr_tab <- lapply(diffexpr_tab, function(df) {
                idx <- str_detect(colnames(df), "stats\\.")
                stats_cols <- colnames(df)[idx]
                const_cols <- colnames(df)[!idx]
                combined_stats <- do.call("cbind", lapply(stats_cols, 
                  function(x) {
                    clust_name <- str_extract(x, "(?<=\\.).$")
                    stat_df <- data.frame(df[[x]]$logFC, df[[x]]$log.FDR, 
                      df[[x]]$log.p.value)
                    colnames(stat_df) <- paste0(c("logFC", "logFDR", 
                      "logpvalue"), ".", clust_name)
                    return(stat_df)
                  }))
                combined_stats <- data.frame(df[, const_cols], 
                  combined_stats)
                rownames(combined_stats) <- rownames(df)
                return(combined_stats)
            })
        }
        diffexpr_tab_labeled <- add_symbol(diffexpr_tab, filter_mito = filter_mito, 
            filter_ribo = filter_ribo, rowData(sce))
        diffexpr_tab_labeled <- plyr::rbind.fill(lapply(names(diffexpr_tab_labeled), 
            function(x) {
                data.frame(contrast = x, diffexpr_tab_labeled[[x]])
            }))
    }
    else {
        stop("Unrecognized DE method")
    }
    return(list(de_table = diffexpr_tab_labeled, universe_ids = unique(rownames(sce))))
}

df_as_map <- function (df, vec, from, to) 
{
    if (length(from) > 1 || length(to) > 1) {
        stop("Please specify ONE index for from and to, each.")
    }
    df <- df[df[, from] != "", ]
    map <- data.frame(from = df[, from], to = df[, to], stringsAsFactors = FALSE)
    map <- unique(map)
    map <- map[complete.cases(map), ]
    map <- data.table::data.table(map)
    data.table::setkey(map, "from")
    return(map[as.character(vec), ]$to)
}

args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# normal b clusters
normal <- c('1','13','12')

sce <- sce[,sce$celltype == 'malignant_b']
sce <- sce[,sce$timepoint %in% c('FL','DLBCL')]
#sce_filtered = sce[,!(sce$cluster %in% normal)]

# Filter by patient
if (!is.null(args$patients)) {
  object <- object %>%
    scater::filter(patient %in% unlist(args$patients))
}

coefficients <- c("timepointFL")

labels <- c("timepoint")

de_res <- lapply(coefficients, function(coef) {
  voom_res <- de_analysis(object, de_method = "voom", formula = ~ timepoint, 
                          cluster = NULL, 
                          filter_mito = TRUE, filter_ribo = TRUE, block = NULL, coef = coef, 
                          min_gene_counts = args$min_gene_counts)
  
  up_genes <- (voom_res$de_table %>% dplyr::filter(logFC > 0, FDR < 0.05) %>% dplyr::arrange(-logFC))$Symbol
  down_genes <- (voom_res$de_table %>% dplyr::filter(logFC < 0, FDR < 0.05) %>% dplyr::arrange(logFC))$Symbol
  background_genes <- mapIds(org.Hs.eg.db, df_as_map(rowData(object), voom_res$universe_ids, from = "ID", to = "Symbol"), 'ENTREZID', 'SYMBOL')
  
  up_genes <- up_genes %>% head(args$ngene)
  down_genes <- down_genes %>% head(args$ngene)
  
  enriched_paths <- lapply(list(up_genes, down_genes), function(genes) {
    if (length(genes) > 0) {
      gene_entrez <- mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
      paths <- enrichPathway(gene=unname(gene_entrez),pvalueCutoff=0.05, readable=TRUE, 
                             universe = unname(background_genes))
    } else {
      paths <- NULL
    }
    
    return(paths)
  })
  names(enriched_paths) <- c("up", "down")
  
  return(list(gene=voom_res$de_table,
              pathway=enriched_paths))
})

names(de_res) <- labels

# Save DE result (pathway and gene tables)
saveRDS(de_res, args$outfname)

cat("Completed.\n")



