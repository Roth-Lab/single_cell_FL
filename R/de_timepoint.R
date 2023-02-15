# Differential expression (grouped by timepoint)

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
library(fgsea)
library(argparse)

parser <- ArgumentParser(description = "Perform DE tests on an SCE for between FL and DLBCL timepoint")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--method_gene', type='character',
                    help="DE method to use", default = "voom")
parser$add_argument('--method_pathway', type='character',
                    help="DE method to use", default = "ReactomePA")
parser$add_argument('--gene_set_file', type='character', metavar = "FILE",
                    help="Gene set file path, if using fgsea", default = NULL)
parser$add_argument('--ngene', type='integer',
                    help="Number of top genes to use", default = 50)
parser$add_argument('--min_gene_counts', type='integer', default = 500,
                    help="Minimum number of gene counts to filter at")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for DE tables.")

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
        #print(fit$coefficients)
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

args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)
method_gene <- args$method_gene

# normal b clusters
normal <- c('1','13','12')

sce <- sce[,sce$celltype == 'malignant_b']
sce <- sce[,sce$timepoint %in% c('FL','DLBCL')]
#sce_filtered = sce[,!(sce$cluster %in% normal)]


## select highly variable genes
#fit <- trendVar(sce, parametric=TRUE, use.spikes = FALSE)
#decomp <- decomposeVar(sce, fit)
decomp <- modelGeneVar(sce)
decomp$Symbol <- rowData(sce)$Symbol
decomp_chosen <- decomp %>% subset(bio > 0)
chosen <- rownames(decomp_chosen)

sce_filtered=sce[chosen,]

# model formula
cont <- "FL"
coef <- "time_pointDLBCL"
#coef <- "time_pointFL_cont"
time_point=relevel(factor(sce_filtered$timepoint),'FL')
formula <- ~batch + patient + time_point
cluster_col <- "timepoint"


if (method_gene == "voom") {
  voom_res <- de_analysis(sce_filtered, de_method = "voom", formula, cluster = NULL, 
                          filter_mito = TRUE, filter_ribo = TRUE, block = NULL, coef = coef, 
                          min_gene_counts = args$min_gene_counts)
  de_table <- voom_res$de_table
  
  up_genes <- (voom_res$de_table %>% dplyr::filter(logFC > 0, FDR < 0.05) %>% dplyr::arrange(-logFC))$Symbol
  down_genes <- (voom_res$de_table %>% dplyr::filter(logFC < 0, FDR < 0.05) %>% dplyr::arrange(logFC))$Symbol
  background_genes <- mapIds(org.Hs.eg.db, df_as_map(rowData(sce), voom_res$universe_ids, from = "ID", to = "Symbol"), 'ENTREZID', 'SYMBOL')
} else if (method_gene == "scran") {
  scran_res <- de_analysis(sce_filtered, de_method = "scran", formula = NULL, cluster = colData(sce_filtered)[,cluster_col], 
                           filter_mito = TRUE, filter_ribo = TRUE, block = NULL, coef = NULL, 
                           min_gene_counts = args$min_gene_counts)
  de_table <- scran_res$de_table
  
  up_genes <- (scran_res$de_table %>% dplyr::filter(contrast == cont, FDR < 0.05, logFC.DLBCL > 0) %>% dplyr::arrange(-logFC.DLBCL))$Symbol
  down_genes <- (scran_res$de_table %>% dplyr::filter(contrast == cont, FDR < 0.05, logFC.DLBCL < 0) %>% dplyr::arrange(logFC.DLBCL))$Symbol
  background_genes <- mapIds(org.Hs.eg.db, df_as_map(rowData(sce), scran_res$universe_ids, from = "ID", to = "Symbol"), 'ENTREZID', 'SYMBOL')
} else {
  stop("Unrecognized DE method.")
}


### Enrichment ###
if (args$method_pathway == "ReactomePA") {
  up_genes <- up_genes %>% head(args$ngene)
  down_genes <- down_genes %>% head(args$ngene)
  
  enriched_paths <- lapply(list(up_genes, down_genes), function(genes) {
    gene_entrez <- mapIds(org.Hs.eg.db, genes, 'ENTREZID', 'SYMBOL')
    paths <- enrichPathway(gene=unname(gene_entrez),pvalueCutoff=0.05, readable=TRUE, 
                           universe = unname(background_genes))
    
    return(paths)
  })
  names(enriched_paths) <- c("up", "down")
} else if (args$method_pathway == "fgsea") {
  # Only allow use of scran with this right now -- can use scran but need to rename columns accordingly
  #stopifnot(method_gene == "voom")
  stopifnot(method_gene == "scran")
  
  if (!is.null(args$gene_set_file)) {
    pathway_list <- gmtPathways(args$gene_set_file)
  } else {
    stop("Must specify a gene set file.")
  }
  
  de_table_summarized <- de_table %>%
    dplyr::filter(contrast == "FL") %>%
    dplyr::select(Symbol, logFC.DLBCL) %>%
    distinct() %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(logFC = mean(logFC.DLBCL))
  
  # de_table_summarized <- de_table %>% 
  #   dplyr::select(Symbol, t) %>% 
  #   na.omit() %>% 
  #   distinct() %>% 
  #   dplyr::group_by(Symbol) %>% 
  #   dplyr::summarise(t = mean(t))
  # 
  fgsea_res <- fgsea(pathways=pathway_list, 
                     stats = deframe(de_table_summarized), 
                     nperm = 10000)
  
  enriched_paths <- fgsea_res %>%
    as_tibble() %>%
    dplyr::arrange(desc(NES))
} else {
  stop("Unrecognized pathway enrichment method.")
}


# Save DE result (pathway and gene tables)
saveRDS(list(gene=de_table, pathway=enriched_paths), args$outfname)
#saveRDS(voom_res, args$outfname)

cat("Completed.\n")



