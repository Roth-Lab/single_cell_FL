#' Formats segment file

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)

#library(scrna.utils)
#library(scrna.sceutils)
#library(cellassign.utils)
library(argparse)

library(feather)

parser <- ArgumentParser(description = "Select genes")
parser$add_argument('--gene_cn', metavar='FILE', type='character',
                    help="Gene CN file (feather)")
parser$add_argument('--pct_pure', type = 'double', default = 0.6,
                    help="Threshold for purity of copy number within each clone, time pair")
parser$add_argument('--diff', type = 'double', default = 1,
                    help="Threshold for copy number difference considered too variable between timepoints")

parser$add_argument('--outfig', type = 'character', metavar = 'FILE',
                    help="Output purity figure")

parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for selected gene-CN file.")
args <- parser$parse_args()

gene_cn_summarized <- feather::read_feather(args$gene_cn)

# Remove genes with multiple CN's that differ by more than 1
bad_genes <- gene_cn_summarized %>%
  dplyr::filter((copy_number_max - copy_number_min) > 1)
gene_cn_summarized_filtered <- gene_cn_summarized %>%
  dplyr::filter(!ensembl_gene_id %in% bad_genes$ensembl_gene_id)

calc_mode <- function(x) {
  keys <- unique(x)
  keys[which.max(tabulate(match(x, keys)))]
}

cluster_time_cn_summarized <- gene_cn_summarized_filtered %>%
  dplyr::group_by(clone_id, ensembl_gene_id) %>%
  dplyr::summarise(mean_cnmean=mean(copy_number_mean),
                   mean_cnmode=mean(copy_number_mode),
                   median_cnmode=median(copy_number_mode),
                   mode_cnmode=calc_mode(copy_number_mode),
                   pct_pure=sum(copy_number_mode == mode_cnmode & copy_number_min == copy_number_max)/n(),
                   n_cell_time=n())

timepoint_variability <- cluster_time_cn_summarized %>% 
  dplyr::group_by(clone_id, ensembl_gene_id) %>%
  dplyr::summarise(min_mediancnmode=min(median_cnmode),
                   max_mediancnmode=max(median_cnmode),
                   diff=max_mediancnmode - min_mediancnmode)

cluster_cn_summarized <- gene_cn_summarized_filtered %>%
  dplyr::group_by(clone_id, ensembl_gene_id) %>%
  dplyr::summarise(mean_cnmean=mean(copy_number_mean),
                   mean_cnmode=mean(copy_number_mode),
                   median_cnmode=median(copy_number_mode),
                   mode_cnmode=calc_mode(copy_number_mode),
                   pct_pure=sum(copy_number_mode == mode_cnmode & copy_number_min == copy_number_max)/n(),
                   n_cell=n())

cluster_variability <- cluster_cn_summarized %>%
  dplyr::mutate(l_median_cnmode=log(median_cnmode+1)) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  dplyr::summarise(mad=mad(l_median_cnmode),
                   sd=sd(l_median_cnmode),
                   max=max(median_cnmode),
                   min=min(median_cnmode),
                   range=max-min,
                   range_log=max(l_median_cnmode)-min(l_median_cnmode))

impure_genes <- unique(cluster_time_cn_summarized$ensembl_gene_id[cluster_time_cn_summarized$pct_pure < args$pct_pure])
intertimepoint_variable_genes <- unique(timepoint_variability$ensembl_gene_id[timepoint_variability$diff > 1])

exclude_genes <- union(impure_genes, intertimepoint_variable_genes)

selected_cn <- cluster_variability %>% 
  dplyr::filter(range >= 0,
                !ensembl_gene_id %in% exclude_genes) %>% 
  dplyr::arrange(-mad)

cluster_cn_summarized <- cluster_cn_summarized %>%
  dplyr::mutate(use_gene = ensembl_gene_id %in% selected_cn$ensembl_gene_id)

write.table(cluster_cn_summarized, args$outfname, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
## TODO: another option: look for long continguous blocks of CN

ggplot(cluster_cn_summarized, aes(x = pct_pure)) + 
  geom_histogram() + 
  facet_wrap(~ clone_id) +
  xlab("% pure")

ggsave(args$outfig, width = 7, height = 5.5)

cat("Completed.\n")
