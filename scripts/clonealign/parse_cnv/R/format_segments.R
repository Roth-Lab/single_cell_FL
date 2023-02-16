#' Formats segment file

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(parallel)

#library(scrna.utils)
#library(scrna.sceutils)
#library(cellassign.utils)
library(argparse)

library(feather)

parser <- ArgumentParser(description = "Format segments")
parser$add_argument('--segment_file', metavar='FILE', type='character', nargs = '+',
                    help="Segment file")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for gene CN table (feather file).")
args <- parser$parse_args()

raw_segments <- dplyr::bind_rows(lapply(args$segment_file, function(f) {
  fread(f)
}))

# add clone information to segments


# Sample max 1000 cells

all_cells <- unique(raw_segments$cell_id)

set.seed(543)
cells_to_sample <- sample(all_cells, min(1000, length(all_cells)), replace = FALSE)

raw_segments <- raw_segments[raw_segments$cell_id %in% cells_to_sample,]

segments <- raw_segments %>%
  dplyr::mutate(chr=paste0("chr", chr))

# Get gene coordinates
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
g <- genes(txdb, single.strand.genes.only=FALSE)

segments_gr <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
overlaps <- findOverlaps(g, segments_gr, ignore.strand = TRUE)
overlapping_ranges <- pintersect(g[queryHits(overlaps)], segments_gr[subjectHits(overlaps)], ignore.strand = TRUE, drop.nohit.ranges = TRUE)

percentage_overlap <- width(overlapping_ranges) / width(segments_gr[subjectHits(overlaps)])
percentage_overlap_sum <- sapply(percentage_overlap, sum)

gene_cn <- tibble(entrezgene = names(g)[queryHits(overlaps)], percent_overlap=percentage_overlap_sum)
gene_cn$ensembl_gene_id <- mapIds(org.Hs.eg.db, 
                                  keys = gene_cn$entrezgene, 
                                  column="ENSEMBL", 
                                  keytype="ENTREZID", 
                                  multiVals="first")

gene_cn <- gene_cn %>%
  cbind(mcols(segments_gr[subjectHits(overlaps)]))

# Filter for complete entries
gene_cn <- gene_cn %>%
  as.data.frame %>%
  drop_na()

# Make mode computation weighted by percent_overlap easier by doing this in two steps

## Parallel this

unique_genes <- unique(gene_cn$ensembl_gene_id)

gc()

cat("Summarizing to percent overlap...\n")

gene_cn_summarized <- mclapply(unique_genes, function(g) {
  dplyr::filter(gene_cn, ensembl_gene_id == g)  %>%
  dplyr::group_by(ensembl_gene_id, cell_id, clone_id, state) %>%
  #dplyr::group_by(ensembl_gene_id, cell_id, time, cluster, copy_number) %>%
  dplyr::summarise(percent_overlap=sum(percent_overlap),
                   n_parts=n()) %>% 
    ungroup()
}, mc.cores = 8) %>% bind_rows()


# gene_cn_summarized <- gene_cn_summarized %>% 
#   dplyr::group_by(ensembl_gene_id, cell_names, time, cluster) %>%
#   dplyr::summarise(copy_number_mean=sum(copy_number*percent_overlap)/sum(percent_overlap),
#                    copy_number_min=min(copy_number),
#                    copy_number_max=max(copy_number),
#                    copy_number_mode=copy_number[which.max(percent_overlap)][1],
#                    n_parts=sum(n_parts))

cat("Garbage collection...\n")

rm(
  list = setdiff(ls(), c("gene_cn_summarized", "unique_genes", "args"))
)

gc()

cat("Summarizing to gene specific CN...\n")

gene_cn_summarized2 <- mclapply(unique_genes, function(g) {
  dplyr::filter(gene_cn_summarized, ensembl_gene_id == g)  %>%
  dplyr::group_by(ensembl_gene_id, cell_id, clone_id) %>%
  #dplyr::group_by(ensembl_gene_id, cell_id, time, cluster) %>%
  dplyr::summarise(copy_number_mean=sum(state*percent_overlap)/sum(percent_overlap),
                   copy_number_min=min(state),
                   copy_number_max=max(state),
                   copy_number_mode=state[which.max(percent_overlap)][1],
                   n_parts=sum(n_parts)) %>% 
    ungroup()
}, mc.cores = 8) %>% bind_rows()


feather::write_feather(gene_cn_summarized2, args$outfname)

cat("Completed.\n")

