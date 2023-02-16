#' Run CloneAlign
library(argparse)
library(annotables)

parser <- ArgumentParser(description = "Run CloneAlign")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="SingleCellExperiment")
parser$add_argument('--cnv', metavar='FILE', type='character',
                   help="CNV gene file")
# parser$add_argument('--clone_prevs', metavar='FILE', type='character',
#                     help="Clonal prevalences file")
parser$add_argument('--gex_var_quantile', type='double',
                    help="Gene expression quantile", default = 0.5)
parser$add_argument('--n_extra_genes', type='integer',
                    help="Number of extra copy-same genes to add", default = 50)
parser$add_argument('--samples', type = 'character', nargs = '+',
                    help="Sample IDs to subset")
parser$add_argument('--initial_shrink', type = 'double', default = 0,
                    help="Initial shrink")
parser$add_argument('--conda_env', type = 'character',
                    help="Name of conda environment with tensorflow", default = "r-tensorflow")
parser$add_argument('--conda_path', type = 'character',
                    help="Conda path", default = "auto")
parser$add_argument('--clones_to_ignore', type = 'character',
                    help="Clones to ignore (comma separated)", default = "auto")
parser$add_argument('--max_copy_number', type='integer',
                    help="Maximum copy number for input to clonealign", default = 11)
parser$add_argument('--mean_counts', type='double',
                    help="Min mean counts per gene to be included after filtering", default = 100)
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for fit.")
args <- parser$parse_args()

# Reset snakemake default
# Sys.setenv(PYTHONPATH='')

#reticulate::use_condaenv("r-reticulate")
reticulate::use_condaenv(args$conda_env, 
                         conda = args$conda_path,
                         required = TRUE)

print(reticulate::py_config())


suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(data.table)
  library(methods)
  library(scran)
  library(yaml)  
  library(clonealign)
})

set.seed(12345)

# devtools::load_all("/cellassign/clonealign/")


sce <- readRDS(args$sce)
rownames(sce) <- rowData(sce)$ID

raw_cnvs <- fread(args$cnv)
#clone_prevalences <- fread(args$clone_prevs)
#gex_var_quantile <- args$gex_var_quantile



## Remove clones to ignore
if(args$clones_to_ignore != "None") {
  clones_to_ignore <- strsplit(args$clones_to_ignore, ",")[[1]]
  raw_cnvs <- raw_cnvs[!(raw_cnvs$clone_id %in% clones_to_ignore),]
}

## Get list of timepoints to evaluate
sce <- sce[, sce$timepoint != 'RLN']
sce <- sce[, sce$celltype == 'malignant_b']

#sce <- sce %>%
#  scater::filter(patient %in% unlist(args$samples))
sce=sce[,sce$patient %in% unlist(args$samples)]
remove(cds)

timepoints <- unique(sce$timepoint)

### get only fl and only malignant B cells ######
#sce <- sce %>%
#  scater::filter(timepoint=='FL')
#sce <- sce %>%
#  scater::filter(celltype_full=='B cells (malignant)')

#colnames(sce) <- paste0(sce$id, "_", sce$Barcode)
colnames(sce) <- paste0(sce$patient, "_", sce$Barcode)


## Important -- select autosomes only

autosomal_genes <- annotables::grch37 %>% 
  dplyr::select(ensembl_gene_id = ensgene, chr) %>% 
  dplyr::filter(chr %in% as.character(1:22)) %>% 
  .$ensembl_gene_id


#save.image(paste0("/cellassign/fitness-scrna/results/dumps/_clonealign-dump-", sce$id[1], ".Rdata"))
save.image(paste0("_clonealign-dump-", sce$patient[1], ".Rdata"))

# clone_prevalences <- clone_prevalences %>%
#   dplyr::filter(time %in% timepoints)

# if(nrow(clone_prevalences) == 0) {
#   stop("No matching clone prevalences!")
# }

# present_clones <- unique(clone_prevalences$cluster)

## Filter CNV data
cnv <- dplyr::filter(raw_cnvs, use_gene) %>%
  dplyr::rename(clone = clone_id,
                copy_number=median_cnmode) %>% 
  dplyr::select(ensembl_gene_id, clone, copy_number) %>% 
  # dplyr::filter(clone %in% present_clones) %>%
  spread(clone, copy_number)

cnv_mat <- cnv %>%
  as.data.frame %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix

common_genes <- intersect(intersect(rownames(cnv_mat), rownames(sce)), autosomal_genes)
sce <- sce[common_genes,]
cnv_mat <- cnv_mat[common_genes,]

processed_data <- preprocess_for_clonealign(sce, 
                                            cnv_mat,
                                            min_counts_per_gene = args$mean_counts)
                                            #mean_counts_per_gene = args$mean_counts)

#saveRDS(processed_data,'~/Desktop/process_data.rds')
## remove sce for memory
#remove(sce)

cnv_mat_filtered <- processed_data$copy_number_data

dist_matrix <- as.matrix(dist(t(cnv_mat_filtered), method = "manhattan"))/nrow(cnv_mat_filtered)
  
clone_remaps <- reshape2::melt(dist_matrix) %>% 
  dplyr::filter(value < 0.005) %>% 
  dplyr::group_by(Var1) %>% 
  dplyr::summarise(merged_id = paste(sort(unique(as.character(Var2))), collapse = "_"))

cnv_mat_remapped <- cnv_mat_filtered[,clone_remaps$Var1[!duplicated(clone_remaps$merged_id)],drop=FALSE]

colnames(cnv_mat_remapped) <- colnames(cnv_mat_remapped) %>%
  plyr::mapvalues(from = clone_remaps$Var1, clone_remaps$merged_id)


processed_data <- preprocess_for_clonealign((processed_data$gene_expression_data), 
                                            cnv_mat_remapped,
                                            min_counts_per_gene = args$mean_counts)
                                          #min_counts_per_cell = 0.05)
                                          #mean_counts_per_gene = args$mean_counts)


print(paste("Final dimensions expression data", dim(processed_data$gene_expression_data)))

# mode(processed_data$gene_expression_data)='integer'
# mode(processed_data$copy_number_data)='integer'

ca <- run_clonealign(processed_data$gene_expression_data, 
                     processed_data$copy_number_data,
                     initial_shrinks = c(2,5,10), #args$initial_shrink,
                     n_repeats = 3,
                     mc_samples = 1,
                     learning_rate = 0.07,
                     max_iter = 5e2,
                     data_init_mu = TRUE,
                     saturation_threshold = args$max_copy_number,
                     clone_call_probability = 0.9
                     )

#sce_tmp=sce[processed_data$retained_genes,]
#sce_tmp$keep=colSums(as.matrix(counts(sce_tmp)))>0
#sce_tmp=sce_tmp %>% scater::filter(keep==TRUE)
#
#cnv_tmp=cnv_mat[processed_data$retained_genes,]
#
#ca <- run_clonealign(sce, 
#                     cnv_mat,
#                     initial_shrinks = c(2,5,10), #args$initial_shrink,
#                     n_repeats = 3,
#                     mc_samples = 1,
#                     learning_rate = 0.07,
#                     max_iter = 5e2,
#                     data_init_mu = TRUE,
#                     saturation_threshold = args$max_copy_number,
#                     clone_call_probability = 0.9,
#                     seed = 12345)


#sce_filtered <- sce[, rownames(processed_data$gene_expression_data)]
#sce_filtered <- sce

#saveRDS(ca,'~/Desktop/ca_tmp.rds')

sce_filtered <- sce[, rownames(processed_data$gene_expression_data)]


clone_fit <- data.frame(
  Sample=sce_filtered$Sample,
  Barcode=sce_filtered$Barcode,
  #cell_id=rownames(processed_data$gene_expression_data),
  cell_id=colnames(sce_filtered),
  id=sce_filtered$dataset,
  clone=ca$clone,
  stringsAsFactors = FALSE
)

ca$clone_fit <- clone_fit
#ca$cnv_mat <- cnv
ca$cnv_mat <- cnv_mat_remapped

# Save results
saveRDS(ca, args$outfname)

cat("Completed.\n")




