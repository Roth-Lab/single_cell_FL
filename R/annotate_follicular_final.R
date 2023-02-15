# Annotate FL SCE with unsupervised cluster assignments and cyclone results

library(tidyverse)
library(tensorflow)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)
library(monocle3)

#library(scrna.utils)
#library(scrna.sceutils)
#library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Annotate SCE with unsupervised clustering and cyclone results")
parser$add_argument('--sce', metavar='FILE', type='character',
                    help="Path to SingleCellExperiment RDS")
parser$add_argument('--cyclone', metavar='FILE', type='character',
                    help="Path to cyclone results")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for annotated CDS")
args <- parser$parse_args()

sce_path <- args$sce
sce <- readRDS(sce_path)

# Read in cell cycle assignments
cyclone_results <- readRDS(args$cyclone)

# Add cell cycle information
sce@colData <- bind_cols(sce@colData %>% as.data.frame, 
                         cyclone_results$normalized.scores) %>% 
  DataFrame(check.names = FALSE)

sce <- sce %>%
  scater::mutate(Cell_Cycle = cyclone_results$phases)

# convert sce object to monocle3 cds object to downstream analysis
  #### construct cds object #####
cell_dat <- data.frame(colData(sce))
rownames(cell_dat) <- cell_dat$sample_barcode

gene_dat <- data.frame(rowData(sce))
rownames(gene_dat) <- gene_dat$ID
gene_dat$gene_short_name <- gene_dat$Symbol

cds <- new_cell_data_set(counts(sce),cell_metadata = cell_dat,gene_metadata = gene_dat)

# store no-batch-correction result into LSI
cds <- preprocess_cds(cds, method='PCA')
cds <- reduce_dimension(cds)
cds <- monocle3::cluster_cells(cds, reduction_method = 'UMAP' ,random_seed = 2)
reducedDims(cds)$LSI <- reducedDims(cds)$UMAP

# add scanorama PCA result
reducedDims(cds)$PCA <- reducedDims(sce)$scanorama_PCA
cds <- reduce_dimension(cds)
cds <- monocle3::cluster_cells(cds, reduction_method = 'UMAP' ,random_seed = 2)
colData(cds)$cluster <- factor(cds@clusters[[1]]$clusters)

# manual annotation result
non_b_assignments <- c(2,15,26,35,5,19,9,10,7,34,4,36,38)
normal_b_assignments <- c(1,12,13)

# add celltype information
colData(cds)$celltype <- ifelse(colData(cds)$cluster %in% non_b_assignments, 'non_b', colData(cds)$cluster)
colData(cds)$celltype <- ifelse(colData(cds)$cluster %in% normal_b_assignments, 'normal_b', colData(cds)$celltype)
colData(cds)$celltype <- ifelse(!(colData(cds)$cluster %in% c(non_b_assignments,normal_b_assignments)), 'malignant_b', colData(cds)$celltype)

# Save annotated SCE
sce$cluster=factor(cds@clusters[[1]]$clusters)
reducedDims(sce)$umap=reducedDims(cds)$UMAP
saveRDS(sce, args$outfname)

# Save annotated CDS
saveRDS(cds, gsub('sce','cds',args$outfname))


cat("Completed.\n")



