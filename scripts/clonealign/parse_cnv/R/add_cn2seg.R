#! /usr/bin/env Rscript
# add clone information to segment files

library(tidyverse)
# library(tensorflow)
# library(SingleCellExperiment)
# library(scater)
# library(data.table)
# library(methods)
# library(scran)
# library(yaml)

#library(scrna.utils)
#library(scrna.sceutils)
#library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Add clone information to segment files")
parser$add_argument('--segment_file', metavar='FILE', type='character',
                    help="Segment file")
parser$add_argument('--clone_file', metavar='FILE', type='character',
                    help="Clone files")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Output path for the new segment file.")
args <- parser$parse_args()


segment <- read.csv(args$segment_file, header=T)
clone <- read.csv(args$clone_file, header=T)

new_segment <- segment %>% 
  dplyr::inner_join(clone, by='cell_id')

# Save the new file
write.csv(new_segment, args$outfname)

cat('Completed.\n')
