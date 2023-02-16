#' Formats segment file

library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(data.table)
library(methods)
library(scran)
library(yaml)
library(reshape2)
#library(scrna.utils)
#library(scrna.sceutils)
#library(cellassign.utils)
library(argparse)

parser <- ArgumentParser(description = "Select genes")
parser$add_argument('--segment_files', metavar='FILE', type='character', nargs ='+',
                    help="Segment files")
parser$add_argument('--outfname', type = 'character', metavar = 'FILE',
                    help="Clone prevalences")
args <- parser$parse_args()

raw_segments <- dplyr::bind_rows(lapply(args$segment_file, function(f) {
  fread(f)
}))

barcode=colsplit(string=raw_segments$cell_id, pattern="-", names=c('time','2','3','4'))
raw_segments$time=barcode$time

segments <- raw_segments %>%
  dplyr::mutate(chr=paste0("chr", chr))

prevalences <- segments %>%
  dplyr::select(cell_id, time, clone_id) %>%
  unique %>%
  dplyr::group_by(clone_id, time) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq=n/sum(n))

write.table(prevalences, file = args$outfname, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

cat("Completed.\n")
