# merge_marker_data.R

library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
#library(VennDiagram)
library(UpSetR)
library(reshape2)
library(tidyverse)
library(grid) 
#futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# IHC panel
panel='IHC_CD70-CD27_panel'
# working.dir should contain folders for each marker
working.dir <- paste0("/path/imaging_data/",panel,"/")

# core info should be a dataframe for the IHC design 
# (with columns: core for core id; study_id for sample id; timepoint for timepoint info)
core_info <- read.table(paste0('/path/imaging_data/core_info_',panel,'.txt'),sep = '\t', header = TRUE)

# sample.data should be a dataframe with one column for the paths for the cell_seg_data for each core
sample.data <- read.delim(paste0(working.dir, "sample_info.txt"), stringsAsFactors = FALSE, header = FALSE)

colnames(sample.data) <- c('path')

# string parsing to get the core id
tp=colsplit(sample.data$path,'/',names = c(1:10))
tp2=colsplit(tp$`10`,'_',names = c(1:7))
#sample.data$core=tp2$`4` # CD10-CD20
sample.data$core=tp2$`3` # CD70-CD27 
sample.data$IHC_ID=gsub('_cell_seg_data.txt','',tp$`10`)

# match by core info
#sample.data <- sample.data %>% left_join(core_info, by = 'core')
#sample.data <- filter(sample.data, core != "", study_id != "")

#markers <- c("CD20", "CD4", "CD8", "FoxP3", "LAG3", "PD1")
#markers <- list.files(paste0(working.dir, "data/Inform_Data"))
markers <- list.files(working.dir)
markers <- markers[!(markers %in% c('merged_data', 'sample_info.txt','plots'))]
#markers <- c("CD4", "CD8", "CD20", "FoxP3", "LAG3", "PD-1")
#markers <- c("CD4", "CD8", "CD20", "FoxP3", "LAG3") # combined_panel no PD1

all.data <- data.frame()
for (i in 1:nrow(sample.data)) {
	print(paste0("Processing sample: ", sample.data$IHC_ID[i]))
	#print(paste0("Processing sample: ", sample.data$core[i]))
	s.data <- data.table()
	m_old <- NULL
	for (m in markers) {
		#m.cellTablePath <- paste0(working.dir, "data/Inform_Data/", m, "/", sample.data$Core_ID[i], "_cell_seg_data.txt")
		m.cellTablePath <- sample.data$path[i]
		
		if (!(is.null(m_old))) {
			m.cellTablePath <- gsub(m_old, m, m.cellTablePath)
			m_old = m
		}

		m.data <- fread(m.cellTablePath)
		m.data <- cbind(Name = sample.data$IHC_ID[i], m.data, Marker_pos = ifelse(m.data$Phenotype == paste0(m, "+"), "POS", "Other"))
		colnames(m.data)[ncol(m.data)] <- m
		m.data <- m.data[,c("Name", "Sample Name", "Tissue Category", "Cell X Position", "Cell Y Position", m), with = FALSE]
		colnames(m.data) <- gsub(" ", "_", colnames(m.data))
		
		# Remove rows with duplicated coordinates, since they can't be distinguished as individual cells
		m.comb <- paste(m.data$Cell_X_Position, m.data$Cell_Y_Position, sep = "_")
		m.dup <- m.comb[which(duplicated(m.comb))]
		if (length(m.dup) > 0) {
			print(paste0("Note: ", length(m.dup), " cells with duplicated coordinates removed"))
			m.data <- m.data[-which(m.comb %in% m.dup),]
		}
		
		if (nrow(s.data) < 1) { 
			s.data <- m.data 
		} else { 
			s.data <- merge(s.data, m.data)
			if (nrow(s.data) != nrow(m.data)) {
				print(paste0("Warning: some non-overlapping x/y coordinates identified, ", nrow(m.data) - nrow(s.data), " cells dropped"))
			}
		}
		
		rm(m.cellTablePath, m.data, m.comb, m.dup)
	}
	s.data <- cbind(s.data, "Cell_ID" = as.character(rownames(s.data)))
	if (nrow(all.data) < 1) { all.data <- s.data } else { all.data <- rbind(all.data, s.data) }
	rm(s.data)
}

all.data <- all.data %>% mutate(cell_label = paste(Name, Cell_ID, sep = "_"))

# Write table for future use
merged.dir <- paste0(working.dir, "merged_data/")
dir.create(merged.dir, recursive = TRUE, showWarnings = FALSE)
write.table(all.data, file = paste0(merged.dir, "all_markers_merged.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

# Plot summary UpSet to show positivity of markers in all cells
all.pos <- lapply(all.data[,markers, with = FALSE], function(x) { return(all.data$cell_label[which(x == "POS")]) })
#all.upset <- upset(fromList(all.pos), nsets = length(markers), nintersects = NA, order.by = "freq", 
#									 mainbar.y.label = "Number of cells", number.angles = 45, number.hjust = 0, number.vjust = -0.5)

all.upset <- upset(fromList(all.pos), nsets = length(markers), nintersects = NA, order.by = "freq", 
									 mainbar.y.label = "Number of cells", number.angles = 45)


plot.dir <- paste0(working.dir, "plots/")
dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)
pdf(paste0(plot.dir, "all_cells_marker_coexpression.pdf"), width = 40, height = 5)
print(all.upset)
grid.text("All cells", x = 0.65, y = 0.95, gp = gpar(fontsize = 8))
dev.off()
