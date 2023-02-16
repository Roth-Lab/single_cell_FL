# Install and load required packages
required = c('ggplot2', 'spatstat', 'tiff', 'raster', 'data.table')
for (lib in required)
{
  if (!require(lib, character.only=TRUE))
  {
    install.packages(lib, repos="http://cran.rstudio.com")
    suppressMessages(library(lib, character.only=TRUE, quietly=TRUE))
  }
}

# Class to represent a phenotype and expression level
PhenExp <- setRefClass("PhenExp",
  fields = list(
    name = "character",
    posMarkers = c("character"),
    negMarkers = c("character"),
    category = "character",
    expression = 'character',
    threshold = 'numeric',
    nickname = 'character',
    color = 'character'
  ),
  methods = list(
    show = function() {
      cat('Phenotype:', name, nickname, '\n')
      cat(expression, '>', threshold, '\n')
    },
    fullname = function() 
      ifelse(name=='All', category, paste(category, name))
  )
)

phenotypeName = 'Phenotype' # Name of the phenotype column in the cell seg data
categoryName = 'Tissue_Category' # Name of the category column
sampleName = 'Sample_Name'
pixelsPerMicron = 2

# Columns that must be in every data set
fixedColumns = c('Cell_ID', 'Cell_X_Position', 'Cell_Y_Position', categoryName, sampleName)
nfixed = length(fixedColumns)

#' Read a cell seg table and compute distance between positively expressing phenotypes.
#' @param cellTablePath Path to the source file
#' @param outPath Path to the result
#' @param pheno1 PhenExp of the 'from' phenotype
#' @param pheno2 PhenExp of the 'to' phenotype
#' @param plot If TRUE, output a plot showing the from and to data with lines to the nearest neigbors
computeExpressionDistance = function(
  cellTablePath, outPath,
  pheno1, pheno2, plot=TRUE
)
{
  # Get the path to the cell seg table and check it
  if (is.null(cellTablePath))
    cellTablePath = file.choose()

  cat('Reading', cellTablePath, '\n')
  data = readCellSegData(cellTablePath)

  phenoData = getPhenotypes(data, pheno1, pheno2)
  
  cat('Computing distances\n')
  nnDist = findExpressionDistance(phenoData)
  
  if (is.null(outPath))
    outPath = sub('.txt', paste0('_dist_', pheno1$name, '_to_', pheno2$name, '.txt'), cellTablePath)
  cat('Writing', outPath, '\n')
  writeCellSegData(nnDist, outPath)
  
  if (plot)
  {
    # Make a plot showing all from and to points with lines connecting nearest neighbors
    plotPath = sub('.txt', '.png', outPath)
    cat('Writing plot to', plotPath, '\n')
    theme_reports()
    p = nnPlot(phenoData, nnDist)
    png(plotPath, width=800, height=600)
    print(p)
    dev.off()
  }
}

#' Split data into expressed phenotypes according to pheno1 and pheno2
getPhenotypes <- function (data, pheno1, pheno2)
{
  return (c(getFromPhenotype(data, pheno1), getToPhenotype(data, pheno2)))
}

getFromPhenotype = function(data, pheno1)
{
  # Check for required data
  stopifnot(all(fixedColumns %in% names(data)))
	markerCols <- c()
	if (all(!is.na(pheno1$posMarkers))) { markerCols <- c(markerCols, pheno1$posMarkers) }
	if (all(!is.na(pheno1$negMarkers))) { markerCols <- c(markerCols, pheno1$negMarkers) }
	stopifnot(all(markerCols %in% colnames(data)))
  if (!is.na(pheno1$expression))
    stopifnot(pheno1$expression %in% names(data))
  
  # Split into one data set for each phenotype and (optional) expression level
  # Keep just the rows corresponding to the desired expression and phenotypes for our data
  # and columns containing fields of interest
  keepRow <- rep(TRUE, nrow(data))
  for (pm in pheno1$posMarkers) {
  	if (is.na(pm)) { next }
  	pm.idx <- data[[pm]] == "POS"
  	keepRow <- keepRow & pm.idx
  }
  for (nm in pheno1$negMarkers) {
  	if (is.na(nm)) { next }
  	nm.idx <- data[[nm]] != "POS"
  	keepRow <- keepRow & nm.idx
  }
  
  if (!is.na(pheno1$category)) {
  	keepRow <- keepRow & data[[categoryName]] == pheno1$category
  }
  
  dataFrom = droplevels(data[keepRow,])
  dataFrom = as.data.frame(dataFrom) # modified

  if (is.na(pheno1$expression)) {
    dataFrom = dataFrom[,c(fixedColumns, markerCols)]
  } else {
    dataFrom = dataFrom[dataFrom[[pheno1$expression]]>=pheno1$threshold, c(fixedColumns, markerCols, pheno1$expression)]
  }

  #names(dataFrom)[1:nfixed] = paste('From', fixedColumns)
  names(dataFrom) = paste('From', names(dataFrom), sep = "_")
  dataFrom = na.omit(dataFrom)
  
  # Convert to point pattern objects. Guess the dimensions...
  ppFrom = with(dataFrom, ppp(From_Cell_X_Position, From_Cell_Y_Position, window=guessWindow(data), 
                              marks=factor(rep(pheno1$fullname(), nrow(dataFrom)))))
  list(dataFrom=dataFrom, ppFrom=ppFrom, pheno1=pheno1)
}

getToPhenotype = function(data, pheno2)
{
  # Check for required data
  stopifnot(all(fixedColumns %in% names(data)))
	markerCols <- c()
	if (all(!is.na(pheno2$posMarkers))) { markerCols <- c(markerCols, pheno2$posMarkers) }
	if (all(!is.na(pheno2$negMarkers))) { markerCols <- c(markerCols, pheno2$negMarkers) }
	stopifnot(all(markerCols %in% colnames(data)))
  if (!is.na(pheno2$expression))
    stopifnot(pheno2$expression %in% names(data))
  
 	keepRow <- rep(TRUE, nrow(data))
	for (pm in pheno2$posMarkers) {
		if (is.na(pm)) { next }
		pm.idx <- data[[pm]] == "POS"
		keepRow <- keepRow & pm.idx
	}
	for (nm in pheno2$negMarkers) {
		if (is.na(nm)) { next }
		nm.idx <- data[[nm]] != "POS"
		keepRow <- keepRow & nm.idx
	}
 	
 	if (!is.na(pheno2$category)) {
 		keepRow <- keepRow & data[[categoryName]] == pheno2$category
 	}
 	
 	dataTo = droplevels(data[keepRow,])
  dataTo = as.data.frame(dataTo) # modified

  if (is.na(pheno2$expression)) {
    dataTo = dataTo[,c(fixedColumns, markerCols)]
  } else {
    dataTo = dataTo[dataTo[[pheno2$expression]]>=pheno2$threshold, c(fixedColumns, markerCols, pheno2$expression)]
  }
  
  #names(dataTo)[1:nfixed] = paste('To', fixedColumns)
  names(dataTo) = paste('To', names(dataTo), sep = "_")
  dataTo = na.omit(dataTo)

  ppTo = with(dataTo, ppp(To_Cell_X_Position, To_Cell_Y_Position, window=guessWindow(data), 
                          marks=factor(rep(pheno2$fullname(), nrow(dataTo)))))
  list(dataTo=dataTo, ppTo=ppTo, pheno2=pheno2)
}

#' Compute nearest neighbor distance between positively expressing phenotypes.
#' @param data
#' @param pheno1 PhenExp of the 'from' phenotype
#' @param pheno2 PhenExp of the 'from' phenotype
#' @return a  data.frame showing ID and location of from and to points with the distance between them
findExpressionDistance = function(phenoData)
{
  d = nncross(phenoData$ppFrom, phenoData$ppTo)
  names(d)[1] = 'Distance'
  result = cbind(phenoData$dataFrom, d, phenoData$dataTo[d$which,])
  result$which = NULL
  droplevels(result)
}

findProportion = function(phenoData, s.data, value)
{

  # remove double positive
  s.data=s.data[!(s.data$CD4=='POS' & s.data$CD8=='POS'),]

  aaa=s.data[!(duplicated(s.data$cell_label)),]
  aaa=as.data.frame(aaa)
  rownames(aaa)=aaa$Cell_ID

  rm(s.data)

  # get CD8+ Cell_ID !!!!!!MODIFIED!!!!!!
  aaa_cd8=aaa[aaa$CD8 == 'POS' & aaa$CD20 != 'POS',]$Cell_ID

  aaa=aaa[,c('Cell_X_Position','Cell_Y_Position')]
  DM= as.matrix(dist(aaa))
  df=melt(DM, varnames = c("source", "target"))
  rm(DM)
  df=df[df$source %in% phenoData$dataFrom$From_Cell_ID,]
  df=df[df$value<=value,]
  df$is_to=ifelse(df$target %in% phenoData$dataTo$To_Cell_ID,1,0)


  # out of CD8+ cells 
  df_tmp=df[df$target %in% aaa_cd8,]
  rm(aaa_cd8)
  df_total=as.data.frame(table(df_tmp$source)) # out of all cells: change df_tmp to df 
  colnames(df_total)=c('source','total')
  
  tmp=aggregate(df_tmp$is_to,by=list(source=df_tmp$source),FUN=sum) # out of all cells: change df_tmp to df 
  colnames(tmp)=c('source','positive')
  df_total$source=as.character(df_total$source)
  tmp$source=as.character(tmp$source)
  df_total=df_total %>% left_join(tmp,by='source')
  df_total$Proportion=df_total$positive / df_total$total

  return(df_total)

}

#' Make a plot of all points in phenoData with lines to nearest neighbors as defined by nnDist
#' @param phenoData result of a call to getPhenotypes
#' @param nnDist result of a call to findExpressionDistance
#' @param lineColor color for the connecting lines
#' @return a ggplot object
nnPlot <- function (phenoData, nnDist, lineColor='gray40') {
  title = bquote('Nearest neighbor from' 
                 ~ italic(.(phenoData$pheno1$name)) 
                 ~ 'to' 
                 ~ italic(.(phenoData$pheno2$name)))
  nnPlotImpl(nnDist, phenoData, title, lineColor)
}

# nnPlot for lists of data
nnPlotMulti <- function (phenoData, nnDist, lineColor='gray40') {
  title = bquote('Nearest neighbor from' 
                 ~ italic(.(phenoData[[1]]$pheno1$fullname())) 
                 ~ 'to' 
                 ~ italic(.(phenoData[[1]]$pheno2$fullname())))
  p = nnPlotBase(title)
  for (i in 1:length(phenoData))
    p = addDistData(p, nnDist[[i]], phenoData[[i]], lineColor)
  p
}

#' Make a plot of mutual nearest neighbors in phenoData - points for which each is the NN of the other
nnPlot2 = function(phenoData, nnDist12, nnDist21, lineColor='darkgray')
{
  # Find mutual pairs by merging nnDist1 and nnDist2
  # First get just the nnDist21 cell ids and rename them to match nnDist12
  nn = nnDist21[, c('From_Cell_ID', 'To_Cell_ID')]
  names(nn) = c('To_Cell_ID', 'From_Cell_ID')
  nn = merge(nn, nnDist12)

  title = bquote('Mutual nearest neighbors for' 
                 ~ italic(.(phenoData$pheno1$fullname())) 
                 ~ 'and\n' 
                 ~ italic(.(phenoData$pheno2$fullname())))
  nnPlotImpl(nn, phenoData, title, lineColor)
}

# Shared implementation for nearest-neighbor plots
nnPlotImpl <- function (nnDist, phenoData, title, lineColor) {
  p = nnPlotBase(title)
  p = addDistData(p, nnDist, phenoData, lineColor)
  p
}

nnPlotBase = function(title)
{
  p = ggplot(data=data.frame(x=0, y=0), aes(x=x, y=y)) # Fake d.f needed to get background to draw...
  p = p + labs(x='Cell X Position', y='Cell Y Position', title=title)
  #addScalesAndBackground(p)
}

addDistData = function(p, nnDist, phenoData, lineColor)
{
  if (nrow(nnDist) > 0)
  {
    p = p + geom_segment(data = nnDist,
                         aes(x=From_Cell_X_Position, y=From_Cell_Y_Position,
                             xend=To_Cell_X_Position, yend=To_Cell_Y_Position), color=lineColor)
    
    p = p + geom_point(data=phenoData$dataFrom, 
                       aes(x=From_Cell_X_Position, y=From_Cell_Y_Position),
                       color=phenoData$pheno1$color,
    									 size=1)
    p = p + geom_point(data=phenoData$dataTo, 
                       aes(To_Cell_X_Position, To_Cell_Y_Position), 
                       color=phenoData$pheno2$color,
    									 size=1)
  }
  p
}

# Add scales, scale line and background image to a ggplot object
# background, xlim, ylim are taken from the environment <cough> <cough>
addScalesAndBackground = function(p)
{
  # Add scales at the image limits. Reverse the y scale to match the image
  p = p + scale_x_continuous(limits=c(0, xlim)) + scale_y_reverse(limits=c(ylim, 0))
  
  # Force square aspect ratio
  p = p + coord_fixed()
  
  # Add background image if we have one
  if (exists('background'))
  {
    p = p + annotation_raster(background, xmin=0, xmax=xlim, ymin=0, ymax=-ylim)
  }
  
  # Add a 200-micron line segment for scale reference
  p = p + geom_segment(aes(x=xlim-50-200, xend=xlim-50, y=ylim-100, yend=ylim-100), color='black', size=1)
  p = p + geom_text(aes(x=xlim-50-200/2, y=ylim-90, label=paste(200, '~mu*m')), 
                    size=3, hjust=0.5, vjust=1, color='black', parse=TRUE)
  p
}

# Guess the actual window size in microns by guessing it is a Vectra / Nuance image
guessWindow = function(data)
    owin(xrange = c(0, guessXMax(data)), yrange = c(0, guessYMax(data)))

guessXMax = function(data)
{
  w = 1392 / pixelsPerMicron
  ceiling(max(data[['Cell_X_Position']])/w)*w
}

guessYMax = function(data)
{
  h = 1040 / pixelsPerMicron
  ceiling(max(data[['Cell_Y_Position']])/h)*h
}

# Compute nearest distance to each phenotype for each cell in an inForm cell seg table
computeAllNearestDistance = function(cellTablePath=NULL, outPath=NULL)
{
  # Get the path to the cell seg table and check it
  if (is.null(cellTablePath))
    cellTablePath = file.choose()
  
  # Read the table
  cat('Reading', cellTablePath, '\n')
  data = readCellSegData(cellTablePath)
  
  # Compute the distances
  cat('Computing distances\n')
  data = cbind(data, findNearestDistance(data))
  
  if (is.null(outPath))
    outPath = sub('.txt', '_dist.txt', cellTablePath)
  cat('Writing', outPath, '\n')
  writeCellSegData(data, outPath)
}

# For each phenotype and each cell, find the minimum distance from
# the cell to a cell in the phenotype. 
# Return a data.frame containing 'Distance to <category>' column
# for each phenotype
findNearestDistance = function(data)
{
  d = makeDistanceMatrix(data)
  
  cats = levels(data[[phenotypeName]])
  result = lapply(cats, FUN=function(cat)
  {
    catCells = data[[phenotypeName]]==cat # Which cells are in the target category?
    catMins = apply(d[,catCells], 1, min)     # Find the minimum distance
    catMins
  })
  names(result) = paste('Distance to', cats)       # The names for the new columns

  data.frame(result, check.names=FALSE)
}
  
# Make a distance matrix - the distance from each cell to each other cell
makeDistanceMatrix = function(data)
{
  stopifnot('Cell_X_Position' %in% names(data),
            'Cell_Y_Position' %in% names(data),
            phenotypeName %in% names(data))
  
  d = as.matrix(dist(data[,c('Cell_X_Position', 'Cell_Y_Position')]))
  rownames(d) = colnames(d) = data[['Cell ID']]
  d
}

# Read a cell seg data file preserving names. Drop empty phenotypes and convert position to microns.
readCellSegData = function(cellTablePath)
{
  stopifnot(grepl('_cell_seg_data', cellTablePath))
  data = fread(cellTablePath, na.strings=c('NA', '#N/A'), check.names=FALSE, stringsAsFactors=TRUE)
  
  # Remove any cells with no phenotype
  droplevels(data[data[[phenotypeName]]!='',])

  # Convert to microns
  data[['Cell_X_Position']] = data[['Cell_X_Position']] / pixelsPerMicron
  data[['Cell_Y_Position']] = data[['Cell_Y_Position']] / pixelsPerMicron
  data
}

# Read a tissue seg summary data file preserving names. Convert size to microns.
readTissueSegSummaryData = function(tissueTablePath)
{
  stopifnot(grepl('_tissue_seg_data_summary', tissueTablePath))
  data = read.delim(tissueTablePath, na.strings=c('NA', '#N/A'), check.names=FALSE)
  
  # All we care about is category name and size
  data = data[,c('Tissue Category', 'Region Area (pixels)')]
  
  # Convert to microns
  data[['Region Area']] = data[['Region Area (pixels)']] / pixelsPerMicron^2
  data
}

# Read a cell seg summary data file preserving names. Convert size to microns.
readCellSegSummaryData = function(cellTablePath)
{
  stopifnot(grepl('_cell_seg_data_summary', cellTablePath))
  data = read.delim(cellTablePath, na.strings=c('NA', '#N/A'), check.names=FALSE)
  
  # All we care about is category name and size
  data = data[,c('Tissue Category', 'Tissue Category Area (pixels)')]
  
  # Convert to microns
  data[['Region Area']] = data[['Tissue Category Area (pixels)']] / pixelsPerMicron^2
  data
}
# Write a cell-seg-like data file
writeCellSegData = function(data, outPath)
  write.table(data, outPath, sep='\t', row.names=FALSE, na='#N/A')

theme_reports = function()
{
  theme_update(
  panel.background = element_rect(fill = "grey95", colour = NA),
  panel.grid.minor = element_line(colour = "grey99", size = 0.25),
  strip.background = element_rect(fill = "grey90", color="grey50"))
}

# Multiple plot function, to print a grid of ggplot objects
# From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
