#--------------------------------------------------------------------------
# NAME: scRNAseq_analysis_demo.R
#
# INPUT FILES:
# - Read count file & R source
#
# ARGUMENTS: "User's input parameters"
# - dirIn: Input directory where input data and R source are located
# - fin  : Input file with read counts
#
# REQUIREMENT:
# - R-library: Seurat, dplyr, topGO, tidyverse, org.Hs.eg.db
#--------------------------------------------------------------------------

#------------------------------------------------------------
# User's input parameters
#------------------------------------------------------------
dirIn  <- "/Input_data_PATH/" 
fin    <- "RNAseq_profile_organoid8M.txt" 
#------------------------------------------------------------

# Install R packages
source("http://bioconductor.org/biocLite.R")
biocLite(c("Seurat","dplyr","topGo","tidyverse","org.Hs.eg.db","org.Mm.eg.db"))

# load library
library(Seurat)
library(dplyr)
library(topGO)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# set working dir
setwd(dirIn)

# sources
source("Utilities_sc.R")


# read data
exprTable  <- read.table(fin,sep="\t",header=T,row.names=1,quote="")
dim(exprTable)

# create a Seurat object with the raw (non-normalized data) Keeping all genes expressed in >= 68 cells (~0.5% of the data),
# and Keep all cells with at least 300 detected genes
geneCutoff <- dim(exprTable)[2] * 0.05
dataObj <- CreateSeuratObject(raw.data = exprTable, min.cells = geneCutoff, min.genes = 300, project = "Organoid")
dataObj

# calculate mitochondria percentage
mito.genes <- grep("^MT-", rownames(dataObj@data), value = T)
percent.mito <- Matrix::colSums(dataObj@raw.data[mito.genes, ])/Matrix::colSums(dataObj@raw.data)
dataObj      <- AddMetaData(object = dataObj, metadata = percent.mito, col.name = "percent.mito")
head(dataObj@meta.data$percent.mito)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
VlnPlot(object = dataObj, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = dataObj, gene1 = "nUMI", gene2 = "nGene",col="black")

# Filter out cells with fewer than 500 genes detected and percent.mito higher than 0.1 (10%)
dataObj <- FilterCells(object = dataObj, subset.names = c("nGene", "percent.mito"), 
                  low.thresholds = c(1000, -Inf), high.thresholds = c(Inf, 0.2))
dataObj

# Normalization using scale factor
dataObj   <- NormalizeData(object= dataObj, normalization.method = "LogNormalize", scale.factor = 10000)

# identifying variable genes across the single cells
dataObj   <- FindVariableGenes(object= dataObj)

# Z-scoring the data and removing unwanted sources of variation
dataObj   <- ScaleData(dataObj)

# Dimensionality Reduction and visualization

# finding 1000 hvg genes, user can use other number of genes
hvgs <- head(rownames(dataObj@hvg.info), "1000")
head(hvgs)

# run PCA
dataObj <- RunPCA(dataObj, genes.use=hvgs)

# draw PCA plot
PCAPlot(object = dataObj, dim.1 = 1, dim.2 = 2)

# generate the heatmap of PCA output
PCHeatmap(object = dataObj, pc.use = 1:12, cells.use = 100, do.balanced = TRUE, label.columns = FALSE)

# determine the number of PCs for t-SNE plot
dataObj <- JackStraw(object = dataObj)
JackStrawPlot(object = dataObj, PCs = 1:18)
PCElbowPlot(object = dataObj)

# find clusters of the cells
dataObj <- FindClusters(object = dataObj,  resolution = 0.6, reduction.type = "pca", dims.use = 1:10)

# run t-sne
dataObj <- RunTSNE(object = dataObj, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = dataObj, pt.size = 2, do.label = T, label.size = 10)

# find differentially expressed gene markers for every cluster compared to all remaining cells
dataObj.markers <- FindAllMarkers(object = dataObj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

# select top 20 genes of each cluster
outTop20 <- dataObj.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

# gererate heatmap
DoHeatmap(object = dataObj, genes.use = outTop20$gene, slim.col.label = TRUE, remove.key = TRUE, group.cex=10, cex.row=8)

# test marker expression
cellMarker <- c("ARR3","PDE6H", "GUCA1C","SAG","GNAT1","NR2E3","NFIB","DKK3","RLBP1")
FeaturePlot(object = dataObj, features.plot = cellMarker, min.cutoff="q9", cols.use= c("lightgrey", "blue"), pt.size=0.5)
VlnPlot(object = dataObj, features.plot = cellMarker, use.raw = TRUE, y.log = TRUE)

#Assigning cell type identity to clusters
# Identify differential expressed genes across conditions
currentClIds <- c(0, 1, 2, 3)
newClIds <- c("Cone", "Rod", "MG", "Mixed")
dataObj@ident <- plyr::mapvalues(x = dataObj@ident, from = currentClIds, to = newClIds)
TSNEPlot(object = dataObj, do.label = TRUE, pt.size = 2, label.size=7)

# perform GO analysis
# select cluster number
nCluster <- 0
outTop50 <- dataObj.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
DEGTable <- outTop50[outTop50$p_val<0.01 & outTop50$cluster==nCluster,]
DEGs <- DEGTable$gene
head(DEGs)

# run GO analysis
tGOterms   <- topGOterms(fg.genes = DEGs, bg.genes = rownames(dataObj@data), organism = "Human") # option: Human, Mouse
outGOterms <- tGOterms$res.table
head(outGOterms)

# generate a graph for significant GO terms
outGO   <- outGOterms[,c(2,6)]
GO_term <- outGO$Term
P_val   <- -log10(as.numeric(outGO$pval))
bplot <- outGO %>% mutate(GO_term = fct_reorder(GO_term, P_val)) %>% ggplot( aes(x=GO_term, y=P_val)) 
bplot + geom_bar(stat="identity", fill="dark blue") +  coord_flip()+ labs(y="- Log10 (P-val)", x=NULL)+ ggtitle("GO analysis for top 50 DEGs of cluster_0")

# SAve a object
saveRDS("dataObj.rds")
