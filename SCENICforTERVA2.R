setwd("/scratch/project_2005050/Rstats")
.libPaths(c("/projappl/project_2005050/project_rpackages_4.1.1", .libPaths()))

library(tidyverse)
library(GENIE3)
library(RcisTarget)
library(AUCell)
library(SCENIC)
library(SCopeLoomR)
library(SeuratDisk)
library(SingleCellExperiment)
library(Seurat)

#"https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
#"https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
 
# Across all samples ####
sce <- as.SingleCellExperiment((DietSeurat(TERVA2_harmony)))
cellInfo <- data.frame(seuratCluster=Idents(TERVA2_harmony))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
                           ))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9462	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

#### SCENIC BY SAMPLE ####

# Read in the data

LDAO <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDAOafterQC2.h5Seurat") #ReQC
PLAO <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLAOafterQC2.h5Seurat") #ReQC
LDeWAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/LD_eWATafterQC.h5Seurat")
PLeWAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLeWATafterQC2.h5Seurat") #ReQC
LDPVAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDPVATafterQC.h5Seurat")
PLPVAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLPVATafterQC2.h5Seurat") #ReQC
LDSpleen <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDSpleenafterQC.h5Seurat")
PLSpleen <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLSpleenafterQC.h5Seurat")

#Late disease aorta ####

sce <- as.SingleCellExperiment((DietSeurat(LDAO)))
cellInfo <- data.frame(seuratCluster=Idents(LDAO))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9731	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

scenicOptions@settings$seed <- 1234
# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

nPcs <- c(5,15,50)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# Using only "high-confidence" regulons (normally similar)
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds

# Prelesion aorta ####

sce <- as.SingleCellExperiment((DietSeurat(PLAO)))
cellInfo <- data.frame(seuratCluster=Idents(PLAO))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9462	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

scenicOptions@settings$seed <- 1234
# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Late disease eWAT ####

sce <- as.SingleCellExperiment((DietSeurat(LDeWAT)))
cellInfo <- data.frame(seuratCluster=Idents(LDeWAT))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9462	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

scenicOptions@settings$seed <- 1234
# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Prelesion eWAT ####

sce <- as.SingleCellExperiment((DietSeurat(TERVA2_harmony)))
cellInfo <- data.frame(seuratCluster=Idents(TERVA2_harmony))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9462	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

scenicOptions@settings$seed <- 1234
# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Late disease PVAT ####

sce <- as.SingleCellExperiment((DietSeurat(LDPVAT)))
cellInfo <- data.frame(seuratCluster=Idents(LDPVAT))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9462	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

scenicOptions@settings$seed <- 1234
# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Prelesion PVAT ####

sce <- as.SingleCellExperiment((DietSeurat(PLPVAT)))
cellInfo <- data.frame(seuratCluster=Idents(PLPVAT))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9462	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

scenicOptions@settings$seed <- 1234
# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Late disease Spleen ####

sce <- as.SingleCellExperiment((DietSeurat(LDSpleen)))
cellInfo <- data.frame(seuratCluster=Idents(LDSpleen))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9462	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

scenicOptions@settings$seed <- 1234
# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

# Prelesion Spleen ####

sce <- as.SingleCellExperiment((DietSeurat(PLSpleen)))
cellInfo <- data.frame(seuratCluster=Idents(PLSpleen))
exprMat <- counts(sce)
loom <- build_loom("/scratch/project_2005050/Rstats/scenic.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

loom <- open_loom("/scratch/project_2005050/Rstats/scenic.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)

colVars <- list(CellType=c("B cells"="forestgreen", 
                           "Fibroblasts"="darkorange", 
                           "Macrophages"="magenta4", 
                           "Cd8+ Tcells"="hotpink", 
                           "Cd4+ Tcells"="red3", 
                           "Macrophages activated"="skyblue", 
                           "Classical and non-classical Monocytes"="darkblue",
                           "Intermediate monocytes "="darkmagenta",
                           "Lgals3+ Macrophages"="cyan2",
                           "NK"="darkviolet",
                           "EC's"="deeppink4",
                           "Cd4+ Foxp3+ Tregs"="darkgoldenrod1",
                           "VSMC's"="darksalmon",
                           "Cd8+ Ccl5+ Teffs"="chartreuse4",
                           "Hdc+ Cpa3+ MAST cells"="darkorchid3",
                           "Plasma cells"="aquamarine3",
                           "Conv DC1"="antiquewhite3",
                           "Conv DC2"="cornsilk3",
                           "ILC's"="deepskyblue3",
                           "Granulocytes"="darkcyan",
                           "Fibroblasts activated"="bisque3",
                           "Pi16+ Fibroblasts"="coral3",
                           "Rgs5+ EC's"="brown2",
                           "Mgp+ Fibroblasts"="darkseagreen4"
))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$seuratCluster)]
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

scenicOptions <- initializeScenic(org="mgi", dbDir="/scratch/project_2005050/Rstats", datasetTitle = "SCENIC TERVA2", nCores=10)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

#9462	genes available in RcisTarget database

interestingGenes <- c("Ly6a", "Trem2", "Slc7a7")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

SCE_corr <- runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered <- log2(exprMat_filtered+1) #Run because the count matrix originates from the raw count data of the seurat object

scenicOptions@settings$seed <- 1234
# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)

