### Trajectory analysis of Fibroblasts in PVAT and eWAT ####
# Author: Lea Mikkola

setwd("~/Desktop/TERVA/TERVA2DATA/")
library(Totem)
library(SingleCellExperiment)
library(cowplot)
library(dyndimred)
library(S4Vectors)
library(Seurat)
library(ggplot2)

load("TERVA_fibro_Totem.RData") #This has been created from the full integrated TERVA2_harmony data (seurat object), by just subsetting the fibroblasts.

Idents(Fibroblasts) <- "tissue_id"
PVAT <- subset(Fibroblasts, idents = "PVAT")
eWAT <- subset(Fibroblasts, idents = "eWAT")

#eWAT

set.seed(42)
counts <- eWAT[["RNA"]]@counts
logcounts <- eWAT[["RNA"]]@data
clusters <- list(as.character(eWAT@meta.data$celltype))
names(clusters[[1]]) <- colnames(eWAT)
fibrodata <- c(counts, logcounts, clusters)
y <- c("counts", "logcounts", "clusters")
names(fibrodata) <- paste0(y)
fibrodata["counts"] <- t(fibrodata$counts) 
fibrodata["logcounts"] <- t(fibrodata$logcounts) # It was suggested that I could use double square brackets and then not double transpose, but I tested this suggestion and dimred_mds was significantly slower (didn't wait it to go through, because the double-transpose is so much quicker)

sce <- SingleCellExperiment(assays = list(counts = t(fibrodata$counts), 
                                          logcounts = t(fibrodata$logcounts)))
sce <- PrepareTotem(sce)

sce <- RunDimRed(object = sce,
                 dim.red.method = "lmds",
                 dim.red.features = NULL,
                 dim.reduction.par.list = list(ndim=5))

dim_red <- dimred_mds(fibrodata$logcounts, ndim = 2)

sce <- RunClustering(sce,
                     k.range = 3:20,
                     min.cluster.size = 5,
                     N.clusterings=10000)

VizCellConnectivity(sce, viz.dim.red = dim_red)

sce <- SelectClusterings(sce, selection.method = 5,
                         selection.N.models = 10,
                         selection.stratified=FALSE,
                         prior.clustering = fibrodata$clusters)

trajnames <- ReturnTrajNames(sce) #Check names, note that these might differ if run with different package versions
clustering_name <- ReturnTrajNames(sce)[1] #If you don't agree with the first trajectory, state here a specific one. If there's only little difference between the estimated clustering results, then it should be fine to take the first one "arbitrarily".
clustering <- ReturnClustering(sce,clustering_name)
VizClustering(sce, clustering = clustering, viz.dim.red = dim_red) 
#VizClustering(sce, clustering = fibrodata$clusters, viz.dim.red = dim_red) #This can take a long time.

#ReturnMSTNetwork(sce,clustering.name = clustering_name) #check the MST network

sce <- RunSmoothing(sce) 
ReturnTrajNames(sce)
VizSmoothedTraj(sce,
                traj.names = ReturnTrajNames(sce),
                viz.dim.red = dim_red,plot.pseudotime = FALSE)

ReturnSmoothedTrajNetwork(sce,clustering.name = clustering_name)

#sce <- ChangeTrajRoot(sce,traj.name=clustering_name, root.cluster = 1) #Change the root cluster if necessary. This should always be checked so that the root makes sense based on the biology.
#VizSmoothedTraj(sce,
#                traj.names = "3.133",
#                viz.dim.red = dim_red,plot.pseudotime = FALSE)

#With pseudotime
VizSmoothedTraj(sce,
                traj.names = clustering_name,
                viz.dim.red = dim_red,plot.pseudotime = TRUE)

# PVAT ####

set.seed(42)
counts <- PVAT[["RNA"]]@counts
logcounts <- PVAT[["RNA"]]@data
clusters <- list(as.character(PVAT@meta.data$celltype))
names(clusters[[1]]) <- colnames(PVAT)
fibrodata <- c(counts, logcounts, clusters)
y <- c("counts", "logcounts", "clusters")
names(fibrodata) <- paste0(y)
fibrodata["counts"] <- t(fibrodata$counts) 
fibrodata["logcounts"] <- t(fibrodata$logcounts) 


sce <- SingleCellExperiment(assays = list(counts = t(fibrodata$counts), 
                                          logcounts = t(fibrodata$logcounts)))
sce <- PrepareTotem(sce)

sce <- RunDimRed(object = sce,
                 dim.red.method = "lmds",
                 dim.red.features = NULL,
                 dim.reduction.par.list = list(ndim=5))


dim_red <- dimred_mds(fibrodata$logcounts, ndim = 2)

sce <- RunClustering(sce,
                     k.range = 3:20,
                     min.cluster.size = 5,
                     N.clusterings=10000)

VizCellConnectivity(sce, viz.dim.red = dim_red)

sce <- SelectClusterings(sce, selection.method = 5,
                         selection.N.models = 10,
                         selection.stratified=FALSE,
                         prior.clustering = fibrodata$clusters)

trajnames <- ReturnTrajNames(sce) #Check names
clustering_name <- ReturnTrajNames(sce)[1] 
clustering <- ReturnClustering(sce,clustering_name)
VizClustering(sce, clustering = clustering, viz.dim.red = dim_red) 
#VizClustering(sce, clustering = fibrodata$clusters, viz.dim.red = dim_red) # This can take a long time.

#ReturnMSTNetwork(sce,clustering.name = "") #

sce <- RunSmoothing(sce) 
ReturnTrajNames(sce)
VizSmoothedTraj(sce,
                traj.names = ReturnTrajNames(sce),
                viz.dim.red = dim_red,plot.pseudotime = FALSE)

ReturnSmoothedTrajNetwork(sce,clustering.name = clustering_name)

sce <- ChangeTrajRoot(sce,traj.name=clustering_name, root.cluster = 3) #Root changed here
VizSmoothedTraj(sce,
                traj.names = clustering _name,
                viz.dim.red = dim_red,plot.pseudotime = FALSE)

#With pseudotime
VizSmoothedTraj(sce,
                traj.names = clustering_name,
                viz.dim.red = dim_red,plot.pseudotime = TRUE)


