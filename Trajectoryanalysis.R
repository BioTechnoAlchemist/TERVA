### Trajectory analysis of Fibroblasts in PVAT and eWAT ####

#remotes::install_github("dynverse/dynfeature")
#remotes::install_github("dynverse/dynplot")
#remotes::install_github("elolab/Totem")

suppressMessages(library(Totem))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(cowplot))
suppressMessages(library(dyndimred))
suppressMessages(library(S4Vectors))

load("TERVA_fibro_Totem.RData")

Idents(Fibroblasts) <- "tissue_id"
PVAT <- subset(Fibroblasts, idents = "PVAT")
eWAT <- subset(Fibroblasts, idents = "eWAT")

set.seed(42)
counts <- Fibroblasts[["RNA"]]@counts
logcounts <- Fibroblasts[["RNA"]]@data
clusters <- list(as.character(Fibroblasts@meta.data$celltype))
names(clusters[[1]]) <- colnames(Fibroblasts)
fibrodata <- c(counts, logcounts, clusters)
y <- c("counts", "logcounts", "clusters")
names(fibrodata) <- paste0(y)
fibrodata["counts"] <- t(fibrodata$counts)
fibrodata["logcounts"] <- t(fibrodata$logcounts)


sce <- SingleCellExperiment(assays = list(counts = t(fibrodata$counts),
                                          logcounts = t(fibrodata$logcounts),
                                          colData = DataFrame()))
sce <- PrepareTotem(sce)

sce <- RunDimRed(object = sce,
                 dim.red.method = "lmds",
                 dim.red.features = NULL,
                 dim.reduction.par.list = list(ndim=5))


#Alternative

set.seed(42)
counts <- Fibroblasts[["RNA"]]@counts
logcounts <- Fibroblasts[["RNA"]]@data
clusters <- list(as.character(Fibroblasts@meta.data$celltype))
names(clusters[[1]]) <- colnames(Fibroblasts)
fibrodata <- c(counts, logcounts, clusters)
y <- c("counts", "logcounts", "clusters")
names(fibrodata) <- paste0(y)
fibrodata["counts"] <- t(fibrodata$counts) # this is fine, but you don't need to transpose if you use double square bracket like this: fibrodata[["counts"]] <- fibrodata$counts
fibrodata["logcounts"] <- t(fibrodata$logcounts) # same here


sce <- SingleCellExperiment(assays = list(counts = t(fibrodata$counts), # this is fine, but if you don't transpose above and use square brackets, you don't need to transpose here
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
VizClustering(sce, clustering = clustering, viz.dim.red = dim_red) # there isn't a 'clustering.names' parameter
VizClustering(sce, clustering = fibrodata$clusters, viz.dim.red = dim_red) 

#eWAT

set.seed(42)
counts <- eWAT[["RNA"]]@counts
logcounts <- eWAT[["RNA"]]@data
clusters <- list(as.character(eWAT@meta.data$celltype))
names(clusters[[1]]) <- colnames(eWAT)
fibrodata <- c(counts, logcounts, clusters)
y <- c("counts", "logcounts", "clusters")
names(fibrodata) <- paste0(y)
fibrodata["counts"] <- t(fibrodata$counts) # this is fine, but you don't need to transpose if you use double square bracket like this: fibrodata[["counts"]] <- fibrodata$counts
fibrodata["logcounts"] <- t(fibrodata$logcounts) # same here ### OWN COMMENT: I actually tested this suggestion and dimred_mds was significantly slower (didn't wait it to go through, because the double-transpose is so much quicker)


sce <- SingleCellExperiment(assays = list(counts = t(fibrodata$counts), # this is fine, but if you don't transpose above and use square brackets, you don't need to transpose here
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
VizClustering(sce, clustering = clustering, viz.dim.red = dim_red) # there isn't a 'clustering.names' parameter
VizClustering(sce, clustering = fibrodata$clusters, viz.dim.red = dim_red) 

#ReturnMSTNetwork(sce,clustering.name = "") #check the MST network

sce <- RunSmoothing(sce) 
ReturnTrajNames(sce)
VizSmoothedTraj(sce,
                traj.names = ReturnTrajNames(sce),
                viz.dim.red = dim_red,plot.pseudotime = FALSE)

ReturnSmoothedTrajNetwork(sce,clustering.name = "3.50")

#sce <- ChangeTrajRoot(sce,traj.name="3.133", root.cluster = 1) #Change the root cluster if necessary
#VizSmoothedTraj(sce,
#                traj.names = "3.133",
#                viz.dim.red = dim_red,plot.pseudotime = FALSE)

#OR with pseudotime
VizSmoothedTraj(sce,
                traj.names = "3.50",
                viz.dim.red = dim_red,plot.pseudotime = TRUE)

#Gene expression

VizFeatureExpression(sce,traj.name = "3.50",
                     feature.names = "Cd248",
                     viz.dim.red = dim_red,
                     plot.traj = TRUE)

