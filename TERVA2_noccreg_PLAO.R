#### TERVA2 scRNA-seq experiment 2021 ####
# Author: Lea Mikkola

setwd("~/Desktop/TERVA/TERVA2DATA/")
library(tidyverse)
library(Seurat)
library(patchwork)
library(SoupX)
library(sctransform)
library(viridis)

#### Ambient RNA removal ####
# SoupX 1.5.2

#Before using load10x, move raw_feature_bc_matrix folder to the per_sample_outs counts-folder, and change the name of sample_feature_bc_matrix to filtered_feature_bc_matrix. Otherwise load10x will end up in an error.
sc <- load10X("~/path", includeFeatures = "Gene Expression") #This will only include the gene expression data. Include ab data later.
sc

#Channel with 32285 genes and 10609 cells (PL AO)

dd <- sc$metaData[colnames(sc$toc), ]
ggplot(dd, aes(tSNE1, tSNE2), group=factor(clusters)) + 
  geom_point(aes(colour = factor(clusters))) +
  scale_color_viridis(discrete = T, option = "C") +
  scale_fill_viridis(discrete = T) +
  theme(legend.position = "bottom") 



sc$Igkc <- as.data.frame(sc$toc["Igkc", ])
p1<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Igkc > 0)) #Could work as a soup marker as is "expressed" all over the place
sc$Ighg1 <- as.data.frame(sc$toc["Ighg1", ])
p2<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Ighg1 > 0)) 
sc$Ighm <- as.data.frame(sc$toc["Ighm", ])
p3<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Ighm > 0))
sc$Iglc2 <- as.data.frame(sc$toc["Iglc2", ])
p4<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Iglc2 > 0))
sc$Iglc3 <- as.data.frame(sc$toc["Iglc3", ])
p5<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Iglc3 > 0))
sc$Iglc1 <- as.data.frame(sc$toc["Iglc1", ])
p6<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Iglc1 > 0))
sc$Ighe <- as.data.frame(sc$toc["Ighe", ])
p7<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Ighe > 0))
sc$Ighg3 <- as.data.frame(sc$toc["Ighg3", ])
p8<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Ighg3 > 0))
sc$Ighd <- as.data.frame(sc$toc["Ighd", ])
p9<-ggplot(dd, aes(tSNE1, tSNE2)) + geom_point(aes(colour = sc$Ighd > 0))

p1+p2+p3+p4+p5+p6+p7+p8+p9

#DR = sc$metaData[,sc$DR]
#DR = as.data.frame(DR)
genes<- as.vector(c("Igkc","Ighg1","Ighm","Iglc2","Iglc3","Iglc1","Ighe","Ighg3","Ighd", "Igkv3-4", "Igkv1-117", "Igkv3-5", "Jchain")) 
#allgenes <- as.vector(sc[["toc"]]@Dimnames[[1]])
#plotMarkerMap(sc, geneSet = genes)

autoEstCont(sc)

#2260 genes passed tf-idf cut-off and 448 soup quantile filter.  Taking the top 100.
#Using 653 independent estimates of rho.
#Estimated global rho of 0.02
#Channel with 32285 genes and 10609 cells

#head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution(sc)

useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = genes)) #Use if AutoEstCont is not stringent enough (< 2%)
plotMarkerMap(sc, geneSet = genes, useToEst = useToEst)
sc <- calculateContaminationFraction(sc, list(IG = genes), useToEst = useToEst) #Estimated global contamination fraction of 5.76%


out <- adjustCounts(sc)
#cntSoggy = rowSums(sc$toc > 0)
#cntStrained = rowSums(out > 0)
#mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
#mostZeroed #What decreased most
#tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20) #Which genes had a quantitative difference
#plotChangeMap(sc, out, "Igkc") #Check some genes of interest for changes in expression after correction
#plotChangeMap(sc, out, "Ighd")
#plotChangeMap(sc, out, "Cd3d")
#plotChangeMap(sc, out, "Cd74")
#plotChangeMap(sc, out, "Cd79a")
#plotChangeMap(sc, out, "Fabp5")
#plotChangeMap(sc, out, "Spp1")


srat <- CreateSeuratObject(out)

#### ADD antibody data ####
tenx <- Read10X("/Users/limikk/Desktop/TERVA/TERVA2DATA/PL_AO_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_AO_GEX_Ab_VDJ_noIntr_SC5PPE/count/filtered_feature_bc_matrix/")
srat[["ADT"]] <- CreateAssayObject(counts = tenx[["Antibody Capture"]])


#### DoubletFinder ####

library(Matrix)
library(fields)
library(KernSmooth)
library(ROCR)
library(parallel)
library(DoubletFinder)

# First some rough filtering, data normalization and scaling - must be done before DoubletFinder

srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-") #Save mito-% as new metadata
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]") #Save ribo-% as new metadata
VlnPlot(srat, features= c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol= 4) #Data before pruning based on no of Features and mt percentage (after SoupX)
sub <- subset(srat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & percent.rb > 0.05) # Here we prune the data with the given thresholds to get rid of some noise.
VlnPlot(sub, features= c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol= 4) #Data after pruning

srat <- NormalizeData(sub)
srat <- FindVariableFeatures(srat)
#s.genes <- str_to_title(cc.genes$s.genes) #If this gives an error, reinstall package stringr 
#g2m.genes <- str_to_title(cc.genes$g2m.genes)
srat <- ScaleData(srat, features = rownames(srat)) #Without features = rownames(srat) only variable features are scaled which may affect the downstream analyses. Scale all!
srat <- RunPCA(srat, verbose = F)
ElbowPlot(srat)

srat <- RunUMAP(srat, dims = 1:20, verbose = F)
srat <- FindNeighbors(srat, dims = 1:20, verbose = F)
srat <- FindClusters(srat, resolution = 0.6, verbose = F)
DimPlot(srat, reduction = "umap", pt.size=1, group.by = "RNA_snn_res.0.6", label= T) 
#VlnPlot(srat, features = c("Cd3d", "Cd79a"))
#DotPlot(srat, features = c("Cd3d", "Cd79a")) + RotatedAxis()
#DoHeatmap(subset(srat, downsample = 100), features = c("Cd3d", "Cd79a", "Acta2", "Spp1", "Jchain", "Dcn"), size = 3)

#No cell cycle regression done for this data
# Continue to DoubletFinder

sweep.res.list <- paramSweep_v3(srat, PCs = 1:20, sct = F) #sct = TRUE when SCTransform used. Results can greatly vary in different data sets with either sct or basic normalization.
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) #GT = FALSE when no ground truth used
bcmvn <- find.pK(sweep.stats)
pK=as.numeric(as.character(bcmvn$pK)) #This script is from WeiZhu1998 from https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/62
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2) 
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.076*nrow(srat@meta.data)) #0.054 is based on the 10X estimate of multiplet formation rate dependent on cell loading. Given that there appears to be more cells than what were assumed to be loaded, the loading number was probably higher. Therefore I use a bit higher doublet rate estimate here.
srat <- doubletFinder_v3(srat, PCs = 1:20, pN = 0.25, pK = 0.17, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
#srat_reuse <- doubletFinder_v3(srat, PCs = 1:20, pN = 0.25, pK = 0.22, nExp = nExp_poi, reuse.pANN = "pANN_0.25_0.22_120", sct = TRUE)

DF.name = colnames(srat@meta.data)[grepl("DF.classification", colnames(srat@meta.data))]
VlnPlot(srat, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

plot1 <- DimPlot(srat, group.by="DF.classifications_0.25_0.17_630", reduction="umap", pt.size=1, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB"))
plot2 <- DimPlot(srat, reduction = "umap", pt.size=1, group.by = "RNA_snn_res.0.6")

plot1 + plot2 

srat_DF <- subset(srat, subset = DF.classifications_0.25_0.17_630 == "Singlet") #Creates a new Seurat Object from the singlets.

plot3 <- DimPlot(srat_DF, group.by="DF.classifications_0.25_0.17_630", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB"))
plot4 <- DimPlot(srat_DF, reduction = "umap", pt.size=0.5, label = T, group.by = "RNA_snn_res.0.6")

plot3 + plot4 # looks same as before so ok



#### Normalize and scale data again for DE analysis ####
#Do this again after creating the new SeuratObject because it has changed from the original data.

srat_N <- NormalizeData(srat_DF, verbose = F)
srat_N <- FindVariableFeatures(srat_N, selection.method = "vst", nfeatures = 3000)
#top10 <- head(VariableFeatures(srat_N), 10)
#plot5 <- VariableFeaturePlot(srat_N)
#plot6 <- LabelPoints(plot = plot5, points = top10, repel = TRUE)
#plot5 + plot6
srat_N <- ScaleData(srat_N, verbose = T, features = rownames(srat_N))
srat_N <- RunPCA(srat_N, verbose = F)
ElbowPlot(srat_N)

srat_N <- RunUMAP(srat_N, dims = 1:20, verbose = F)
srat_N <- FindNeighbors(srat_N, dims = 1:20, verbose = F)
srat_N <- FindClusters(srat_N, resolution = 0.7, verbose = F)
plot7 <- DimPlot(srat_N, reduction = "umap", pt.size=0.5, label = T, group.by = "RNA_snn_res.0.7")

#plot some genes of interest
FeaturePlot(srat_N, features = c("Myh11", "Cnn1", "Tcf21", "C1qa", "Lyz2", "H2-Ab1", "Vwf", "Flt1", "Kdr", "Ltb", "Pdgfra", "Lum", "Pecam1", "Cdh5", "Cd79a", "Cd3d", "Cnp"), pt.size = 0.2,
            ncol = 4)
FeaturePlot(srat_N, features = c("Il7r", "Icos", "Kit","Rag1"), pt.size = 0.2,
            ncol = 3) #Il7r = CD127, Icos = CD278, Kit = CD117
VlnPlot(srat_N, features = c("Il7r", "Icos", "Kit"), pt.size = 0.2, ncol = 3)
DotPlot(srat_N, features = c("Il7r", "Icos", "Kit")) + RotatedAxis()

VlnPlot(srat_N, features = c("Myh11", "Cnn1", "C1qa", "Lyz2", "Flt1", "Pdgfra", "Lum", "Pecam1", "Cdh5", "Cd79a", "Cd3d"), pt.size = 0.2, ncol = 3)
plot9 <- DotPlot(srat_N, features = c("Myh11", "Cnn1", "C1qa", "Lyz2", "Flt1", "Pdgfra", "Lum", "Pecam1", "Cdh5", "Cd79a", "Cd3d")) + RotatedAxis()
plot10 <- DoHeatmap(subset(srat_N, downsample = 100), features = c("Myh11", "Cnn1", "C1qa", "Lyz2", "Flt1", "Pdgfra", "Lum", "Pecam1", "Cdh5", "Cd79a", "Cd3d"), size = 3)

plot7 + plot9 + plot10

#### SAVE H5-object to use in subsequent analyses ####

library(SeuratDisk)
library(SeuratData)

SaveH5Seurat(object = srat_N, overwrite = T, verbose = T)

