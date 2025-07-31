#### TERVA2 data integration with harmony ####
setwd("/scratch/project_2005050/Rstats")
.libPaths(c("/projappl/project_2005050/project_rpackages_4.1.1/", .libPaths()))
library(tidyverse) 
library(patchwork) 
library(cowplot)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(factoextra)
library(clustree)
library(harmony)
library(colorBlindness)
library(EnhancedVolcano)
library(SeuratDisk, lib.loc = "/appl/soft/math/r-env/421/421-rpackages")
library(Seurat, lib.loc = "/appl/soft/math/r-env/421/421-rpackages") #v4.1.1 (SeuratObject v4.1.0)

#### Load in H5 data ####

OBAO <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDAOafterQC.h5Seurat") 
NOBAO <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLAOafterQC.h5Seurat") 
OBeWAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/LD_eWATafterQC.h5Seurat")
NOBeWAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLeWATafterQC.h5Seurat") 
OBPVAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDPVATafterQC.h5Seurat")
NOBPVAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLPVATafterQC.h5Seurat")
OBSpleen <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDSpleenafterQC.h5Seurat")
NOBSpleen <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLSpleenafterQC.h5Seurat")

New_idents <- rep("OBAO", times = length(OBAO$orig.ident))
OBAO@meta.data$Sample <- as.factor(New_idents)
New_idents <- rep("NOBAO", times = length(NOBAO$orig.ident))
NOBAO@meta.data$Sample <- as.factor(New_idents)

New_idents <- rep("OBeWAT", times = length(OBeWAT$orig.ident))
OBeWAT@meta.data$Sample <- as.factor(New_idents)
New_idents <- rep("NOBeWAT", times = length(NOBeWAT$orig.ident))
NOBeWAT@meta.data$Sample <- as.factor(New_idents)

New_idents <- rep("OBPVAT", times = length(OBPVAT$orig.ident))
OBPVAT@meta.data$Sample <- as.factor(New_idents)
New_idents <- rep("NOBPVAT", times = length(NOBPVAT$orig.ident))
NOBPVAT@meta.data$Sample <- as.factor(New_idents)

New_idents <- rep("OBSpleen", times = length(OBSpleen$orig.ident))
OBSpleen@meta.data$Sample <- as.factor(New_idents)
New_idents <- rep("NOBSpleen", times = length(NOBSpleen$orig.ident))
NOBSpleen@meta.data$Sample <- as.factor(New_idents)

#### Merge the separate data sets ####
# merge.data=T 
TERVA2 <-  merge(OBAO, y = c(NOBAO, OBeWAT, NOBeWAT, OBPVAT, NOBPVAT, OBSpleen, NOBSpleen), add.cell.ids = c("OBAO", "NOBAO", "OBeWAT", "NOBeWAT", "OBPVAT", "NOBPVAT", "OBSpleen", "NOBSpleen"), project = "TERVA2", merge.data = T)

TERVA2

table(TERVA2$Sample)

#NOBAO   NOBPVAT NOBSpleen   NOBeWAT      OBAO    OBPVAT  OBSpleen    OBeWAT 
#7830      5166      5264      6420      3048      7396      4234      8076 

#### Normalize, scale etc. RNA assay ####
# See https://github.com/immunogenomics/harmony/issues/41 for why I decided to use NormalizeData instead of SCTransform.
# See https://portals.broadinstitute.org/harmony/SeuratV3.html. Before running Harmony, make a Seurat object and following the standard pipeline through PCA.
#IMPORTANT DIFFERENCE: In the Seurat integration tutorial, you need to define a Seurat object for each dataset. With Harmony integration, create only one Seurat object with all cells.

remove("OBAO","OBeWAT","NOBeWAT","OBSpleen", "NOBAO", "OBPVAT", "NOBPVAT", "NOBSpleen", "New_idents")

# Remove some genes that commonly cause technical noise

genes <- GetAssayData(TERVA2, assay = "RNA")
kept.genes <- genes[-(which(rownames(genes) %in% c('Gm42418','AY036118'))),]
TERVA2[["RNA"]] <- CreateAssayObject(counts = kept.genes)
rm(genes, kept.genes)

TERVA2 <- NormalizeData(TERVA2) #We need to re-normalize because we removed the two genes above.
TERVA2 <- ScaleData(TERVA2, features = rownames(TERVA2), verbose = FALSE)
TERVA2 <- FindVariableFeatures(TERVA2, selection.method = "vst", nfeatures = 3000)
#saveRDS(TERVA2, "TERVA2_merged_normalized_scaled.rds")
TERVA2 <- RunPCA(TERVA2, verbose = FALSE)
#Create mat and pca objects for calculating the total variance:
mat <- GetAssayData(TERVA2, assay = "RNA", slot = "scale.data")
pca <- TERVA2[["pca"]]
# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
props <- eigValues / sum(eigValues)
plot(props, ylab="Proportion of variance", xlab="Principal Component",main = "Proportion of variance explained by the PC's")
rm(mat, pca)
gc()

TERVA2 <- RunUMAP(TERVA2, dims = 1:20)
TERVA2 <- FindNeighbors(TERVA2, dims = 1:20)
TERVA2 <- FindClusters(TERVA2, resolution = 0.5, random.seed = 42) 

# Check how the non-integrated data looks after merging and normalization.

cbcols <- c("#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", "#490092", "#006ddb", "#b66dff")
clustersbysample <- DimPlot(TERVA2, reduction = "umap", shuffle = T, repel = T, pt.size = 2, group.by = "Sample", cols = cbcols) +
  labs(title = "UMAP projection of the merged data") +
  theme(plot.title = element_text(size=30, hjust= 0.5),
    legend.text = element_text(size = 30),
    axis.title = element_text(size = 30))

clustersbysample 

#### Normalization of ADT data ####

abs <- GetAssayData(TERVA2, assay = "ADT")
kept.abs <- abs[-(which(rownames(abs) %in% c('CD11a'))),] # Remove one antibody that BioLegend reported to be of poor quality
TERVA2[["ADT"]] <- CreateAssayObject(counts = kept.abs)
rm(abs, kept.abs)

TERVA2 <- NormalizeData(TERVA2, assay = "ADT", normalization.method = "CLR")
TERVA2 <- ScaleData(TERVA2, assay = "ADT")
TERVA2 <- RunPCA(TERVA2, assay = "ADT", reduction.name = "apca", features = rownames(TERVA2@assays$ADT))

#### Running Harmony ####

TERVA2_harmony <- TERVA2 %>% RunHarmony(group.by.vars = "Sample", plot_convergence = T)

TERVA2_harmony <- TERVA2_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)

rm(TERVA2)

#### Check the Ab data ####

# Antibody data had varying success rates for the different samples.

ADTs <- as.list(TERVA2_harmony@assays$ADT@counts@Dimnames)
ADTnames <- ADTs[[1]]

for(i in 1:length(ADTnames)){ #Print and save (137) ADT images to a specific folder
  FeaturePlot(TERVA2_harmony, ADTnames[i])
  ggsave(paste0("/scratch/project_2005050/Rstats/Manuscriptimages/ABimages/", i, ".jpg"))
  while (!is.null(dev.list()))  dev.off()
}


#### Clustree ####

resolution.range <- seq(from = 0, to = 1.5, by = 0.1)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = resolution.range, random.seed = 42)

plot_cls1 <- clustree(TERVA2_harmony)
plot_cls1 +
  labs(title = "Cluster tree with resolutions from 0 to 1.5")

#### Check how everything looks at this point, some basic plots ####

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "Sample")
clbysample_harmony <- DimPlot(TERVA2_harmony,reduction = "umap", shuffle = T, cols = cbcols) + labs(title = "Samples after integration (Harmony)")
clbysample_harmony

nb.cols <- 26
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.0.5")
clbycl_harmony_0.5 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 0.5")
nb.cols <- 29
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.0.6")
clbycl_harmony_0.6 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 0.6")
nb.cols <- 33
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.0.7")
clbycl_harmony_0.7 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 0.7")
nb.cols <- 36
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.0.8")
clbycl_harmony_0.8 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 0.8")
nb.cols <- 40
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.0.9")
clbycl_harmony_0.9 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 0.9")
nb.cols <- 41
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.1")
clbycl_harmony_1 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 1")
nb.cols <- 44
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.1.1")
clbycl_harmony_1.1 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 1.1")
nb.cols <- 48
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.1.2")
clbycl_harmony_1.2 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 1.2")
nb.cols <- 49
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.1.3")
clbycl_harmony_1.3 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors) +
  labs(title = "Clustering with resolution 1.3")

(clbycl_harmony_0.5 + clbycl_harmony_0.6 + clbycl_harmony_0.7) / (clbycl_harmony_0.8 + clbycl_harmony_0.9 + clbycl_harmony_1) / (clbycl_harmony_1.1 + clbycl_harmony_1.2 + clbycl_harmony_1.3) 

clbysample_harmony + clbycl_harmony_1

# Feature plots of B and T cell markers (to inspect the double positive cluster)

TB <- FeaturePlot(TERVA2_harmony, features= c("Cd79a", "Cd19", "Cd3d", "Cd3e", "Cd3g"))
TB

TBcells <- WhichCells(TERVA2_harmony, expression = Cd79a > 0  &  Cd19 > 0 & Cd3d > 0 & Cd3e > 0 & Cd3g > 0)
DimPlot(TERVA2_harmony, cells.highlight = TBcells, order = T)
length(TBcells) / length(TERVA2_harmony$RNA_snn_res.1) * 100 #2.72 % of all cells
TBsub <- subset(TERVA2_harmony, subset = RNA_snn_res.1 == "12")
length(TBcells) / length(TBsub$nCount_RNA) * 100 # 83.09 % of cluster 12

DefaultAssay(TERVA2_harmony) <- "ADT"
TBsub <- subset(TERVA2_harmony, subset = RNA_snn_res.1 == "12")
TBcells <- WhichCells(TBsub, expression = CD19 > 0 & CD3 > 0)
length(TBcells) / length(TBsub$nCount_ADT) * 100 # 88.75 % of cluster 12

DefaultAssay(TERVA2_harmony) <- "RNA"
TB_avgexpr <- as.data.frame(AverageExpression(TERVA2_harmony, assays = "RNA", features= c("Cd79a", "Cd19", "Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd8b1", "Cd4", "Tcf7"), return.seurat = FALSE, slot = "data"))
TB_avgexpr <- TB_avgexpr %>% dplyr::select(RNA.0, RNA.3, RNA.5, RNA.12, RNA.20, RNA.22)
VlnPlot(TERVA2_harmony, slot = "counts", idents = c("0", "3", "5", "12", "20", "22"), features= c("Cd79a", "Cd19", "Cd3d", "Cd3e", "Cd3g", "Cd8a", "Cd8b1", "Cd4", "Tcf7"), cols = c("goldenrod","green","blue", "gold", "purple", "orange"))
VlnPlot(TERVA2_harmony, assay = "ADT", slot = "counts", idents = c("0", "3", "5", "12", "20", "22"), features= c("CD19", "CD3", "TCRbetachain", "TCRgamma-delta", "TCRVgamma1.1-Cr4", "TCRVgamma2", "TCRVgamma3","TCRVbeta8.1-8.2","TCRVbeta5.1-5.2", "TCRValpha2", "TCRValpha8.3-cloneB21.14", "TCRValpha8.3-cloneKT50", "TCRValpha11.1-11.2", "TCRgammadelta"), cols = c("goldenrod","green","blue", "gold", "purple", "orange"))

# Create a Violin plot and UMAP visualization of clusters with a specific resolution (1.0).

nb.cols <- 41
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.1")
clbycl_harmony_1 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors, label.box = T, label.color = "white") +
  labs(title = "Clustering with resolution 1")
gene_counts <- VlnPlot(TERVA2_harmony, features = "nFeature_RNA", group.by = "RNA_snn_res.1", cols = mycolors)
gene_counts + clbycl_harmony_1

# Remove the CD79a and Cd3d positive cluster. After removing this, the cell numbers per sample drop, especially in the aorta samples.
TERVA2_harmony <- subset(TERVA2_harmony, subset = RNA_snn_res.1 %in% c(0:11,13:40)) #Remove cluster 9 that's CD79a and Cd3d positive and might bias subsequent workflows
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.1")
clbycl_harmony_1 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors, label.box = T, label.color = "white") +
  labs(title = "Clustering with resolution 1 after removing Cd79a/Cd3d double-positive cluster")
clbycl_harmony_1

# Look at cell markers
Fcellannot <- FeaturePlot(TERVA2_harmony, features= c("Myh11", "Nkg7", "Rora", "Dcn", "Vcam1", "Rgs5", "Cd34", "Upk3b", "Pecam1", "Cd79a", "Cd3d", "Cd68", "Spp1", "S100a8", "S100a9", "Top2a")) 
Fadip <- FeaturePlot(TERVA2_harmony, features= c("Adipoq", "Pdgfra", "Ces1f", "Btc", "Apoe", "Cacna1a", "Prune2", "Mt2", "Tcf21"))
Fangio <- FeaturePlot(TERVA2_harmony, features= c("Cdh5", "Sdpr", "Egfl7", "Ptprb", "Ecscr", "Cldn5", "Icam2", "Slc9a3r2", "Myh11", "Cnn1", "Vcam1", "Cd36"))
CD4_CD8 <- FeaturePlot(TERVA2_harmony, order = T, features= c("Cd8a", "Cd8b1", "Cd4", "Ccl5", "Foxp3","Vps37b", "Ramp3", "Tcf7", "Rag1"))
ILC <-FeaturePlot(TERVA2_harmony, order = T, features= c("Icos", "Kit", "Il23r", "Rora", "Gata3"))
Plasmacells <- FeaturePlot(TERVA2_harmony, features = c("Ighd", "IgD", "CD138-Syndecan-1", "CD45R-B220", "Tnfrsf17", "Cxcr4"))
Monocytes <- FeaturePlot(TERVA2_harmony, order=T, features = c("Ly6c1", "Ly6c2", "Ly-6C", "Cx3cr1", "Ccr2", "Sell", "CD62L", "Cd209a", "I-A-I-E", "Spn", "CD43", "Treml4"), ncol = 4)
Macrophages <- FeaturePlot(TERVA2_harmony, order =T, features = c("Cd163","Cd80", "Cd86", "Ccr5", "CD11b", "CD11c", "Cd14", "CD15-SSEA-1", "Cd68", "Cd36", "Ptgs2", "Irf5", "Stat1", "Nos2", "Cxcr1", "Cxcr2", "Mrc1", "Irf4", "Stat6", "Socs3", "Sphk1", "Tlr8"))
Fcrs <- FeaturePlot(TERVA2_harmony, features = c("B2m","Ero1l", "Fcamr", "Fcer1a", "Fcer1g", "Fcer2a", "Fcgrt", "Fcrla", "Fcrlb", "Fcrls", "Fcrl1", "Fcrl5", "Fcrl6", "Fcgbp", "Fcgrt", "Fcgr2b", "Fcgr4", "Fcgr3", "Fcgr1"))
Igs_lc <- FeaturePlot(TERVA2_harmony, split.by = "Sample", features = c("Igkv3-4", "Igkv4-55", "Igkv1-135",  "Jchain", "Ighm"))
Pi16 <-FeaturePlot(TERVA2_harmony, features = c("Pi16"))
Cd74 <- FeaturePlot(TERVA2_harmony, features = c("Cd74", "CD44", "Col1a1", "Lum"))
Calca <-FeaturePlot(TERVA2_harmony, order =T, features = c("Calca", "Igf1"))
Meso <-FeaturePlot(TERVA2_harmony, order =T, features = c("Krt19", "Upk3b", "Msln"))

#T cells
DotPlot(TERVA2_harmony, scale = F, idents = c("3", "5", "20"), features= c("Cd8a", "Cd8b1", "Cd4", "Cd3d", "Cd3e", "Cd3g", "Tcf7", "Rag1", "Ccr9"), cols = c("blue", "gold"))
VlnPlot(TERVA2_harmony, slot = "counts", idents = c("3", "5", "20"), features= c("Cd8a", "Cd8b1", "Cd4", "Cd3d", "Cd3e", "Cd3g", "Tcf7", "Rag1", "Ccr9"), cols = c("blue", "gold", "purple"))
VlnPlot(TERVA2_harmony, assay = "ADT", slot = "counts", idents = c("3", "5", "20"), features= c("CD8a", "CD4"), cols = c("blue", "gold", "purple"))
FeaturePlot(TERVA2_harmony, order = T, features= c("Cd8a", "CD8a", "Cd8b1", "Cd4", "CD4", "Cd3d", "Cd3e", "Cd3g"))

DPcells_Cd8a <- WhichCells(TERVA2_harmony, expression = Cd8a > 0 & Cd4 > 0)
DPcells_Cd8b1 <- WhichCells(TERVA2_harmony, expression = Cd8b1 > 0 & Cd4 > 0)
DPcells_Cd8b1_Cd8a <- WhichCells(TERVA2_harmony, expression = Cd8b1 > 0 & Cd8a > 0 & Cd4 > 0)

DimPlot(TERVA2_harmony, cells.highlight = DPcells_Cd8a, order = T)

DP_avgexpr <- as.data.frame(AverageExpression(TERVA2_harmony, assays = "RNA", features= c("Cd8a", "Cd8b1", "Cd4", "Cd3d", "Cd3e", "Cd3g", "Tcf7"), return.seurat = FALSE, slot = "data"))
DP_avgexpr <- DP_avgexpr %>% dplyr::select(RNA.3, RNA.5, RNA.20)

#Notch3 expr in VSMCs
VlnPlot(TERVA2_harmony, slot = "counts", idents = c("17", "36", "30", "39", "28"), features= c("Notch3"))
DP_avgexpr_Notch3 <- as.data.frame(AverageExpression(TERVA2_harmony, assays = "RNA", features= c("Notch3"), return.seurat = FALSE, slot = "data"))
DP_avgexpr_Notch3 <- DP_avgexpr_Notch3 %>% dplyr::select(RNA.17, RNA.36, RNA.30, RNA.39, RNA.28)

#Fibroblast / EC markers
VlnPlot(TERVA2_harmony, slot = "counts", idents = c("4","28","3"), features= c("Pecam1", "Cdh5", "Col1a1", "Lum"))
Cd74fibro <- WhichCells(TERVA2_harmony, expression = Cd74 > 0 & Col1a1 > 0 & Lum > 0)
Ptprcendo <- WhichCells(TERVA2_harmony, expression = Ptprc > 0 & Pecam1 > 0 & Kdr > 0)
DimPlot(TERVA2_harmony, cells.highlight = Cd74fibro, order = T)
DimPlot(TERVA2_harmony, cells.highlight = Ptprcendo, order = T)

#Dendritic cells
VlnPlot(TERVA2_harmony, slot = "counts", idents = c("7","9","21"), features= c("Cd209a", "H2-DMb1", "Itgam", "Itgae", "Notch2", "Btla", "Clec9a", "Clec10a", "Sirpa", "Ccr7", "Mrc1"))

#Monocytes / Macrophages
VlnPlot(TERVA2_harmony, slot = "data", idents = c("1","17",32), features= c("Cd68", "Socs3", "Tgm2", "Mrc1", "Retnla", "Cd14", "Fcgr3"))

molecule_counts <- VlnPlot(TERVA2_harmony, features = "nCount_RNA", group.by = "RNA_snn_res.1", cols = mycolors)
gene_counts <- VlnPlot(TERVA2_harmony, features = "nFeature_RNA", group.by = "RNA_snn_res.1", cols = mycolors)
TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.1")
nb.cols <- 41
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
clbycl_harmony_1 <- DimPlot(TERVA2_harmony,label = T, shuffle = T, cols = mycolors, pt.size = 2, label.box = T, label.color = "white") +
  labs(title = "Clustering with resolution 1")
molecule_counts + gene_counts | clbycl_harmony_1

clbycl_harmony_1_tissue <- DimPlot(TERVA2_harmony,label = T, shuffle = T, split.by = "Sample", cols = mycolors, pt.size = 2, label.box = T, label.color = "white") +
  labs(title = "Clustering with resolution 1")
clbycl_harmony_1_tissue

# Check CD45 across the data

T2_harmony_ptprc <- FeaturePlot(TERVA2_harmony, features= c("Ptprc"), order = T)
T2_harmony_ptprc

#### Add disease and tissue group id ####

coldata <- as.data.frame(TERVA2_harmony@meta.data$Sample)

coldata <- mutate(coldata,
                  Status = case_when(
                    startsWith(TERVA2_harmony@meta.data$Sample, "OB") ~ "Obese",
                    startsWith(TERVA2_harmony@meta.data$Sample, "NOB") ~ "Non-Obese"
                  ))

group_id <- coldata$Status
TERVA2_harmony@meta.data$group_id <- as.factor(group_id)

#Create tissue ID's

tissuedata <- as.data.frame(TERVA2_harmony@meta.data$Sample)

tissuedata <- mutate(tissuedata,
                     Tissue = case_when(
                       endsWith(TERVA2_harmony@meta.data$Sample, "AO") ~ "Aorta",
                       endsWith(TERVA2_harmony@meta.data$Sample, "PVAT") ~ "PVAT",
                       endsWith(TERVA2_harmony@meta.data$Sample, "eWAT") ~ "eWAT",
                       endsWith(TERVA2_harmony@meta.data$Sample, "Spleen") ~ "Spleen"
                     ))

tissue_id <- tissuedata$Tissue
TERVA2_harmony@meta.data$tissue_id <- as.factor(tissue_id)

#### Predicting cell types from reference ####
# Note that also manual annotation was used based on the marker expression above

library(celldex)
library(SingleR)

mouse.ref <- celldex::MouseRNAseqData()
saveRDS(mouse.ref, "mouseref.rds")
mouse.ref <- readRDS("/scratch/project_2005050/Rstats/mouseref.rds")
imm.ref <- celldex::ImmGenData()

sce <- as.SingleCellExperiment((DietSeurat(TERVA2_harmony)))
sce

mouse.main <- SingleR(test= sce, assay.type.test = 1, ref = mouse.ref, labels = mouse.ref$label.main)
mouse.fine <- SingleR(test= sce, assay.type.test = 1, ref = mouse.ref, labels = mouse.ref$label.fine)

table(mouse.main$pruned.labels)
table(mouse.fine$pruned.labels)

TERVA2_harmony@meta.data$mouse.main <- mouse.main$pruned.labels
TERVA2_harmony@meta.data$mouse.fine <- mouse.fine$pruned.labels

TERVA2_harmony <- SetIdent(TERVA2_harmony, value = "mouse.main")
m_main <- DimPlot(TERVA2_harmony, label = T , repel = T, label.size = 3) + NoLegend()

TERVA2_harmony <- SetIdent(TERVA2_harmony, value = "mouse.fine")
m_fine <- DimPlot(TERVA2_harmony, label = T , repel = T, label.size = 3) + NoLegend()

m_fine + m_main

#### Differential gene expression between clusters ####

TERVA2_harmony <- SetIdent(TERVA2_harmony, value = "RNA_snn_res.1")
TERVA2_harmony_allmarkers_wilcoxon <- FindAllMarkers(TERVA2_harmony, verbose = T, min.cells.group = 10)
#saveRDS(TERVA2_harmony_allmarkers_wilcoxon, "TERVA2_harmony_allmarkers_wilcoxon_reQC.rds")
Cluster_markers <- readRDS("TERVA2_harmony_allmarkers_wilcoxon_reQC.rds")

split <- split(Cluster_markers, Cluster_markers$cluster)
list2env(split, envir = globalenv())

top6 <- Cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 6, wt = avg_log2FC)

clusterheatmap <- DoHeatmap(subset(TERVA2_harmony, downsample=200), slot = "scale.data", features = top6$gene, size = 4, angle = 90, label = T, group.colors = mycolors) + NoLegend() + 
  theme(axis.text.y = element_text(size=12)) 

#### Rename clusters ####

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.1")
TERVA2_harmony <- RenameIdents(TERVA2_harmony, 
                               "0" = "B cells", "1" = "Folr2+ Lyve1+ M2 Macrophages", "2" = "Ccl11+ Fibroblasts", "3" = "Lef1+ Tcf7+ Cd4+ T cells", 
                               "4" = "Gpihbp1+ Fabp4+ Endothelial cells", "5" = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells", "6" = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages", 
                               "7" = "Cd248+ Pi16+ Fibroblasts", "8" = "Mfap4+ Fibroblasts",
                               "9" = "Conventional Dendritic cells DC1", "10" = "Mgp+ Aebp1+ Activated fibroblasts", "11" = "Classical and non-classical monocytes", 
                               "13" = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells", "14" = "Foxp3+ Regulatory T cells", 
                               "15" = "Gpihbp1+ Fabp4+ Endothelial cells", "16" = "Innate lymphoid cells", 
                               "17" = "Vascular smooth muscle cells", "18" = "s100a9+/a8+ Granulocytes", "19" = "Intermediate monocytes",
                               "20" = "Vps37b+ Ramp3+ Cd8+ T memory cells", "21" = "Trem2+ Lgals3+ Macrophages", "22" = "Mki67+ Top2a+ Proliferating cells", 
                               "23" = "Natural killer cells", "24" = "Conventional Dendritic cells DC1", "25" = "Fscn1+ Apol7c+ Dendritic cells",
                               "26" = "Pf4+ Retnla+ Macrophages", "27" = "Plasma cells", "28" = "Vascular smooth muscle cells", "29" = "Mesothelial cells",
                               "30" = "Notch3 low VSMCs", "31" = "Il1b+ Fibroblasts", "32" = "Rgs5+ Endothelial cells", "33" = "Gpihbp1+ Fabp4+ Endothelial cells",
                               "34" = "Cd248+ Pi16+ Fibroblasts", "35" = "Gpihbp1+ Fabp4+ Endothelial cells", "36" = "Notch3 high VSMCs", "37" = "Pecam1+ Cd5+ Col1a1+ Lum+ cells",
                               "38" = "MAST cells", "39" = "Notch3 low VSMCs", "40" = "s100a9+/a8+ Granulocytes")

TERVA2_harmony$celltype.group <- paste(Idents(TERVA2_harmony), TERVA2_harmony$group_id, sep = "_")
TERVA2_harmony$celltype <- Idents(TERVA2_harmony)
celltypes <- as.vector(unique(TERVA2_harmony$celltype))
celltypes

#[1] "Vps37b+ Ramp3+ Cd8+ T memory cells"           "Natural killer cells"                         "Lef1+ Tcf7+ Cd4+ T cells"                    
#[4] "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells"           "Innate lymphoid cells"                        "Gpihbp1+ Fabp4+ Endothelial cells"           
#[7] "B cells"                                      "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages" "Classical and non-classical monocytes"       
#[10] "Mki67+ Top2a+ Proliferating cells"            "Vascular smooth muscle cells"                 "Trem2+ Lgals3+ Macrophages"                  
#[13] "Conventional Dendritic cells DC1"             "Cd248+ Pi16+ Fibroblasts"                     "Mfap4+ Fibroblasts"                          
#[16] "Folr2+ Lyve1+ M2 Macrophages"                 "Pf4+ Retnla+ Macrophages"                     "Foxp3+ Regulatory T cells"                   
#[19] "Fscn1+ Apol7c+ Dendritic cells"               "Intermediate monocytes"                       "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells"            
#[22] "s100a9+/a8+ Granulocytes"                     "Ccl11+ Fibroblasts"                           "Plasma cells"                                
#[25] "Mgp+ Aebp1+ Activated fibroblasts"            "Il1b+ Fibroblasts"                            "Mesothelial cells"                           
#[28] "Notch3 high VSMCs"                            "MAST cells"                                   "Rgs5+ Endothelial cells"                     
#[31] "Notch3 low VSMCs"                             "Pecam1+ Cd5+ Col1a1+ Lum+ cells"             

Idents(TERVA2_harmony) <- "celltype"
nb.cols <- 35
mycolors <- colorRampPalette(colorBlindness::paletteMartin)(nb.cols)
reannot <- DimPlot(TERVA2_harmony, pt.size = 1.2, cols = mycolors, order = F)
reannot + 
  labs(title = "Cell annotations of the integrated data") +
  theme(plot.title = element_text(size=40, hjust= 0.5),
        axis.title = element_text(size=20),
        axis.text = element_text(size=20)
        ) +
  guides(col = guide_legend(ncol = 1))

# Count module score and visualize it
TERVA2_harmony <- AddModuleScore(TERVA2_harmony, features = list("Apoe",    "Lyz2",    "Cd74",    "C1qb",    "C1qc",    "C1qa",    "Ccl6",    "Ctsc",    "Gm46603", "H2-Aa",   "Fcer1g",  "F13a1",   "Tyrobp",  "H2-Ab1",  "Mrc1",    "H2-Eb1", 
                                                                 "Retnla",  "C4b",     "Wfdc17",  "Mgp", "Cd248", "Pi16"), name = "Pi16fibro_topgenes")

moduleplot <- FeaturePlot(TERVA2_harmony,
            features = "Pi16fibro_topgenes22", label = F, repel = TRUE, pt.size = 1.2) +
  theme(legend.text = element_text(size=15)) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + 
  labs(title = "Module expression of top DE genes of PVAT-derived activated fibroblasts")

reannot + labs(title = "Main cell populations and subclustering of the integrated data") + NoLegend() + moduleplot

#The same as above but without the legend
Idents(TERVA2_harmony) <- "celltype"
reannot <- DimPlot(TERVA2_harmony, label = T, repel = T, label.box = T, label.color = "white", pt.size = 1.2, cols = mycolors, order = T) + NoLegend()
reannot + 
  labs(title = "Clustering of the integrated data") +
  theme(plot.title = element_text(size=40, hjust= 0.5),
        axis.title = element_text(size=20),
        axis.text = element_text(size=20))

reannot2 <- DimPlot(TERVA2_harmony, split.by = "tissue_id", pt.size = 1.2, cols = mycolors, order = T, label.color = "white", label.size = 1.5) + NoLegend()

#Plot by group

diseasestates <- DimPlot(TERVA2_harmony, reduction = "umap", group.by = "group_id", pt.size = 1.2, label = T, label.box =  T, label.color = "white", label.size = 4, repel = T, cols = c("#009292", "#490092")) + NoLegend()
diseasestates

((reannot + diseasestates) / (reannot2)) | clusterheatmap

#### Tissue specific DE nalyses ####

Idents(TERVA2_harmony) <- "tissue_id"
AO <- subset(TERVA2_harmony, idents = "Aorta")
PVAT <- subset(TERVA2_harmony, idents = "PVAT")
eWAT <- subset(TERVA2_harmony, idents = "eWAT")
Spleen <- subset(TERVA2_harmony, idents = "Spleen")

#### Aorta ####

# Aosta obese vs non-obese DE's per cell type ####
Idents(AO) <- "celltype.group"
AO_Vps37bRamp3Cd8Tmem <- FindMarkers(AO, only.pos = F, ident.1 = "Vps37b+ Ramp3+ Cd8+ T memory cells_Obese", ident.2 = "Vps37b+ Ramp3+ Cd8+ T memory cells_Non-Obese", verbose = FALSE)
AO_NK <- FindMarkers(AO, only.pos = F, ident.1 = "Natural killer cells_Obese", ident.2 = "Natural killer cells_Non-Obese", verbose = FALSE)
AO_Lef1Tcf7Cd4Tcells <- FindMarkers(AO, only.pos = F, ident.1 = "Lef1+ Tcf7+ Cd4+ T cells_Obese", ident.2 = "Lef1+ Tcf7+ Cd4+ T cells_Non-Obese", verbose = FALSE)
AO_Lef1Tcf7Cd8Tcells <- FindMarkers(AO, only.pos = F, ident.1 = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells_Obese", ident.2 = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells_Non-Obese", verbose = FALSE)
AO_Cd8Ccl5Nkg7CytotoxicTcells <- FindMarkers(AO, only.pos = F, ident.1 = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells_Obese", ident.2 = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells_Non-Obese", verbose = FALSE)
AO_Foxp3RegTcells <- FindMarkers(AO, only.pos = F, ident.1 = "Foxp3+ Regulatory T cells_Obese", ident.2 = "Foxp3+ Regulatory T cells_Non-Obese", verbose = FALSE)
AO_ILC <- FindMarkers(AO, only.pos = F, ident.1 = "Innate lymphoid cells_Obese", ident.2 = "Innate lymphoid cells_Non-Obese", verbose = FALSE)
AO_Bcells <- FindMarkers(AO, only.pos = F, ident.1 = "B cells_Obese", ident.2 = "B cells_Non-Obese", verbose = FALSE)
AO_Plasmacells <- FindMarkers(AO, only.pos = F, ident.1 = "Plasma cells_Obese", ident.2 = "Plasma cells_Non-Obese", verbose = FALSE)
AO_Macro_Folr2Lyve1 <- FindMarkers(AO, only.pos = F, ident.1 = "Folr2+ Lyve1+ M2 Macrophages_Obese", ident.2 = "Folr2+ Lyve1+ M2 Macrophages_Non-Obese", verbose = FALSE)
AO_Macro_Trem2Lgals3 <- FindMarkers(AO, only.pos = F, ident.1 = "Trem2+ Lgals3+ Macrophages_Obese", ident.2 = "Trem2+ Lgals3+ Macrophages_Non-Obese", verbose = FALSE)
AO_Macro_Pf4Retnla <- FindMarkers(AO, only.pos = F, ident.1 = "Pf4+ Retnla+ Macrophages_Obese", ident.2 = "Pf4+ Retnla+ Macrophages_Non-Obese", verbose = FALSE)
AO_Macro_infl <- FindMarkers(AO, only.pos = F, ident.1 = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages_Obese", ident.2 = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages_Non-Obese", verbose = FALSE)
AO_Mono_inter <- FindMarkers(AO, only.pos = F, ident.1 = "Intermediate monocytes_Obese", ident.2 = "Intermediate monocytes_Non-Obese", verbose = FALSE)
AO_Mono_clandnc <- FindMarkers(AO, only.pos = F, ident.1 = "Classical and non-classical monocytes_Obese", ident.2 = "Classical and non-classical monocytes_Non-Obese", verbose = FALSE)
AO_Conv_DC1 <- FindMarkers(AO, only.pos = F, ident.1 = "Conventional Dendritic cells DC1_Obese", ident.2 = "Conventional Dendritic cells DC1_Non-Obese", verbose = FALSE)
AO_Migr_DC <- FindMarkers(AO, only.pos = F, ident.1 = "Fscn1+ Apol7c+ Dendritic cells_Obese", ident.2 = "Fscn1+ Apol7c+ Dendritic cells_Non-Obese", verbose = FALSE)
AO_Mki67Top2acells <- FindMarkers(AO, only.pos = F, ident.1 = "Mki67+ Top2a+ Proliferating cells_Obese", ident.2 = "Mki67+ Top2a+ Proliferating cells_Non-Obese", verbose = FALSE)
AO_Fibro_Pi16 <- FindMarkers(AO, only.pos = F, ident.1 = "Cd248+ Pi16+ Fibroblasts_Obese", ident.2 = "Cd248+ Pi16+ Fibroblasts_Non-Obese", verbose = FALSE)
AO_Fibro_Ccl11 <- FindMarkers(AO, only.pos = F, ident.1 = "Ccl11+ Fibroblasts_Obese", ident.2 = "Ccl11+ Fibroblasts_Non-Obese", verbose = FALSE)
AO_Fibro_Mgp <- FindMarkers(AO, only.pos = F, ident.1 = "Mgp+ Aebp1+ Activated fibroblasts_Obese", ident.2 = "Mgp+ Aebp1+ Activated fibroblasts_Non-Obese", verbose = FALSE)
AO_Fibro_Mfap <- FindMarkers(AO, only.pos = F, ident.1 = "Mfap4+ Fibroblasts_Obese", ident.2 = "Mfap4+ Fibroblasts_Non-Obese", verbose = FALSE)
AO_s100a9_a8Granulocytes <- FindMarkers(AO, only.pos = F, ident.1 = "s100a9+/a8+ Granulocytes_Obese", ident.2 = "s100a9+/a8+ Granulocytes_Non-Obese", verbose = FALSE)
AO_Gpihbp1Fabp4EC <- FindMarkers(AO, only.pos = F, ident.1 = "Gpihbp1+ Fabp4+ Endothelial cells_Obese", ident.2 = "Gpihbp1+ Fabp4+ Endothelial cells_Non-Obese", verbose = FALSE)
AO_VSMCs <- FindMarkers(AO, only.pos = F, ident.1 = "Vascular smooth muscle cells_Obese", ident.2 = "Vascular smooth muscle cells_Non-Obese", verbose = FALSE)


sctlist <- c("AO_Vps37bRamp3Cd8Tmem", "AO_NK", "AO_Lef1Tcf7Cd4Tcells", "AO_Lef1Tcf7Cd8Tcells", "AO_Cd8Ccl5Nkg7CytotoxicTcells", "AO_Foxp3RegTcells", "AO_ILC", "AO_Bcells", "AO_Plasmacells",
             "AO_Macro_Folr2Lyve1", "AO_Macro_Trem2Lgals3", "AO_Macro_Pf4Retnla", "AO_Macro_infl", "AO_Mono_inter", "AO_Mono_clandnc", "AO_Conv_DC1", "AO_Migr_DC", "AO_Mki67Top2acells", "AO_Fibro_Pi16",
             "AO_Fibro_Ccl11", "AO_Fibro_Mgp", "AO_Fibro_Mfap", "AO_s100a9_a8Granulocytes", "AO_Gpihbp1Fabp4EC", "AO_VSMCs")

for(i in 1:length(sctlist)) {                
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}

DimPlot(AO, reduction = "umap", label = T , repel = T, label.size = 3) + NoLegend()

top10 <- AO_Cd4posT_sub %>% filter(AO_Cd4posT_sub$p_val_adj < 5e-2) %>% top_n(n = 10, wt = avg_log2FC) %>% rownames()
DoHeatmap(AO, features = top10)

cells <- c("Hdc+ Cpa3+ MAST cells", "Mgp+ Fibroblasts EC's with shear stress markers", "Plasma cells", "Conv DC2", "NK's", "Conv DC1", "Ras5+ EC's", 
"Dividing cells", "Granulocytes", "Rpl-Rps+ cells", "Lgals3+ Macrophages", "ILC's", "Cd8+ Ccl5+ Teffs", "Fibroblasts activated", "VSMC's", 
"Cd4+ Foxp3+ Tregs", "Pi16+ Fibroblasts", "Classical and non-classical Monocytes", "Intermediate monocytes", "Macrophages activated",
"EC's", "Cd4+ Tcells", "Cd8+ Tcells", "Macrophages", "Fibroblasts", "B cells")

AO1 <-  EnhancedVolcano(AO_Vps37bRamp3Cd8Tmem,
                          lab = rownames(AO_Vps37bRamp3Cd8Tmem),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "Cd8+ Tmem cells",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)

AO2 <-  EnhancedVolcano(AO_Migr_DC,
                          lab = rownames(AO_Migr_DC),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "Fscn1+ Apol7c+ Dendritic cells",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)

AO1 + AO2


#### eWAT ####

Idents(eWAT) <- "celltype.group"
DimPlot(eWAT, reduction = "umap", group.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()
p2 <- DimPlot(eWAT, reduction = "umap", split.by = "group_id", group.by = "seurat_clusters", pt.size = 1.2, label = T , repel = T, label.size = 4) + NoLegend() + ggtitle("eWAT")
l2 <- FeaturePlot(eWAT, features = c("Ly6a", "Ly-6A-E-Sca-1", "Lgals3")) +
  ggtitle("Ly-6A-E-Sca-1", subtitle = "eWAT")


Idents(eWAT) <- "celltype"
macros_eWAT <- subset(eWAT, idents = c("Macrophages", "Macrophages activated", "Lgals3+ Macrophages"))
Idents(macros_eWAT) <- "group_id"

#Look at some macrophage markers to differentiate the populations
FeaturePlot(macros_eWAT, split.by = "group_id", features = c("Trem2", "Lgals3"), blend = TRUE, cols = c("navy", "darkgoldenrod1"))
FeaturePlot(macros_eWAT, split.by = "group_id", features = c("Lamp2", "Mrc1"), blend = TRUE, cols = c("navy", "darkgoldenrod1"))
RidgePlot(macros_eWAT, group.by = "group_id", features = c("Trem2", "Lgals3"))


Idents(eWAT) <- "celltype.group"
eWAT_Vps37bRamp3Cd8Tmem <- FindMarkers(eWAT, only.pos = F, ident.1 = "Vps37b+ Ramp3+ Cd8+ T memory cells_Obese", ident.2 = "Vps37b+ Ramp3+ Cd8+ T memory cells_Non-Obese", verbose = FALSE)
eWAT_NK <- FindMarkers(eWAT, only.pos = F, ident.1 = "Natural killer cells_Obese", ident.2 = "Natural killer cells_Non-Obese", verbose = FALSE)
eWAT_Lef1Tcf7Cd4Tcells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Lef1+ Tcf7+ Cd4+ T cells_Obese", ident.2 = "Lef1+ Tcf7+ Cd4+ T cells_Non-Obese", verbose = FALSE)
eWAT_Lef1Tcf7Cd8Tcells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells_Obese", ident.2 = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells_Non-Obese", verbose = FALSE)
eWAT_Cd8Ccl5Nkg7CytotoxicTcells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells_Obese", ident.2 = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells_Non-Obese", verbose = FALSE)
eWAT_Foxp3RegTcells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Foxp3+ Regulatory T cells_Obese", ident.2 = "Foxp3+ Regulatory T cells_Non-Obese", verbose = FALSE)
eWAT_ILC <- FindMarkers(eWAT, only.pos = F, ident.1 = "Innate lymphoid cells_Obese", ident.2 = "Innate lymphoid cells_Non-Obese", verbose = FALSE)
eWAT_Bcells <- FindMarkers(eWAT, only.pos = F, ident.1 = "B cells_Obese", ident.2 = "B cells_Non-Obese", verbose = FALSE)
eWAT_Plasmacells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Plasma cells_Obese", ident.2 = "Plasma cells_Non-Obese", verbose = FALSE)
eWAT_Macro_Folr2Lyve1 <- FindMarkers(eWAT, only.pos = F, ident.1 = "Folr2+ Lyve1+ M2 Macrophages_Obese", ident.2 = "Folr2+ Lyve1+ M2 Macrophages_Non-Obese", verbose = FALSE)
eWAT_Macro_Trem2Lgals3 <- FindMarkers(eWAT, only.pos = F, ident.1 = "Trem2+ Lgals3+ Macrophages_Obese", ident.2 = "Trem2+ Lgals3+ Macrophages_Non-Obese", verbose = FALSE)
eWAT_Macro_Pf4Retnla <- FindMarkers(eWAT, only.pos = F, ident.1 = "Pf4+ Retnla+ Macrophages_Obese", ident.2 = "Pf4+ Retnla+ Macrophages_Non-Obese", verbose = FALSE)
eWAT_Macro_infl <- FindMarkers(eWAT, only.pos = F, ident.1 = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages_Obese", ident.2 = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages_Non-Obese", verbose = FALSE)
eWAT_Mono_inter <- FindMarkers(eWAT, only.pos = F, ident.1 = "Intermediate monocytes_Obese", ident.2 = "Intermediate monocytes_Non-Obese", verbose = FALSE)
eWAT_Mono_clandnc <- FindMarkers(eWAT, only.pos = F, ident.1 = "Classical and non-classical monocytes_Obese", ident.2 = "Classical and non-classical monocytes_Non-Obese", verbose = FALSE)
eWAT_Conv_DC1 <- FindMarkers(eWAT, only.pos = F, ident.1 = "Conventional Dendritic cells DC1_Obese", ident.2 = "Conventional Dendritic cells DC1_Non-Obese", verbose = FALSE)
eWAT_Migr_DC <- FindMarkers(eWAT, only.pos = F, ident.1 = "Fscn1+ Apol7c+ Dendritic cells_Obese", ident.2 = "Fscn1+ Apol7c+ Dendritic cells_Non-Obese", verbose = FALSE)
eWAT_Mki67Top2acells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Mki67+ Top2a+ Proliferating cells_Obese", ident.2 = "Mki67+ Top2a+ Proliferating cells_Non-Obese", verbose = FALSE)
eWAT_Fibro_Pi16 <- FindMarkers(eWAT, only.pos = F, ident.1 = "Cd248+ Pi16+ Fibroblasts_Obese", ident.2 = "Cd248+ Pi16+ Fibroblasts_Non-Obese", verbose = FALSE)
eWAT_Fibro_Ccl11 <- FindMarkers(eWAT, only.pos = F, ident.1 = "Ccl11+ Fibroblasts_Obese", ident.2 = "Ccl11+ Fibroblasts_Non-Obese", verbose = FALSE)
eWAT_Fibro_Mgp <- FindMarkers(eWAT, only.pos = F, ident.1 = "Mgp+ Aebp1+ Activated fibroblasts_Obese", ident.2 = "Mgp+ Aebp1+ Activated fibroblasts_Non-Obese", verbose = FALSE)
eWAT_Fibro_Il1b <- FindMarkers(eWAT, only.pos = F, ident.1 = "Il1b+ Fibroblasts_Obese", ident.2 = "Il1b+ Fibroblasts_Non-Obese", verbose = FALSE)
eWAT_Fibro_Mfap <- FindMarkers(eWAT, only.pos = F, ident.1 = "Mfap4+ Fibroblasts_Obese", ident.2 = "Mfap4+ Fibroblasts_Non-Obese", verbose = FALSE)
eWAT_s100a9_a8Granulocytes <- FindMarkers(eWAT, only.pos = F, ident.1 = "s100a9+/a8+ Granulocytes_Obese", ident.2 = "s100a9+/a8+ Granulocytes_Non-Obese", verbose = FALSE)
eWAT_Gpihbp1Fabp4EC <- FindMarkers(eWAT, only.pos = F, ident.1 = "Gpihbp1+ Fabp4+ Endothelial cells_Obese", ident.2 = "Gpihbp1+ Fabp4+ Endothelial cells_Non-Obese", verbose = FALSE)
eWAT_Rgs5EC <- FindMarkers(eWAT, only.pos = F, ident.1 = "Rgs5+ Endothelial cells_Obese", ident.2 = "Rgs5+ Endothelial cells_Non-Obese", verbose = FALSE)
eWAT_VSMCs <- FindMarkers(eWAT, only.pos = F, ident.1 = "Vascular smooth muscle cells_Obese", ident.2 = "Vascular smooth muscle cells_Non-Obese", verbose = FALSE)
eWAT_Mesothelialcells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Mesothelial cells_Obese", ident.2 = "Mesothelial cells_Non-Obese", verbose = FALSE)
eWAT_Pecam1Cd5Col1a1Lumcells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Pecam1+ Cd5+ Col1a1+ Lum+ cells_Obese", ident.2 = "Pecam1+ Cd5+ Col1a1+ Lum+ cells_Non-Obese", verbose = FALSE)

sctlist <- c("eWAT_Vps37bRamp3Cd8Tmem", "eWAT_NK", "eWAT_Lef1Tcf7Cd4Tcells", "eWAT_Lef1Tcf7Cd8Tcells", "eWAT_Cd8Ccl5Nkg7CytotoxicTcells", "eWAT_Foxp3RegTcells", "eWAT_ILC", "eWAT_Bcells", "eWAT_Plasmacells",
             "eWAT_Macro_Folr2Lyve1", "eWAT_Macro_Trem2Lgals3", "eWAT_Macro_Pf4Retnla", "eWAT_Macro_infl", "eWAT_Mono_inter", "eWAT_Mono_clandnc", "eWAT_Conv_DC1", "eWAT_Migr_DC", "eWAT_Mki67Top2acells",
             "eWAT_Fibro_Pi16", "eWAT_Fibro_Ccl11", "eWAT_Fibro_Mgp", "eWAT_Fibro_Il1b", "eWAT_Fibro_Mfap", "eWAT_s100a9_a8Granulocytes", "eWAT_Gpihbp1Fabp4EC", "eWAT_Rgs5EC", "eWAT_VSMCs", "eWAT_Mesothelialcells", 
             "eWAT_Pecam1Cd5Col1a1Lumcells")

for(i in 1:length(sctlist)) {                
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}

eWAT1 <-  EnhancedVolcano(eWAT_Fibro_Mgp,
                          lab = rownames(eWAT_Fibro_Mgp),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "activated Mgp+ fibroblasts",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)

eWAT2 <-  EnhancedVolcano(eWAT_Macro_Trem2Lgals3,
                          lab = rownames(eWAT_Macro_Trem2Lgals3),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "Trem2+ macrophages",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)

eWAT3 <-  EnhancedVolcano(eWAT_VSMCs,
                          lab = rownames(eWAT_VSMCs),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "VSMCs",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)


eWAT1 + eWAT2 + eWAT3

#### PVAT ####
Idents(PVAT) <- "celltype.group"
DimPlot(PVAT, reduction = "umap", group.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()

Idents(PVAT) <- "celltype.group"
PVAT_Vps37bRamp3Cd8Tmem <- FindMarkers(PVAT, only.pos = F, ident.1 = "Vps37b+ Ramp3+ Cd8+ T memory cells_Obese", ident.2 = "Vps37b+ Ramp3+ Cd8+ T memory cells_Non-Obese", verbose = FALSE)
PVAT_NK <- FindMarkers(PVAT, only.pos = F, ident.1 = "Natural killer cells_Obese", ident.2 = "Natural killer cells_Non-Obese", verbose = FALSE)
PVAT_Lef1Tcf7Cd4Tcells <- FindMarkers(PVAT, only.pos = F, ident.1 = "Lef1+ Tcf7+ Cd4+ T cells_Obese", ident.2 = "Lef1+ Tcf7+ Cd4+ T cells_Non-Obese", verbose = FALSE)
PVAT_Lef1Tcf7Cd8Tcells <- FindMarkers(PVAT, only.pos = F, ident.1 = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells_Obese", ident.2 = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells_Non-Obese", verbose = FALSE)
PVAT_Cd8Ccl5Nkg7CytotoxicTcells <- FindMarkers(PVAT, only.pos = F, ident.1 = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells_Obese", ident.2 = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells_Non-Obese", verbose = FALSE)
PVAT_Foxp3RegTcells <- FindMarkers(PVAT, only.pos = F, ident.1 = "Foxp3+ Regulatory T cells_Obese", ident.2 = "Foxp3+ Regulatory T cells_Non-Obese", verbose = FALSE)
PVAT_ILC <- FindMarkers(PVAT, only.pos = F, ident.1 = "Innate lymphoid cells_Obese", ident.2 = "Innate lymphoid cells_Non-Obese", verbose = FALSE)
PVAT_Bcells <- FindMarkers(PVAT, only.pos = F, ident.1 = "B cells_Obese", ident.2 = "B cells_Non-Obese", verbose = FALSE)
PVAT_Plasmacells <- FindMarkers(PVAT, only.pos = F, ident.1 = "Plasma cells_Obese", ident.2 = "Plasma cells_Non-Obese", verbose = FALSE)
PVAT_Macro_Folr2Lyve1 <- FindMarkers(PVAT, only.pos = F, ident.1 = "Folr2+ Lyve1+ M2 Macrophages_Obese", ident.2 = "Folr2+ Lyve1+ M2 Macrophages_Non-Obese", verbose = FALSE)
PVAT_Macro_Trem2Lgals3 <- FindMarkers(PVAT, only.pos = F, ident.1 = "Trem2+ Lgals3+ Macrophages_Obese", ident.2 = "Trem2+ Lgals3+ Macrophages_Non-Obese", verbose = FALSE)
PVAT_Macro_Pf4Retnla <- FindMarkers(PVAT, only.pos = F, ident.1 = "Pf4+ Retnla+ Macrophages_Obese", ident.2 = "Pf4+ Retnla+ Macrophages_Non-Obese", verbose = FALSE)
PVAT_Macro_infl <- FindMarkers(PVAT, only.pos = F, ident.1 = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages_Obese", ident.2 = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages_Non-Obese", verbose = FALSE)
PVAT_Mono_inter <- FindMarkers(PVAT, only.pos = F, ident.1 = "Intermediate monocytes_Obese", ident.2 = "Intermediate monocytes_Non-Obese", verbose = FALSE)
PVAT_Mono_clandnc <- FindMarkers(PVAT, only.pos = F, ident.1 = "Classical and non-classical monocytes_Obese", ident.2 = "Classical and non-classical monocytes_Non-Obese", verbose = FALSE)
PVAT_Conv_DC1 <- FindMarkers(PVAT, only.pos = F, ident.1 = "Conventional Dendritic cells DC1_Obese", ident.2 = "Conventional Dendritic cells DC1_Non-Obese", verbose = FALSE)
PVAT_Migr_DC <- FindMarkers(PVAT, only.pos = F, ident.1 = "Fscn1+ Apol7c+ Dendritic cells_Obese", ident.2 = "Fscn1+ Apol7c+ Dendritic cells_Non-Obese", verbose = FALSE)
PVAT_Mki67Top2acells <- FindMarkers(PVAT, only.pos = F, ident.1 = "Mki67+ Top2a+ Proliferating cells_Obese", ident.2 = "Mki67+ Top2a+ Proliferating cells_Non-Obese", verbose = FALSE)
PVAT_Fibro_Pi16 <- FindMarkers(PVAT, only.pos = F, ident.1 = "Cd248+ Pi16+ Fibroblasts_Obese", ident.2 = "Cd248+ Pi16+ Fibroblasts_Non-Obese", verbose = FALSE)
PVAT_Fibro_Ccl11 <- FindMarkers(PVAT, only.pos = F, ident.1 = "Ccl11+ Fibroblasts_Obese", ident.2 = "Ccl11+ Fibroblasts_Non-Obese", verbose = FALSE)
PVAT_Fibro_Mgp <- FindMarkers(PVAT, only.pos = F, ident.1 = "Mgp+ Aebp1+ Activated fibroblasts_Obese", ident.2 = "Mgp+ Aebp1+ Activated fibroblasts_Non-Obese", verbose = FALSE)
PVAT_Fibro_Mfap <- FindMarkers(PVAT, only.pos = F, ident.1 = "Mfap4+ Fibroblasts_Obese", ident.2 = "Mfap4+ Fibroblasts_Non-Obese", verbose = FALSE)
PVAT_s100a9_a8Granulocytes <- FindMarkers(PVAT, only.pos = F, ident.1 = "s100a9+/a8+ Granulocytes_Obese", ident.2 = "s100a9+/a8+ Granulocytes_Non-Obese", verbose = FALSE)
PVAT_Gpihbp1Fabp4EC <- FindMarkers(PVAT, only.pos = F, ident.1 = "Gpihbp1+ Fabp4+ Endothelial cells_Obese", ident.2 = "Gpihbp1+ Fabp4+ Endothelial cells_Non-Obese", verbose = FALSE)
PVAT_VSMCs <- FindMarkers(PVAT, only.pos = F, ident.1 = "Vascular smooth muscle cells_Obese", ident.2 = "Vascular smooth muscle cells_Non-Obese", verbose = FALSE)
PVAT_Notch3lowVSMCs <- FindMarkers(PVAT, only.pos = F, ident.1 = "Notch3 low VSMCs_Obese", ident.2 = "Notch3 low VSMCs_Non-Obese", verbose = FALSE)
PVAT_Notch3highVSMCs <- FindMarkers(PVAT, only.pos = F, ident.1 = "Notch3 high VSMCs_Obese", ident.2 = "Notch3 high VSMCs_Non-Obese", verbose = FALSE)
PVAT_Mesothelialcells <- FindMarkers(PVAT, only.pos = F, ident.1 = "Mesothelial cells_Obese", ident.2 = "Mesothelial cells_Non-Obese", verbose = FALSE)
PVAT_MASTcells <- FindMarkers(PVAT, only.pos = F, ident.1 = "MAST cells_Obese", ident.2 = "MAST cells_Non-Obese", verbose = FALSE)


sctlist <- c("PVAT_Vps37bRamp3Cd8Tmem", "PVAT_NK", "PVAT_Lef1Tcf7Cd4Tcells", "PVAT_Lef1Tcf7Cd8Tcells", "PVAT_Cd8Ccl5Nkg7CytotoxicTcells", "PVAT_Foxp3RegTcells", "PVAT_ILC", "PVAT_Bcells", "PVAT_Plasmacells", 
             "PVAT_Macro_Folr2Lyve1", "PVAT_Macro_Trem2Lgals3", "PVAT_Macro_Pf4Retnla", "PVAT_Macro_infl", "PVAT_Mono_inter", "PVAT_Mono_clandnc", "PVAT_Conv_DC1", "PVAT_Migr_DC", "PVAT_Mki67Top2acells", "PVAT_Fibro_Pi16",
             "PVAT_Fibro_Ccl11", "PVAT_Fibro_Mgp", "PVAT_Fibro_Mfap", "PVAT_s100a9_a8Granulocytes", "PVAT_Gpihbp1Fabp4EC", "PVAT_VSMCs", "PVAT_Notch3lowVSMCs", "PVAT_Notch3highVSMCs", "PVAT_Mesothelialcells", "PVAT_MASTcells")

for(i in 1:length(sctlist)) {                
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}

PVAT1 <-  EnhancedVolcano(PVAT_Fibro_Mgp,
                          lab = rownames(PVAT_Fibro_Mgp),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "activated Mgp+ fibroblasts",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)

PVAT2 <-  EnhancedVolcano(PVAT_Fibro_Pi16,
                          lab = rownames(PVAT_Fibro_Pi16),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "Pi16+ progenitor fibroblasts",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)

PVAT3 <-  EnhancedVolcano(PVAT_Gpihbp1Fabp4EC,
                          lab = rownames(PVAT_Gpihbp1Fabp4EC),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "ECs",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)

PVAT4 <-  EnhancedVolcano(PVAT_Mono_inter,
                          lab = rownames(PVAT_Mono_inter),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Differentially expressed genes in',
                          subtitle = "intermediate monocytes",
                          pCutoff = 1e-2,
                          FCcutoff = 0.4,
                          pointSize = 3.0,
                          col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                          labSize = 6.0,
                          titleLabSize = 15.0)


(PVAT1 + PVAT2) / (PVAT3 + PVAT4)

#### Spleen ####

Idents(Spleen) <- "celltype.group"
DimPlot(Spleen, reduction = "umap", group.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()
p4 <-DimPlot(Spleen, reduction = "umap", split.by = "group_id", group.by = "seurat_clusters", pt.size = 1.2, label = T , repel = T, label.size = 4) + NoLegend() + ggtitle("Spleen")
l4 <- FeaturePlot(Spleen, features = c("Ly6a", "Ly-6A-E-Sca-1")) + 
  ggtitle("Ly-6A-E-Sca-1", subtitle = "Spleen")


Idents(Spleen) <- "celltype.group"
SP_NK <- FindMarkers(Spleen, only.pos = F, ident.1 = "Natural killer cells_Obese", ident.2 = "Natural killer cells_Non-Obese", verbose = FALSE)
SP_Lef1Tcf7Cd4Tcells <- FindMarkers(Spleen, only.pos = F, ident.1 = "Lef1+ Tcf7+ Cd4+ T cells_Obese", ident.2 = "Lef1+ Tcf7+ Cd4+ T cells_Non-Obese", verbose = FALSE)
SP_Lef1Tcf7Cd8Tcells <- FindMarkers(Spleen, only.pos = F, ident.1 = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells_Obese", ident.2 = "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells_Non-Obese", verbose = FALSE)
SP_Cd8Ccl5Nkg7CytotoxicTcells <- FindMarkers(Spleen, only.pos = F, ident.1 = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells_Obese", ident.2 = "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells_Non-Obese", verbose = FALSE)
SP_Foxp3RegTcells <- FindMarkers(Spleen, only.pos = F, ident.1 = "Foxp3+ Regulatory T cells_Obese", ident.2 = "Foxp3+ Regulatory T cells_Non-Obese", verbose = FALSE)
SP_ILC <- FindMarkers(Spleen, only.pos = F, ident.1 = "Innate lymphoid cells_Obese", ident.2 = "Innate lymphoid cells_Non-Obese", verbose = FALSE)
SP_Bcells <- FindMarkers(Spleen, only.pos = F, ident.1 = "B cells_Obese", ident.2 = "B cells_Non-Obese", verbose = FALSE)
SP_Plasmacells <- FindMarkers(Spleen, only.pos = F, ident.1 = "Plasma cells_Obese", ident.2 = "Plasma cells_Non-Obese", verbose = FALSE)
SP_Macro_Folr2Lyve1 <- FindMarkers(Spleen, only.pos = F, ident.1 = "Folr2+ Lyve1+ M2 Macrophages_Obese", ident.2 = "Folr2+ Lyve1+ M2 Macrophages_Non-Obese", verbose = FALSE)
SP_Macro_Pf4Retnla <- FindMarkers(Spleen, only.pos = F, ident.1 = "Pf4+ Retnla+ Macrophages_Obese", ident.2 = "Pf4+ Retnla+ Macrophages_Non-Obese", verbose = FALSE)
SP_Macro_infl <- FindMarkers(Spleen, only.pos = F, ident.1 = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages_Obese", ident.2 = "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages_Non-Obese", verbose = FALSE)
SP_Mono_inter <- FindMarkers(Spleen, only.pos = F, ident.1 = "Intermediate monocytes_Obese", ident.2 = "Intermediate monocytes_Non-Obese", verbose = FALSE)
SP_Mono_clandnc <- FindMarkers(Spleen, only.pos = F, ident.1 = "Classical and non-classical monocytes_Obese", ident.2 = "Classical and non-classical monocytes_Non-Obese", verbose = FALSE)
SP_Conv_DC1 <- FindMarkers(Spleen, only.pos = F, ident.1 = "Conventional Dendritic cells DC1_Obese", ident.2 = "Conventional Dendritic cells DC1_Non-Obese", verbose = FALSE)
SP_Migr_DC <- FindMarkers(Spleen, only.pos = F, ident.1 = "Fscn1+ Apol7c+ Dendritic cells_Obese", ident.2 = "Fscn1+ Apol7c+ Dendritic cells_Non-Obese", verbose = FALSE)
SP_Mki67Top2acells <- FindMarkers(Spleen, only.pos = F, ident.1 = "Mki67+ Top2a+ Proliferating cells_Obese", ident.2 = "Mki67+ Top2a+ Proliferating cells_Non-Obese", verbose = FALSE)
SP_s100a9_a8Granulocytes <- FindMarkers(Spleen, only.pos = F, ident.1 = "s100a9+/a8+ Granulocytes_Obese", ident.2 = "s100a9+/a8+ Granulocytes_Non-Obese", verbose = FALSE)
SP_Gpihbp1Fabp4EC <- FindMarkers(Spleen, only.pos = F, ident.1 = "Gpihbp1+ Fabp4+ Endothelial cells_Obese", ident.2 = "Gpihbp1+ Fabp4+ Endothelial cells_Non-Obese", verbose = FALSE)
SP_MASTcells <- FindMarkers(Spleen, only.pos = F, ident.1 = "MAST cells_Obese", ident.2 = "MAST cells_Non-Obese", verbose = FALSE)

sctlist <- c("SP_NK", "SP_Lef1Tcf7Cd4Tcells", "SP_Lef1Tcf7Cd8Tcells", "SP_Cd8Ccl5Nkg7CytotoxicTcells", "SP_Foxp3RegTcells", "SP_ILC", "SP_Bcells", "SP_Plasmacells", "SP_Macro_Folr2Lyve1", "SP_Macro_Pf4Retnla", 
             "SP_Macro_infl", "SP_Mono_inter", "SP_Mono_clandnc", "SP_Conv_DC1", "SP_Migr_DC", "SP_Mki67Top2acells", "SP_s100a9_a8Granulocytes", "SP_Gpihbp1Fabp4EC", "SP_MASTcells")

for(i in 1:length(sctlist)) {                
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}

### Save whole data as H5 object ####

SaveH5Seurat(object = TERVA2_harmony, overwrite = T, verbose = T)

# Get cell labels from the TERVA data

cell_identities <- Idents(TERVA2_harmony)
write.csv(cell_identities, file = "cell_identities_TERVA2_LM.csv")

