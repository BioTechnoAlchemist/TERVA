#### TERVA2 data integration with harmony ####
setwd("/scratch/project_2005050/Rstats")

.libPaths(c("/projappl/project_2005050/project_rpackages_4.1.1", .libPaths()))

#library(harmony)
#library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(viridis)
library(gridExtra)
library(clustree)
library(RColorBrewer)
#library(EnhancedVolcano)

#devtools::install_github('satijalab/seurat-data')
#BiocManager::install("EnhancedVolcano")

#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#remotes::install_github("mojaveazure/seurat-disk")

#### Load in H5 data ####

LDAO <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDAOafterQC2.h5Seurat") #ReQC
PLAO <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLAOafterQC2.h5Seurat") #ReQC
LDeWAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/LD_eWATafterQC.h5Seurat")
PLeWAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLeWATafterQC2.h5Seurat") #ReQC
LDPVAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDPVATafterQC.h5Seurat")
PLPVAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLPVATafterQC2.h5Seurat") #ReQC
LDSpleen <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDSpleenafterQC.h5Seurat")
PLSpleen <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLSpleenafterQC.h5Seurat")

New_idents <- rep("LDAO", times = 3048)
LDAO@meta.data$Sample <- as.factor(New_idents)
New_idents <- rep("PLAO", times = 7661)
PLAO@meta.data$Sample <- as.factor(New_idents)

New_idents <- rep("LDeWAT", times = 8076)
LDeWAT@meta.data$Sample <- as.factor(New_idents)
New_idents <- rep("PLeWAT", times = 6420)
PLeWAT@meta.data$Sample <- as.factor(New_idents)

New_idents <- rep("LDPVAT", times = 7396)
LDPVAT@meta.data$Sample <- as.factor(New_idents)
New_idents <- rep("PLPVAT", times = 5396)
PLPVAT@meta.data$Sample <- as.factor(New_idents)

New_idents <- rep("LDSpleen", times = 4234)
LDSpleen@meta.data$Sample <- as.factor(New_idents)
New_idents <- rep("PLSpleen", times = 5264)
PLSpleen@meta.data$Sample <- as.factor(New_idents)

#### Merge the separate data sets ####
# merge.data=T 
TERVA2 <-  merge(LDAO, y = c(PLAO, LDeWAT, PLeWAT, LDPVAT, PLPVAT, LDSpleen, PLSpleen), add.cell.ids = c("LDAO", "PLAO", "LDeWAT", "PLeWAT", "LDPVAT", "PLPVAT", "LDSpleen", "PLSpleen"), project = "TERVA2", merge.data = T)

remove("LDAO","LDeWAT","PLeWAT","LDSpleen", "PLAO", "LDPVAT", "PLPVAT", "PLSpleen", "New_idents")

TERVA2 
#An object of class Seurat 
#32423 features across 47495 samples within 2 assays 
#Active assay: RNA (32285 features, 0 variable features)
#1 other assay present: ADT

#unique(TERVA2@meta.data$Sample)
table(TERVA2$Sample)
#LDAO   LDeWAT   LDPVAT LDSpleen     PLAO   PLeWAT   PLPVAT PLSpleen 
#3048     8076     7396     4234     7661     6420     5396     5264 

#### Normalize, scale etc. RNA assay ####
# See https://github.com/immunogenomics/harmony/issues/41 for why I decided to use NormalizeData instead of SCTransform.
# See https://portals.broadinstitute.org/harmony/SeuratV3.html. Before running Harmony, make a Seurat object and following the standard pipeline through PCA.
#IMPORTANT DIFFERENCE: In the Seurat integration tutorial, you need to define a Seurat object for each dataset. With Harmony integration, create only one Seurat object with all cells.

#no need to re-normalize as the merged objects have normalized data. Just do rescaling.
TERVA2 <- ScaleData(TERVA2, features = rownames(TERVA2), verbose = FALSE)
TERVA2 <- FindVariableFeatures(TERVA2, selection.method = "vst", nfeatures = 3000)
TERVA2 <- RunPCA(TERVA2, verbose = FALSE)
#ElbowPlot(TERVA2)

TERVA2 <- RunUMAP(TERVA2, dims = 1:20)
TERVA2 <-  FindNeighbors(TERVA2, dims = 1:20)
TERVA2 <-  FindClusters(TERVA2, resolution = 0.5) 

# Check how the non-integrated data looks after merging and normalization.

#clustersbysample <- DimPlot(TERVA2, reduction = "umap", pt.size=0.5, label = T, group.by = "Sample")
#clustersbyfindcluster <- DimPlot(TERVA2, reduction = "umap", pt.size=0.5, label = T, group.by = "RNA_snn_res.0.5")

#clustersbysample + clustersbyfindcluster



#### Normalization of ADT data ####

TERVA2 <- NormalizeData(TERVA2, assay = "ADT", normalization.method = "CLR")
TERVA2 <- ScaleData(TERVA2, assay = "ADT")
TERVA2 <- RunPCA(TERVA2, assay = "ADT", reduction.name = "apca", features = rownames(TERVA2@assays$ADT))

#saveRDS(TERVA2, "TERVA2_reQC.rds")

#### Running Harmony ####

TERVA2_harmony <- TERVA2 %>% RunHarmony("Sample", plot_convergence = T)
#saveRDS(TERVA2_harmony, "TERVA2_harmony2.rds")

#### Read in saved harmonized data ####

#Use the one further down with the name TERVA2_harmony_2.
#TERVA2_harmony <- readRDS("/scratch/project_2005050/Rstats/TERVA2_harmony.rds") #use this to skip all the preliminary steps before harmony. Md5sum checked: ok.


TERVA2_harmony <- TERVA2_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()

#### Clustree ####

TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.1)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.2)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.3)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.4)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.5) 
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.6) #This seems to be optimal resolution
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.7)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.8)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 0.9)
TERVA2_harmony <- FindClusters(TERVA2_harmony, resolution = 1)

library(clustree)
plot_cls1 <- clustree(TERVA2_harmony)
plot_cls2 <- clustree(TERVA2_harmony, node_colour= "Cd79a", node_colour_aggr= "median")
plot_cls3 <- clustree(TERVA2_harmony, node_colour= "Cd3d", node_colour_aggr= "median")
plot_cls4 <- clustree(TERVA2_harmony, node_colour= "Ighm", node_colour_aggr= "median")

plot_cls1




# Plot how integrated data sets are presented within clusters


plot_integrated_clusters(TERVA2_harmony)

#### Cell cycle analysis ####
#exp.mat <- read.table(file = "/scratch/project_2005050/Rstats/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
#                      as.is = TRUE, row.names = 1)

s.genes <- str_to_title(cc.genes$s.genes) #If this gives an error, reinstall package stringr 
g2m.genes <- str_to_title(cc.genes$g2m.genes)

TERVA2_harmony <- CellCycleScoring(object = TERVA2_harmony, g2m.features = g2m.genes, 
                                   s.features = s.genes)
TERVA2_harmony@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")


#### Check how everything looks at this point, some basic plots ####

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "Sample")
clbysample_harmony <- DimPlot(TERVA2_harmony,reduction = "umap") + plot_annotation(title = "TERVA2 samples after integration (Harmony)")

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "RNA_snn_res.0.6")
clbycl_harmony <- DimPlot(TERVA2_harmony,label = T) + NoLegend()

clbysample_harmony + clbycl_harmony

VlnPlot(TERVA2_harmony, features = c("Igkv3-4", "Igkv1-135",	"Ighv1-72",	"Igkv4-55",	"Jchain", "Ighm")) #,	"Igkv5-43",	"Ighv9-3",	"Igkv3-12",	"Igkv12-46",	"Ighv10-3",	"Igkv6-20",	"Igkv12-44",	"Igkv5-39",	"Iglv1",	"Igkv4-57",	"Ighv1-26",	"Igkv9-120",	"Iglc2",	"Igkv6-15",	"Igkv3-10",	"Igkv19-93",	"Igkv4-91",	"Igkc",	"Ighm",	"Igkv10-96",	"Ighv9-4",	"Igtp",	"Igkv12-89",	"Igkv1-99",	"Igkv4-59",	"Iglc3",	"Igf2r",	"Ighv1-78",	"Igf1",	"Ighv5-16",	"Ighv1-72",	"Igsf6",	"Igbp1",	"Ighv1-62-2",	"Igkv14-111",	"Igkv1-110",	"Iglv3",	"Ighv1-81",	"Ighg3",	"Ighg2c",	"Iglc1",	"Ighg2b",	"Ighv1-53",	"Igfbp6",	"Igsf9b",	"Ighv2-9",	"Ighv14-2",	"Igf1r",	"Igsf8",	"Igfbp-7",	"Igfbp3",	"Igfbp5",	"Igsf1",	"Igdcc4",	"Igsf3",	"Igkv16-104",	"Igha",	"Igkv4-79",	"Ighv5-9",	"Ighv9-2",	"Ighv1-50",	"Ighv1-85",	"Ighv3-6",	"Ighv6-3",	"Igfbp7",	"Iglv2",	"Igsf10",	"Ighv1-15",	"Igfbp4",	"Ighd", "Igflr1",	"Ighv1-69",	"Igkv17-121",	"Ighv1-9",	"Ighv1-18",	"Ighv1-64",	"Igkv6-17",	"Ighv1-55",	"Igkv6-25",	"Igkv17-127",	"Ighj2",	"Ighv1-80",	"Igkj2"))



FOI <- FeaturePlot(TERVA2_harmony, features= c("Slc7a7", "Lgals3", "Ly6a"), split.by = "Sample", pt.size = 1)
Fcellannot <- FeaturePlot(TERVA2_harmony, features= c("Myh11", "Nkg7", "Rora", "Gata3", "Dcn", "Vcam1", "Col2a1", "Rgs5", "Cd34", "Upk3b", "Pecam1", "Cd79a", "Cd3d", "Cd68", "Spp1", "S100a8", "S100a9", "Top2a")) 
Fadip <- FeaturePlot(TERVA2_harmony, features= c("Adipoq", "Pdgfra", "Ces1f", "Btc", "Apoe", "Cacna1a", "Prune2", "Mt2", "Tcf21"))
Fangio <- FeaturePlot(TERVA2_harmony, features= c("Cdh5", "Sdpr", "Egfl7", "Ptprb", "Ecscr", "Cldn5", "Icam2", "Slc9a3r2", "Myh11", "Cnn1", "Vcam1", "Cd36"))
TB <- FeaturePlot(TERVA2_harmony, features= c("Cd79a", "Cd3d"))
CD4_CD8 <- FeaturePlot(TERVA2_harmony, features= c("Cd8a", "CD8a", "Ccl5", "Cd4", "CD4", "Foxp3", "Icos", "Kit", "Il23r"))
Plasmacells <- FeaturePlot(TERVA2_harmony, features = c("Ighd", "IgD", "CD138-Syndecan-1", "CD45R-B220", "Tnfrsf17", "Cxcr4"))
Monocytes <- FeaturePlot(TERVA2_harmony, features = c("Ly-6C", "Cx3cr1", "Ccr2", "Sell", "Spn", "Cd209a", "I-A-I-E", "Treml4"))
DCs <- FeaturePlot(TERVA2_harmony, features = c("Ptprc", "CD11c", "I-A-I-E", "Ly-6C", "Bst2", "SiglecH", "CD45R-B220", "Cd8a", "Xcr1", "CD24", "CD370-CLEC9A-DNGR1", "Cd4", "Sirpa", "CD11b"))
Macrophages <- FeaturePlot(TERVA2_harmony, features = c("Cd163","Cd80", "Cd86", "Ccr5", "CD11b", "CD11c", "Cd14", "CD15-SSEA-1", "Cd68", "Cd36", "Ptgs2", "Irf5", "Stat1", "Nos2", "Cxcr1", "Cxcr2", "Mrc1", "Irf4", "Stat6", "Socs3", "Sphk1", "Tlr8"))
Fcrs <- FeaturePlot(TERVA2_harmony, features = c("B2m","Ero1l", "Fcamr", "Fcer1a", "Fcer1g", "Fcer2a", "Fcgrt", "Fcrla", "Fcrlb", "Fcrls", "Fcrl1", "Fcrl5", "Fcrl6", "Fcgbp", "Fcgrt", "Fcgr2b", "Fcgr4", "Fcgr3", "Fcgr1"))
Igs_lc <- FeaturePlot(TERVA2_harmony, split.by = "Sample", features = c("Igkv3-4", "Igkv4-55", "Igkv1-135",  "Jchain", "Ighm"))
Lgals3 <- FeaturePlot(TERVA2_harmony, split.by = "Sample", features = c("Lgals3") , pt.size = 1)
Rgs5 <- FeaturePlot(TERVA2_harmony, features = c("Rgs5") , pt.size = 1)


Idents(TERVA2_harmony) <- "celltype"
VlnPlot(TERVA2_harmony, features = c("Lgals3"), idents = c("Macrophages","Macrophages activated","Intermediate monocytes", "Classical and non-classical monocytes", "Lgals3+ Macrophages", "Granulocytes", "Vascular smooth muscle cells", "Fibroblasts", "Fibroblasts activated"), split.by = "Sample", cols = c("darkorange", "magenta", "cyan4", "tomato", "turquoise2", "plum3", "slateblue2", "yellow")) # + theme(legend.position = 'none')
Idents(TERVA2_harmony) <- "celltype.group"
RidgePlot(TERVA2_harmony, features = c("Lgals3"), idents = c("Macrophages_Late_disease", "Macrophages_Prelesion","Macrophages activated_Late_disease","Macrophages activated_Prelesion", "Lgals3+ Macrophages_Late_disease", "Lgals3+ Macrophages_Prelesion", "Vascular smooth muscle cells_Late_disease", "Vascular smooth muscle cells_Prelesion", "Fibroblasts_Late_disease", "Fibroblasts_Prelesion"), sort = "increasing") + NoLegend()
DoHeatmap(TERVA2_harmony, features = c("Lgals3"))

Generalmarkers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top10
heatmap <- DoHeatmap(TERVA2_harmony, features = top10$gene) + NoLegend()

FOI
Fcellannot
Fadip | Fangio
TB
CD4_CD8
Plasmacells
Monocytes
DCs
Macrophages
Fcrs
Igs_lc
Lgals3


FeaturePlot(TERVA2_harmony, features = c("CD4", "CD8a", "CD366-Tim-3", "CD279-PD-1", "CD117-c-kit", "Ly-6C", "CD11b","Ly-6G", "CD49f", "CD44", "CD54", "CD90.2", "CD15-SSEA-1", "CD73"))

molecule_counts <- VlnPlot(TERVA2_harmony, features = "nCount_RNA", group.by = "RNA_snn_res.0.6")
gene_counts <- VlnPlot(TERVA2_harmony, features = "nFeature_RNA", group.by = "RNA_snn_res.0.6")

molecule_counts + gene_counts | clbycl_harmony

# Subset the CD79a and Cd3d positive cluster. After removing this, the cell numbers per sample drop, especially in the aorta samples.
#TERVA2_harmony <- subset(TERVA2_harmony, subset = RNA_snn_res.0.6 %in% c(0:8,10:35)) #Remove cluster 9 that's CD79a and Cd3d positive and might bias subsequent workflows


# Check CD45 and eWAT markers across the data

T2_harmony_ptprcvsTcf21 <- FeaturePlot(TERVA2_harmony, features= c("Ptprc", "Tcf21"))

T2_harmony_ptprcvsTcf21 + T2_merged_ptprcvsTcf21 + 
  plot_layout(ncol = 2, widths = c(1,1))

#### Crude check of TCR/BCR cells across the data ####
library(dittoSeq)
v1 <- dittoDimPlot(TERVA2_harmony, var = "TCR_chains", size = 2,
                   cells.use = grepl("TRA", TERVA2$TCR_chains, fixed = TRUE), 
                   show.others = TRUE)

v2 <- dittoDimPlot(TERVA2_harmony, var = "BCR_isotype", size = 2,
                   cells.use = grepl("IGHM", TERVA2$BCR_isotype, fixed = TRUE), 
                   show.others = TRUE)
v1 + v2 + TB + TCRplot + BCRplot

FeaturePlot(TERVA2_harmony, features = c("CD278-ICOS-clone7E.17G9", "Icos", "CD117-c-kit", "Kit", "CD127-IL-7Ralpha", "Ilr7"))
#CD278, CD117, CD127 ab and ge data
#Il7r = CD127, Icos = CD278, Kit = CD117 

#TCR markers from the ab data                      

TCRplot <- FeaturePlot(TERVA2_harmony, features = c("TCRValpha8.3-cloneB21.14", "TCRValpha8.3-cloneKT50", "TCRValpha11.1-11.2", "TCRgammadelta", "TCRValpha2"))

# Immunoglobulins and other B cells markers + Cd3d for control

BCRplot <- FeaturePlot(TERVA2_harmony, features = c("IgG1-kappa-isotypeCtrl", "IgG2a-kappa-isotypeCtrl", "IgG2b-unknown-isotypeCtrl", "IgM", "IgD", "CD93-AA4.1-earlyBlineage"))

TCRplot
BCRplot

#### Check the Ab data ####

# Antibody data had varying success rates for the different samples.
# LD eWAT: high fraction of antibod reads in aggregate barcodes (25.80 %, ideal < 5 %)

ADTs <- as.list(TERVA2_harmony@assays$ADT@counts@Dimnames)
ADTnames <- ADTs[[1]]

for(i in 1:length(ADTnames)){ #Print and save (138) ADT images to a specific folder
  FeaturePlot(TERVA2_harmony, ADTnames[i])
  ggsave(paste0("/scratch/project_2005050/Rstats/", i, ".jpg"))
  while (!is.null(dev.list()))  dev.off()
}

#### Predicting cell types from reference ####

# Tiit's cell types (+ other markers)

#At a very crude cell type level, can see we have cluster of:
#B cell (Cd79a) (Cd97b, Ebf1, Ms4a1, Ly-6d, Fcmr, Siglecg,H2-dmb2, Vpreb3)
#T cell (Cd3d) (mem T cells: Cd3g, Lck, Thy1, Ms4a4b, Cd8b1, Ccl5)
#Natural killer (Nkg7 pos, Cd3d neg)
#ILC (Rora, Gata3)
#macrophage/monocyte/myeloid-dendritic (Cd68) (monocytes; S100a9, S100a8, Retnlg, Msrb1, Slpi, Wfdc21, Il1b, G0s2, Lcn2, Cxcl2, macrophages; Retnla, Mt1, C1qb, Fcrls, C1qc, C1qa, Pf4, Sepp1, Lpl, Lyz1; Gngt2, Cst3, Clec4a3, Ifitm3, Plac8, Hes1, Ms4a6c, Lyz2,Atf3, Ear2)
#endothelial (Pecam1)
#smooth muscle (Myh11)
#Pericyte (Rgs5)
#Possibly neutrophil (S100a9, Csf3r)
#Fibroblasts and pre-adipocytes (Dcn, Cd34)
#Epithelial (Upk3b)
#Dividing cells of various types (Top2a)
#(some other small ones as well, like basophils/mast cells)
# principal hematopoietic markers (Cd3e, Adgre1, Ly-6c, Cd19, Cd68, H2-ab1)


#BiocManager::install("SingleR")
#BiocManager::install("celldex")

library(celldex)
library(SingleR)

#mouse.ref <- celldex::MouseRNAseqData()
#saveRDS(mouse.ref, "mouseref.rds")
mouse.ref <- readRDS("/scratch/project_2005050/Rstats/mouseref.rds")

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

#### Rename unannotated clusters ####

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "seurat_clusters")
TERVA2_harmony <- RenameIdents(TERVA2_harmony, 
                               "0" = "B cells", "1" = "Fibroblasts", "2" = "Macrophages", "3" = "Cd8+ Tcells", 
                               "4" = "Cd4+ Tcells", "5" = "Endothelial cells", "6" = "Macrophages activated", "7" = "Intermediate Monocytes", 
                               "8" = "Classical and non-classical Monocytes", "10" = "Pi16+ adventitial Fibroblasts", "11" = "Cd4+ Foxp3+ Tregs", 
                               "12" = "Vascular smooth muscle cells", "13" = "Fibroblasts activated", "14" = "Cd8+ Ccl5+ Teffs", 
                               "15" = "Innate lymphoid cells", "16" = "Lgals3+ Macrophages", 
                               "17" = "Rpl-Rps+ cells", "18" = "Granulocytes", "19" = "Dividing cells",
                               "20" = "Rgs5+ EC's", "21" = "Conventional Dendritic cells 1", "22" = "Natural killer cells", 
                               "23" = "Conventional Dendritic cells 2", "24" = "Plasma cells", "25" = "Pi16+ adventitial Fibroblasts",
                               "26" = "EC's with shear stress markers", "27" = "Mgp+ Fibroblasts", 
                               "28" = "Hdc+ Cpa3+ MAST cells", "29" = "Mgp+ Fibroblasts", "30" = "Adipose tissue Mesothelial cells",
                               "31" = "Endothelial cells", "32" = "Igf2+ Endothelial tip cells", "33" = "Mmp3+ Adipocytes",
                               "34" = "Fibroblasts", "35" = "Granulocytes")

nb.cols <- 31
mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)
reannot <- DimPlot(TERVA2_harmony, label = T , repel = T, pt.size = 3, label.size = 9, cols = mycolors)
reannot + 
  labs(title = "Clustering of the integrated data") +
  theme(plot.title = element_text(size=30, hjust= 0.5),
        axis.title = element_text(size=20),
        axis.text = element_text(size=20)) +
  guides(col = guide_legend(ncol = 1)) + NoLegend()


  

saveRDS(TERVA2_harmony, "TERVA2_harmony_2.rds")
TERVA2_harmony <- readRDS("/scratch/project_2005050/Rstats/TERVA2_harmony_2.rds")

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "seurat_clusters")
TERVA2_harmony_allmarkers_wilcoxon <- FindAllMarkers(TERVA2_harmony, verbose = T, min.cells.group = 10)
saveRDS(TERVA2_harmony_allmarkers_wilcoxon, "TERVA2_harmony_allmarkers_wilcoxon_reQC.rds")


#### LD vs PL analysis ####

coldata <- as.data.frame(TERVA2_harmony@meta.data$Sample)

coldata <- mutate(coldata,
                  Status = case_when(
                    startsWith(TERVA2_harmony@meta.data$Sample, "LD") ~ "Late_disease",
                    startsWith(TERVA2_harmony@meta.data$Sample, "PL") ~ "Prelesion"
                  ))

group_id <- coldata$Status
TERVA2_harmony@meta.data$group_id <- as.factor(group_id)

#Plot by group

DimPlot(TERVA2_harmony, reduction = "umap", group.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()
DimPlot(TERVA2_harmony, reduction = "umap", split.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()

FeaturePlot(TERVA2_harmony, split.by = "group_id", features = "Rag1")
FeaturePlot(TERVA2_harmony, split.by = "group_id", features = c("Icos", "Kit", "Il-r7"))

#Find conserved markers for cells of interest between PL and LD (This is across tissues)

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "seurat_clusters") #When all clusters separately 
B.markers <- FindConservedMarkers(TERVA2_harmony, ident.1 = c(0,25), grouping.var = "group_id", verbose = TRUE)

TERVA2_harmony <- SetIdent(TERVA2_harmony, value = "mouse.fine") #When mouse reference annotated
B.mousemarkers <- FindConservedMarkers(TERVA2_harmony, ident.1 = "B cells", grouping.var = "group_id", verbose = TRUE)

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "seurat_clusters") #When all clusters separately 
T.markers <- FindConservedMarkers(TERVA2_harmony, ident.1 = c(2,5,9,11,23,24), grouping.var = "group_id", verbose = TRUE)

TERVA2_harmony <- SetIdent(TERVA2_harmony, value = "mouse.fine") #When mouse reference annotated
T.mousemarkers <- FindConservedMarkers(TERVA2_harmony, ident.1 = "T cells", grouping.var = "group_id", verbose = TRUE)

TERVA2_harmony <- SetIdent(TERVA2_harmony, value = "mouse.fine") #When mouse reference annotated
Macrophage.mousemarkers <- FindConservedMarkers(TERVA2_harmony, ident.1 = "Macrophages", grouping.var = "group_id", verbose = TRUE)
Macrophage_act.mousemarkers <- FindConservedMarkers(TERVA2_harmony, ident.1 = "Macrophages activated", grouping.var = "group_id", verbose = TRUE)



#Find de genes between PL and LD (across tissues)
#TERVA2_harmony <- subset(TERVA2_harmony, subset = seurat_clusters %in% c(0:9,11:50)) #Remove the CD79a + Cd3d pos cluster
#TERVA2_harmony <- SetIdent(TERVA2_harmony, value = "mouse.fine")
TERVA2_harmony$celltype.group <- paste(Idents(TERVA2_harmony), TERVA2_harmony$group_id, sep = "_")
TERVA2_harmony$celltype <- Idents(TERVA2_harmony)
celltypes <- as.vector(unique(TERVA2_harmony$celltype))
celltypes
#[1] "Cd4+ Tcells"                           "Natural killer cells"                  "Cd8+ Tcells"                           "Cd8+ Ccl5+ Teffs"                     
#[5] "Cd4+ Foxp3+ Tregs"                     "B cells"                               "Classical and non-classical Monocytes" "Lgals3+ Macrophages"                  
#[9] "Plasma cells"                          "Dividing cells"                        "Vascular smooth muscle cells"          "Intermediate Monocytes"               
#[13] "Macrophages activated"                 "Innate lymphoid cells"                 "Rpl-Rps+ cells"                        "Pi16+ adventitial Fibroblasts"        
#[17] "Mgp+ Fibroblasts"                      "Macrophages"                           "Conventional Dendritic cells 2"        "Fibroblasts"                          
#[21] "Granulocytes"                          "Conventional Dendritic cells 1"        "Hdc+ Cpa3+ MAST cells"                 "EC's with shear stress markers"       
#[25] "Fibroblasts activated"                 "Rgs5+ EC's"                            "Adipose tissue Mesothelial cells"      "Endothelial cells"                    
#[29] "Mmp3+ Adipocytes"                      "Igf2+ Endothelial tip cells"  

####  DE from the finely reannotated cell types DE between states ####
Idents(TERVA2_harmony) <- "celltype"
Allmarkers_macrophages_wilcoxon <- FindMarkers(TERVA2_harmony, ident.1 = "Lgals3+ Macrophages", ident.2 = c("Macrophages", "Macrophages activated"), verbose = T, min.cells.group = 30)
saveRDS(TERVA2harmony_allmarkers_wilcoxon, "TERVA2harmony_allmarkers_Macrophages_wilcoxon_reclustered.rds")


Idents(TERVA2_harmony) <- "celltype.group"
Cd4posT <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Cd4+ Tcells_Late_disease", ident.2 = "Cd4+ Tcells_Prelesion", verbose = FALSE)
NK <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Natural killer cells_Late_disease", ident.2 = "Natural killer cells_Prelesion", verbose = FALSE)
Cd8posT <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Cd8+ Tcells_Late_disease", ident.2 = "Cd8+ Tcells_Prelesion", verbose = FALSE)
Cd8posCcl5posTeff <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Cd8+ Ccl5+ Teffs_Late_disease", ident.2 = "Cd8+ Ccl5+ Teffs_Prelesion", verbose = FALSE)
Cd4Foposxp3posTreg <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Cd4+ Foxp3+ Tregs_Late_disease", ident.2 = "Cd4+ Foxp3+ Tregs_Prelesion", verbose = FALSE)
Bcells_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "B cells_Late_disease", ident.2 = "B cells_Prelesion", verbose = FALSE)
Macrophagesact_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Macrophages activated_Late_disease", ident.2 = "Macrophages activated_Prelesion", verbose = FALSE)
Plasmacells_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Plasma cells_Late_disease", ident.2 = "Plasma cells_Prelesion", verbose = FALSE)
Dividing <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Dividing cells_Late_disease", ident.2 = "Dividing cells_Prelesion", verbose = FALSE)
Fibro_Pi16 <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Pi16+ adventitial Fibroblasts_Late_disease", ident.2 = "Pi16+ adventitial Fibroblasts_Prelesion", verbose = FALSE)
Fibroblasts_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Fibroblasts_Late_disease", ident.2 = "Fibroblasts_Prelesion", verbose = FALSE)
Fibro_Mgp_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Mgp+ Fibroblasts_Late_disease", ident.2 = "Mgp+ Fibroblasts_Prelesion", verbose = FALSE)
Fibroact_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Fibroblasts activated_Late_disease", ident.2 = "Fibroblasts activated_Prelesion", verbose = FALSE)
Macrophages_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Macrophages_Late_disease", ident.2 = "Macrophages_Prelesion", verbose = FALSE)
Macrophages_Lgals3_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Lgals3+ Macrophages_Late_disease", ident.2 = "Lgals3+ Macrophages_Prelesion", verbose = FALSE)
Int_monocytes <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Intermediate monocytes_Late_disease", ident.2 = "Intermediate monocytes_Prelesion", verbose = FALSE)
ILC <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Innate lymphoid cells_Late_disease", ident.2 = "Innate lymphoid cells_Prelesion", verbose = FALSE)
Conv_DC2 <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Conventional Dendritic cells 2_Late_disease", ident.2 = "Conventional Dendritic cells 2_Prelesion", verbose = FALSE)
Class_monocytes <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Classical and non-classical Monocytes_Late_disease", ident.2 = "Classical and non-classical Monocytes_Prelesion", verbose = FALSE)
Granulocytes_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Granulocytes_Late_disease", ident.2 = "Granulocytes_Prelesion", verbose = FALSE)
Conv_DC1 <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Conventional Dendritic cells 1_Late_disease", ident.2 = "Conventional Dendritic cells 1_Prelesion", verbose = FALSE)
Endothelialcells_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Endothelial cells_Late_disease", ident.2 = "Endothelial cells_Prelesion", verbose = FALSE)
EC_rgs5_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Rgs5+ EC's_Late_disease", ident.2 = "Rgs5+ EC's_Prelesion", verbose = FALSE)
VSMCs_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Vascular smooth muscle cells_Late_disease", ident.2 = "Vascular smooth muscle cells_Prelesion", verbose = FALSE)
MAST_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Hdc+ Cpa3+ MAST cells_Late_disease", ident.2 = "Hdc+ Cpa3+ MAST cells_Prelesion", verbose = FALSE)
ECtip_Igf2_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Igf2+ Endothelial tip cells_Late_disease", ident.2 = "Igf2+ Endothelial tip cells_Prelesion", verbose = FALSE)
Mmp3_Adipocytes <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Mmp3+ Adipocytes_Late_disease", ident.2 = "Mmp3+ Adipocytes_Prelesion", verbose = FALSE)
Granulocytes_reannot <- FindMarkers(TERVA2_harmony, only.pos = F, ident.1 = "Granulocytes_Late_disease", ident.2 = "Granulocytes_Prelesion", verbose = FALSE)


sctlist <- c("Cd4posT", "Cd8posT", "NK", "Cd8posCcl5posTeff", "Cd4Foposxp3posTreg", "Bcells_reannot", "Macrophagesact_reannot", "Plasmacells_reannot", "Dividing", "Fibro_Pi16", 
             "Fibroblasts_reannot", "Fibro_Mgp_reannot", "Fibroact_reannot", "Macrophages_reannot", "Macrophages_Lgals3_reannot", "Int_monocytes", "ILC", "Conv_DC2", "Class_monocytes", "Granulocytes_reannot", "Conv_DC1", 
             "Endothelialcells_reannot", "EC_rgs5_reannot", "VSMCs_reannot", "MAST_reannot")

for(i in 1:length(sctlist)) {                
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}





#### Tissue specific analyses ####

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

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "tissue_id") #When all clusters separately 

AO <- subset(TERVA2_harmony, idents = "Aorta")
PVAT <- subset(TERVA2_harmony, idents = "PVAT")
eWAT <- subset(TERVA2_harmony, idents = "eWAT")
Spleen <- subset(TERVA2_harmony, idents = "Spleen")

# FindAllMarkers per tissue type ####

               
AO <- SetIdent(AO,value = "celltype")               
AO_allmarkers_macrophages_wilcoxon <- FindMarkers(AO, ident.1 = "Lgals3+ Macrophages", ident.2 = c("Macrophages", "Macrophages activated"), verbose = T, min.cells.group = 30)
saveRDS(AO_allmarkers_wilcoxon, "AO_allmarkers_Macrophages_wilcoxon_reclustered.rds")

PVAT <- SetIdent(PVAT,value = "seurat_clusters")
PVAT_allmarkers_wilcoxon <- FindAllMarkers(PVAT, verbose = T, min.cells.group = 10)
saveRDS(PVAT_allmarkers_wilcoxon, "PVAT_allmarkers_wilcoxon_reQC.rds")

eWAT <- SetIdent(eWAT,value = "seurat_clusters")
eWAT_allmarkers_wilcoxon <- FindAllMarkers(eWAT, verbose = T, min.cells.group = 10)
saveRDS(eWAT_allmarkers_wilcoxon, "eWAT_allmarkers_wilcoxon_reQC.rds")

Spleen <- SetIdent(Spleen,value = "seurat_clusters")
Spleen_allmarkers_wilcoxon <- FindAllMarkers(Spleen, verbose = T, min.cells.group = 10)
saveRDS(Spleen_allmarkers_wilcoxon, "Spleen_allmarkers_wilcoxon_reQC.rds")


#### Aorta ####

Idents(AO) <- "celltype.group"
DimPlot(AO, label = T , repel = T, label.size = 3) + NoLegend()
DimPlot(AO, reduction = "umap", group.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()
FeaturePlot(AO, features = c("Cd3d", "Cd79a", "Igkv3-4"))
p1 <- DimPlot(AO, reduction = "umap", split.by = "group_id", group.by = "seurat_clusters", pt.size = 1.2, label = T , repel = T, label.size = 4) + NoLegend() + ggtitle("Aorta")
l1 <- FeaturePlot(AO, features = c("Ly6a", "Ly-6A-E-Sca-1")) +
  ggtitle("Ly-6A-E-Sca-1", subtitle = "Aorta")

Idents(AO) <- "celltype"
macros_AO <- subset(AO, idents = c("Macrophages", "Macrophages activated", "Lgals3+ Macrophages"))
Idents(macros_AO) <- "group_id"
FeaturePlot(macros_AO, split.by = "group_id", features = c("Trem2", "Lgals3"), blend = TRUE, cols = c("navy", "darkgoldenrod1"))
FeaturePlot(macros_AO, split.by = "group_id", features = c("Lamp2", "Mrc1"), blend = TRUE, cols = c("navy", "darkgoldenrod1"))
RidgePlot(macros_AO, group.by = "group_id", features = c("Trem2", "Lgals3"))

#AO[["percent.ig"]] <- PercentageFeatureSet(AO, pattern = "^Ig-") +  test.use= "LR", latent.vars = c("percent.mt", "percent.rb","percent.ig") in FindMarkers did not change anything

# AO LD vs PL DE's per cell type ####
Idents(AO) <- "celltype.group"
AO_Cd4posT <- FindMarkers(AO, only.pos = F, ident.1 = "Cd4+ Tcells_Late_disease", ident.2 = "Cd4+ Tcells_Prelesion", verbose = FALSE)
AO_NK <- FindMarkers(AO, only.pos = F, ident.1 = "Natural killer cells_Late_disease", ident.2 = "Natural killer cells_Prelesion", verbose = FALSE)
AO_Cd8posT <- FindMarkers(AO, only.pos = F, ident.1 = "Cd8+ Tcells_Late_disease", ident.2 = "Cd8+ Tcells_Prelesion", verbose = FALSE)
AO_Cd8posCcl5posTeff <- FindMarkers(AO, only.pos = F, ident.1 = "Cd8+ Ccl5+ Teffs_Late_disease", ident.2 = "Cd8+ Ccl5+ Teffs_Prelesion", verbose = FALSE)
AO_Cd4Foposxp3posTreg <- FindMarkers(AO, only.pos = F, ident.1 = "Cd4+ Foxp3+ Tregs_Late_disease", ident.2 = "Cd4+ Foxp3+ Tregs_Prelesion", verbose = FALSE)
AO_Bcells_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "B cells_Late_disease", ident.2 = "B cells_Prelesion", verbose = FALSE)
AO_Macrophagesact_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Macrophages activated_Late_disease", ident.2 = "Macrophages activated_Prelesion", verbose = FALSE)
AO_Plasmacells_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Plasma cells_Late_disease", ident.2 = "Plasma cells_Prelesion", verbose = FALSE)
AO_Dividing <- FindMarkers(AO, only.pos = F, ident.1 = "Dividing cells_Late_disease", ident.2 = "Dividing cells_Prelesion", verbose = FALSE)
AO_Fibro_Pi16 <- FindMarkers(AO, only.pos = F, ident.1 = "Pi16+ adventitial Fibroblasts_Late_disease", ident.2 = "Pi16+ adventitial Fibroblasts_Prelesion", verbose = FALSE)
AO_Fibroblasts_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Fibroblasts_Late_disease", ident.2 = "Fibroblasts_Prelesion", verbose = FALSE)
AO_Fibro_Mgp_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Mgp+ Fibroblasts_Late_disease", ident.2 = "Mgp+ Fibroblasts_Prelesion", verbose = FALSE)
AO_Fibroact_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Fibroblasts activated_Late_disease", ident.2 = "Fibroblasts activated_Prelesion", verbose = FALSE)
AO_Macrophages_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Macrophages_Late_disease", ident.2 = "Macrophages_Prelesion", verbose = FALSE)
AO_Macrophages_Lgals3_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Lgals3+ Macrophages_Late_disease", ident.2 = "Lgals3+ Macrophages_Prelesion", verbose = FALSE)
AO_Int_monocytes <- FindMarkers(AO, only.pos = F, ident.1 = "Intermediate Monocytes_Late_disease", ident.2 = "Intermediate Monocytes_Prelesion", verbose = FALSE)
AO_ILC <- FindMarkers(AO, only.pos = F, ident.1 = "Innate lymphoid cells_Late_disease", ident.2 = "Innate lymphoid cells_Prelesion", verbose = FALSE)
AO_Conv_DC2 <- FindMarkers(AO, only.pos = F, ident.1 = "Conventional Dendritic cells 2_Late_disease", ident.2 = "Conventional Dendritic cells 2_Prelesion", verbose = FALSE)
AO_Class_monocytes <- FindMarkers(AO, only.pos = F, ident.1 = "Classical and non-classical Monocytes_Late_disease", ident.2 = "Classical and non-classical Monocytes_Prelesion", verbose = FALSE)
AO_Granulocytes_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Granulocytes_Late_disease", ident.2 = "Granulocytes_Prelesion", verbose = FALSE)
AO_Conv_DC1 <- FindMarkers(AO, only.pos = F, ident.1 = "Conventional Dendritic cells 1_Late_disease", ident.2 = "Conventional Dendritic cells 1_Prelesion", verbose = FALSE)
AO_EC_rgs5_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Rgs5+ EC's_Late_disease", ident.2 = "Rgs5+ EC's_Prelesion", verbose = FALSE)
AO_VSMCs_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Vascular smooth muscle cells_Late_disease", ident.2 = "Vascular smooth muscle cells_Prelesion", verbose = FALSE)
AO_MAST_reannot <- FindMarkers(AO, only.pos = F, ident.1 = "Hdc+ Cpa3+ MAST cells_Late_disease", ident.2 = "Hdc+ Cpa3+ MAST cells_Prelesion", verbose = FALSE)

sctlist <- c("AO_Cd4posT", "AO_Cd8posT", "AO_NK", "AO_Cd8posCcl5posTeff", "AO_Cd4Foposxp3posTreg", "AO_Bcells_reannot", "AO_Macrophagesact_reannot", "AO_Plasmacells_reannot", "AO_Dividing", "AO_Fibro_Pi16", 
             "AO_Fibroblasts_reannot", "AO_Fibro_Mgp_reannot", "AO_Fibroact_reannot", "AO_Macrophages_reannot", "AO_Macrophages_Lgals3_reannot", "AO_Int_monocytes", "AO_ILC", "AO_Conv_DC2", "AO_Class_monocytes", 
             "AO_Granulocytes_reannot", "AO_Conv_DC1", "AO_EC_rgs5_reannot", "AO_VSMCs_reannot", "AO_MAST_reannot")

for(i in 1:length(sctlist)) {                
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}

AO_Macrophages_reannot$gene <- rownames(AO_Macrophages_reannot)
AO_Macrophages_reannot_sub <- filter(AO_Macrophages_reannot, !str_detect(gene, "^Ig")) %>%
  filter(!str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "Gm8797")) %>%
  filter(!str_detect(gene, "Jchain"))

AO_ILC$gene <- rownames(AO_ILC)
AO_ILC_sub <- filter(AO_ILC, !str_detect(gene, "^Ig")) %>%
  filter(!str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "Gm8797")) %>%
  filter(!str_detect(gene, "Jchain"))

AO_Conv_DC2$gene <- rownames(AO_Conv_DC2)
AO_Conv_DC2_sub <- filter(AO_Conv_DC2, !str_detect(gene, "^Ig")) %>%
  filter(!str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "Gm8797")) %>%
  filter(!str_detect(gene, "Jchain"))

AO_VSMCs_reannot$gene <- rownames(AO_VSMCs_reannot)
AO_VSMCs_reannot_sub <- filter(AO_VSMCs_reannot, !str_detect(gene, "^mt")) %>%
  filter(!str_detect(gene, "^Gm"))

AO_Cd8posT$gene <- rownames(AO_Cd8posT)
AO_Cd8posT_sub <- filter(AO_Cd8posT, !str_detect(gene, "^Ig")) %>%
  filter(!str_detect(gene, "Jchain")) %>%
  filter(!str_detect(gene, "Gm8797"))

AO_Cd4posT$gene <- rownames(AO_Cd4posT)
AO_Cd4posT_sub <- filter(AO_Cd4posT, !str_detect(gene, "^Ig")) %>%
  filter(!str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "^Rp")) %>%
  filter(!str_detect(gene, "Jchain"))
  

DimPlot(AO, reduction = "umap", label = T , repel = T, label.size = 3) + NoLegend()

top10 <- AO_Cd4posT_sub %>% filter(AO_Cd4posT_sub$p_val_adj < 5e-2) %>% top_n(n = 10, wt = avg_log2FC) %>% rownames()
DoHeatmap(AO, features = top10)

cells <- c("Hdc+ Cpa3+ MAST cells", "Mgp+ Fibroblasts EC's with shear stress markers", "Plasma cells", "Conv DC2", "NK's", "Conv DC1", "Ras5+ EC's", 
"Dividing cells", "Granulocytes", "Rpl-Rps+ cells", "Lgals3+ Macrophages", "ILC's", "Cd8+ Ccl5+ Teffs", "Fibroblasts activated", "VSMC's", 
"Cd4+ Foxp3+ Tregs", "Pi16+ Fibroblasts", "Classical and non-classical Monocytes", "Intermediate monocytes", "Macrophages activated",
"EC's", "Cd4+ Tcells", "Cd8+ Tcells", "Macrophages", "Fibroblasts", "B cells")

topgenes <- c("Ccl4", "Cd74", "H2-Aa", "Grn", "Lgals3", "Ly6a", "Jun", "Fos", "Apoe", "Lyz2", "Grn", "Myc", "Trem2", "Lcp1")
topgenes <- unique(topgenes)
Idents(AO) <- "celltype"
D1 <- DotPlot(AO, features = topcells, idents = cells) + 
  RotatedAxis() +
  ggtitle("Aorta")

V1<-EnhancedVolcano(AO_Bcells_reannot,
                lab = rownames(AO_Bcells_reannot),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Late disease vs Prelesion in aorta B cells',
                pCutoff = 5e-2,
                FCcutoff = 0.4,
                pointSize = 3.0,
                col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                labSize = 6.0,
                titleLabSize = 15.0)

V2<-EnhancedVolcano(AO_Macrophages_Lgals3_reannot,
                    lab = rownames(AO_Macrophages_Lgals3_reannot),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in aorta Lgals3+ macrophages',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V3<-EnhancedVolcano(AO_Macrophages_reannot_sub,
                    lab = rownames(AO_Macrophages_reannot_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in aorta macrophages',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V4<-EnhancedVolcano(AO_ILC_sub,
                    lab = rownames(AO_ILC_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in aorta ILCs',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V5<-EnhancedVolcano(AO_Conv_DC2_sub,
                    lab = rownames(AO_Conv_DC2_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in aorta conventional DC2s',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V6<-EnhancedVolcano(AO_Cd8posT_sub,
                    lab = rownames(AO_Cd8posT_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in aorta Cd8+ T cells',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V4<-EnhancedVolcano(AO_Cd4posT_sub,
                              lab = rownames(AO_Cd4posT_sub),
                              x = 'avg_log2FC',
                              y = 'p_val_adj',
                              title = 'Late disease vs Prelesion in aorta Cd4+ T cells',
                              pCutoff = 5e-2,
                              FCcutoff = 0.4,
                              pointSize = 3.0,
                              col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                              labSize = 6.0,
                              titleLabSize = 15.0)

grid.arrange(V1,V2,V3,V4,V5,V6, ncol=3)

#### eWAT ####
Idents(eWAT) <- "celltype.group"
DimPlot(eWAT, reduction = "umap", group.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()
p2 <- DimPlot(eWAT, reduction = "umap", split.by = "group_id", group.by = "seurat_clusters", pt.size = 1.2, label = T , repel = T, label.size = 4) + NoLegend() + ggtitle("eWAT")
l2 <- FeaturePlot(eWAT, features = c("Ly6a", "Ly-6A-E-Sca-1", "Lgals3")) #+
  #ggtitle("Ly-6A-E-Sca-1", subtitle = "eWAT")

Idents(eWAT) <- "celltype"
macros_eWAT <- subset(eWAT, idents = c("Macrophages", "Macrophages activated", "Lgals3+ Macrophages"))
Idents(macros_eWAT) <- "group_id"
FeaturePlot(macros_eWAT, split.by = "group_id", features = c("Trem2", "Lgals3"), blend = TRUE, cols = c("navy", "darkgoldenrod1"))
FeaturePlot(macros_eWAT, split.by = "group_id", features = c("Lamp2", "Mrc1"), blend = TRUE, cols = c("navy", "darkgoldenrod1"))
RidgePlot(macros_eWAT, group.by = "group_id", features = c("Trem2", "Lgals3"))

Idents(eWAT) <- "celltype"
eWAT_ECs <- subset(eWAT, idents = c("Endothelial cells", "EC's with shear stress markers", "Rgs5+ EC's"))
eWAT_ECs <- RunPCA(eWAT_ECs, verbose = FALSE)
ElbowPlot(eWAT_ECs)

eWAT_ECs <- eWAT_ECs %>% 
  RunUMAP(dims = 1:10) %>%
  FindNeighbors(dims = 1:10)

resolution.range <- seq(from = 0, to = 1, by = 0.1)
eWAT_ECs <- FindClusters(eWAT_ECs, resolution = resolution.range)
clustree(eWAT_ECs)
DimPlot(eWAT_ECs, reduction = "umap", group.by = "RNA_snn_res.0.4", pt.size = 1.2, label = T , repel = T, label.size = 4) + NoLegend() + ggtitle("eWAT EC's reclustered")
Idents(eWAT_ECs) <- "RNA_snn_res.0.2"
eWAT_EC_markers_res0.2 <- FindAllMarkers(eWAT_ECs)
Idents(eWAT_ECs) <- "RNA_snn_res.0.1"
eWAT_EC_markers_res0.1 <- FindAllMarkers(eWAT_ECs)
Idents(eWAT_ECs) <- "RNA_snn_res.0.3"
eWAT_EC_markers_res0.3 <- FindAllMarkers(eWAT_ECs)
Idents(eWAT_ECs) <- "RNA_snn_res.0.4"
eWAT_EC_markers_res0.4 <- FindAllMarkers(eWAT_ECs)


Pericytes <- FeaturePlot(eWAT_ECs, features = c("Cacna1c", "Myo1b","Stac","Notch3", "Kcnj8", "Abcc9", "Cspg4", "Rgs5"))
ECs <- FeaturePlot(eWAT_ECs, features = c("Cldn5","Cdh5","Pecam1", "Jam2","Mecom","Ptprb", "Cdh13", "Acta2", "Myh11"))
Pericytes | ECs

eWAT_EC_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10 
DoHeatmap(eWAT_ECs, features = top10$gene) + NoLegend()

eWAT_EC_markers_res0.1 %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10 
DoHeatmap(eWAT_ECs, features = top10$gene) + NoLegend()

VlnPlot(eWAT_ECs, features = "nFeature_RNA", group.by = "RNA_snn_res.0.2")

#### Mouse WAT atlas ####

WAT_meta <- read_tsv("/scratch/project_2005050/Rstats/mouse_WAT_atlas/metadata.tsv")
WAT_atlas <- Read10X("/scratch/project_2005050/Rstats/mouse_WAT_atlas/")

WAT_atlas <- AddMetaData(WAT_atlas, WAT_meta$cell_type__custom, col.name = "celltype")
WAT_atlas <- AddMetaData(WAT_atlas, WAT_meta$cell_subtype__custom, col.name = "celltype_sub")
unique(WAT_atlas@meta.data$celltype)

#The data set is massive and Puhti will run out of memory when scaling so subset here
Idents(WAT_atlas) <- "celltype"
WAT_sub <- subset(WAT_atlas, idents = c("endothelial", "pericyte", "male_epithelial"))

Idents(WAT_sub) <- "orig.ident"
WAT_sub <- WAT_sub %>% 
  NormalizeData() %>%
  ScaleData() %>%
  FindVariableFeatures(nfeatures = 3000, selection.method = "vst")
WAT_sub <- RunPCA(WAT_sub)
ElbowPlot(WAT_sub)

WAT_sub <- RunUMAP(WAT_sub, reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(dims = 1:20)
resolution.range <- seq(from = 0, to = 1, by = 0.1)
WAT_sub <- FindClusters(WAT_sub, resolution = resolution.range)
clustree(WAT_sub) #Based on this I'd use res 0.3

WAT_sub <- FindClusters(WAT_sub, res = 0.3)

WAT_clusters <- DimPlot(WAT_sub, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.2, label = T , repel = T, label.size = 4) + ggtitle("WAT atlas")
WAT_celltypes <- DimPlot(WAT_sub, reduction = "umap", group.by = "celltype", pt.size = 1.2, label = T , repel = T, label.size = 4) + ggtitle("WAT atlas")
WAT_subtypes <- DimPlot(WAT_sub, reduction = "umap", group.by = "celltype_sub", pt.size = 1.2, label = T , repel = T, label.size = 4) + ggtitle("WAT atlas")
WAT_clusters | WAT_celltypes | WAT_subtypes


WAT_atlas.anchors <- FindTransferAnchors(reference = WAT_sub, query = eWAT_ECs, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = WAT_atlas.anchors, refdata = WAT_sub$celltype_sub, 
                            dims = 1:30)
eWAT.query <- AddMetaData(eWAT_ECs, metadata = predictions)

eWAT.query <- FindClusters(eWAT.query, resolution = 0.4)

Idents(eWAT.query) <- "RNA_snn_res.0.2"
orig <- DimPlot(eWAT.query, reduction = "umap", group.by = "RNA_snn_res.0.2", pt.size = 1.2, label = T , repel = T, label.size = 4) + ggtitle("eWAT EC's reclustered")
Idents(eWAT.query) <- "predicted.id"
celltype <- DimPlot(eWAT.query, reduction = "umap", group.by = "predicted.id", pt.size = 1.2, label = T , repel = T, label.size = 4)  + ggtitle("eWAT EC's reclustered")
Idents(eWAT.query) <- "celltype"
celltype_orig <- DimPlot(eWAT.query, reduction = "umap", group.by = "celltype", pt.size = 1.2, label = T , repel = T, label.size = 4)  + ggtitle("eWAT EC's reclustered")


orig | celltype | celltype_orig

## With eWAT EC data that hasn't been reclustered

WAT_atlas.anchors <- FindTransferAnchors(reference = WAT_sub, query = eWAT_ECs, 
                                         dims = 1:30)
predictions <- TransferData(anchorset = WAT_atlas.anchors, refdata = WAT_sub$celltype_sub, 
                            dims = 1:30)
eWAT_ECs.query <- AddMetaData(eWAT_ECs, metadata = predictions)

Idents(eWAT_ECs.query) <- "predicted.id"
celltype_T <- DimPlot(eWAT_ECs.query, reduction = "umap", group.by = "predicted.id", pt.size = 1.2, label = T , repel = T, label.size = 4)  + ggtitle("eWAT EC's")

celltype_T








####

Idents(eWAT) <- "celltype.group"
eWAT_Tcellresponse <- FindMarkers(eWAT, only.pos = T, ident.1 = "T cells_Late_disease", ident.2 = "T cells_Prelesion", verbose = FALSE)

eWAT_Bcellresponse <- FindMarkers(eWAT, only.pos = T, ident.1 = "B cells_Late_disease", ident.2 = "B cells_Prelesion", verbose = FALSE)

eWAT_NKcresponse <- FindMarkers(eWAT, only.pos = T, ident.1 = "NK cells_Late_disease", ident.2 = "NK cells_Prelesion", verbose = FALSE)

eWAT_Macrophageresponse <- FindMarkers(eWAT, only.pos = T, ident.1 = "Macrophages activated_Late_disease", ident.2 = "Macrophages activated_Prelesion", verbose = FALSE)
eWAT_Macrophageresponse2 <- FindMarkers(eWAT, only.pos = T, ident.1 = "Macrophages_Late_disease", ident.2 = "Macrophages_Prelesion", verbose = FALSE)

eWAT_Monocyteresponse <- FindMarkers(eWAT, only.pos = T, ident.1 = "Monocytes_Late_disease", ident.2 = "Monocytes_Prelesion", verbose = FALSE)

eWAT_Granulocyteresponse <- FindMarkers(eWAT, only.pos = T, ident.1 = "Granulocytes_Late_disease", ident.2 = "Granulocytes_Prelesion", verbose = FALSE)

FeaturePlot(eWAT, features = c("Slc7a7", "Trem2", "Ly6a"), split.by = "group_id") #Rag1 not expressed in eWAT cells 

Idents(eWAT) <- "celltype.group"
eWAT_Cd4posT <- FindMarkers(eWAT, only.pos = F, ident.1 = "Cd4+ Tcells_Late_disease", ident.2 = "Cd4+ Tcells_Prelesion", verbose = FALSE)
eWAT_NK <- FindMarkers(eWAT, only.pos = F, ident.1 = "NK's_Late_disease", ident.2 = "NK's_Prelesion", verbose = FALSE)
eWAT_Cd8posT <- FindMarkers(eWAT, only.pos = F, ident.1 = "Cd8+ Tcells_Late_disease", ident.2 = "Cd8+ Tcells_Prelesion", verbose = FALSE)
eWAT_Cd8posCcl5posTeff <- FindMarkers(eWAT, only.pos = F, ident.1 = "Cd8+ Ccl5+ Teffs_Late_disease", ident.2 = "Cd8+ Ccl5+ Teffs_Prelesion", verbose = FALSE)
eWAT_Cd4Foposxp3posTreg <- FindMarkers(eWAT, only.pos = F, ident.1 = "Cd4+ Foxp3+ Tregs_Late_disease", ident.2 = "Cd4+ Foxp3+ Tregs_Prelesion", verbose = FALSE)
#MemBcells <- FindMarkers(eWAT, only.pos = F, ident.1 = "Memory B cells ?_Late_disease", ident.2 = "Memory B cells ?_Prelesion", verbose = FALSE)
eWAT_Bcells_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "B cells_Late_disease", ident.2 = "B cells_Prelesion", verbose = FALSE)
eWAT_Macrophagesact_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Macrophages activated_Late_disease", ident.2 = "Macrophages activated_Prelesion", verbose = FALSE)
eWAT_Plasmacells_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Plasma cells ?_Late_disease", ident.2 = "Plasma cells ?_Prelesion", verbose = FALSE)
eWAT_Dividing <- FindMarkers(eWAT, only.pos = F, ident.1 = "Dividing cells_Late_disease", ident.2 = "Dividing cells_Prelesion", verbose = FALSE)
eWAT_Fibro_Pi16 <- FindMarkers(eWAT, only.pos = F, ident.1 = "Pi16+ Fibroblasts_Late_disease", ident.2 = "Pi16+ Fibroblasts_Prelesion", verbose = FALSE)
eWAT_Fibroblasts_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Fibroblasts_Late_disease", ident.2 = "Fibroblasts_Prelesion", verbose = FALSE)
#eWAT_Fibro_Mgp_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Mgp+ Fibroblasts_Late_disease", ident.2 = "Mgp+ Fibroblasts_Prelesion", verbose = FALSE)
eWAT_Fibroact_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Fibroblasts activated_Late_disease", ident.2 = "Fibroblasts activated_Prelesion", verbose = FALSE)
eWAT_Macrophages_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Macrophages_Late_disease", ident.2 = "Macrophages_Prelesion", verbose = FALSE)
eWAT_Macrophages_Lgals3_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Lgals3+ Macrophages_Late_disease", ident.2 = "Lgals3+ Macrophages_Prelesion", verbose = FALSE)
eWAT_Int_monocytes <- FindMarkers(eWAT, only.pos = F, ident.1 = "Intermediate monocytes_Late_disease", ident.2 = "Intermediate monocytes_Prelesion", verbose = FALSE)
eWAT_ILC <- FindMarkers(eWAT, only.pos = F, ident.1 = "ILC's_Late_disease", ident.2 = "ILC's_Prelesion", verbose = FALSE)
#NKcellsandorILC2 <- FindMarkers(eWAT, only.pos = F, ident.1 = "NKcells and ILC2's?_Late_disease", ident.2 = "NKcells and ILC2's?_Prelesion", verbose = FALSE)
eWAT_Conv_DC2 <- FindMarkers(eWAT, only.pos = F, ident.1 = "Conv DC2_Late_disease", ident.2 = "Conv DC2_Prelesion", verbose = FALSE)
eWAT_Class_monocytes <- FindMarkers(eWAT, only.pos = F, ident.1 = "Classical and non-classical Monocytes_Late_disease", ident.2 = "Classical and non-classical Monocytes_Prelesion", verbose = FALSE)
eWAT_Granulocytes_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Granulocytes_Late_disease", ident.2 = "Granulocytes_Prelesion", verbose = FALSE)
eWAT_Conv_DC1 <- FindMarkers(eWAT, only.pos = F, ident.1 = "Conv DC1_Late_disease", ident.2 = "Conv DC1_Prelesion", verbose = FALSE)
#Nonclass_monocytes <- FindMarkers(eWAT, only.pos = F, ident.1 = "Non-classical monocytes_Late_disease", ident.2 = "Non-classical monocytes_Prelesion", verbose = FALSE)
eWAT_Endothelialcells_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Endothelial cells_Late_disease", ident.2 = "Endothelial cells_Prelesion", verbose = FALSE)
eWAT_EC_rgs5_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Rgs5+ EC's_Late_disease", ident.2 = "Rgs5+ EC's_Prelesion", verbose = FALSE)
eWAT_VSMCs_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "VSMC's_Late_disease", ident.2 = "VSMC's_Prelesion", verbose = FALSE)
eWAT_MAST_reannot <- FindMarkers(eWAT, only.pos = F, ident.1 = "Hdc+ Cpa3+ MAST cells_Late_disease", ident.2 = "Hdc+ Cpa3+ MAST cells_Prelesion", verbose = FALSE)


sctlist <- c("eWAT_Cd4posT", "eWAT_Cd8posT", "eWAT_NK", "eWAT_Cd8posCcl5posTeff", "eWAT_Cd4Foposxp3posTreg", "eWAT_Bcells_reannot", "eWAT_Macrophagesact_reannot", "eWAT_Plasmacells_reannot", "eWAT_Dividing", "eWAT_Fibro_Pi16", 
             "eWAT_Fibroblasts_reannot", "eWAT_Fibroact_reannot", "eWAT_Macrophages_reannot", "eWAT_Macrophages_Lgals3_reannot", "eWAT_Int_monocytes", "eWAT_ILC", "eWAT_Conv_DC2", "eWAT_Class_monocytes", 
             "eWAT_Granulocytes_reannot", "eWAT_Conv_DC1", "eWAT_Endothelialcells_reannot", "eWAT_EC_rgs5_reannot", "eWAT_VSMCs_reannot", "eWAT_MAST_reannot")


for(i in 1:length(sctlist)) {                
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}

eWAT_Int_monocytes$gene <- rownames(eWAT_Int_monocytes)
eWAT_Int_monocytes_sub <- filter(eWAT_Int_monocytes, !str_detect(gene, "^mt")) %>%
  filter(!str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "Gm42418"))

eWAT_Macrophages_Lgals3_reannot$gene <- rownames(eWAT_Macrophages_Lgals3_reannot)
eWAT_Macrophages_Lgals3_reannot_sub <- filter(eWAT_Macrophages_Lgals3_reannot, !str_detect(gene, "^mt")) %>%
  filter(!str_detect(gene, "^Rp")) %>%
  filter(!str_detect(gene, "^Gm"))

eWAT_EC_rgs5_reannot$gene <- rownames(eWAT_EC_rgs5_reannot)
eWAT_EC_rgs5_reannot_sub <- filter(eWAT_EC_rgs5_reannot, !str_detect(gene, "^Rp"))


topgenes <- c("Plac8", "Ccl4", "Fabp4", "Ly6e", "Mfge8", "Gpnmb", "Trem2", "Cd74", "Apoe", "Lyz2", "Retnla", "H2-Ab1", "H2-Eb1", "Cxcl2")
topgenes <- unique(topgenes)

Idents(eWAT) <- "celltype"
D2 <- DotPlot(eWAT, features = topgenes, idents = cells) + 
  RotatedAxis() +
  ggtitle("eWAT")

V1<-EnhancedVolcano(eWAT_Int_monocytes_sub,
                    lab = rownames(eWAT_Int_monocytes_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in eWAT intermediate monocytes',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V2<-EnhancedVolcano(eWAT_Macrophages_Lgals3_reannot_sub,
                    lab = rownames(eWAT_Macrophages_Lgals3_reannot_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in eWAT Lgals3+ macrophages',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V3<-EnhancedVolcano(eWAT_EC_rgs5_reannot_sub,
                    lab = rownames(eWAT_EC_rgs5_reannot_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in eWAT Rgs5+ endothelial cells',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

grid.arrange(V1,V2,V3, ncol=3)

#### PVAT ####
Idents(PVAT) <- "celltype.group"
DimPlot(PVAT, reduction = "umap", group.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()
p3 <- DimPlot(PVAT, reduction = "umap", split.by = "group_id", group.by = "seurat_clusters", pt.size = 1.2, label = T , repel = T, label.size = 4) + NoLegend() + ggtitle("PVAT")
l3 <- FeaturePlot(PVAT, features = c("Ly6a", "Ly-6A-E-Sca-1")) +
  ggtitle("Ly-6A-E-Sca-1", subtitle = "PVAT")

Idents(PVAT) <- "celltype"
macros_PVAT <- subset(PVAT, idents = c("Macrophages", "Macrophages activated", "Lgals3+ Macrophages"))
Idents(macros_PVAT) <- "group_id"
FeaturePlot(macros_PVAT, split.by = "group_id", features = c("Trem2", "Lgals3"), blend = TRUE, cols = c("navy", "darkgoldenrod1"))
FeaturePlot(macros_PVAT, split.by = "group_id", features = c("Lamp2", "Mrc1"), blend = TRUE, cols = c("navy", "darkgoldenrod1"))
RidgePlot(macros_PVAT, group.by = "group_id", features = c("Trem2", "Lgals3"))


Idents(PVAT) <- "celltype.group"

PVAT_Cd4posT <- FindMarkers(PVAT, only.pos = F, ident.1 = "Cd4+ Tcells_Late_disease", ident.2 = "Cd4+ Tcells_Prelesion", verbose = FALSE)
PVAT_NK <- FindMarkers(PVAT, only.pos = F, ident.1 = "NK's_Late_disease", ident.2 = "NK's_Prelesion", verbose = FALSE)
PVAT_Cd8posT <- FindMarkers(PVAT, only.pos = F, ident.1 = "Cd8+ Tcells_Late_disease", ident.2 = "Cd8+ Tcells_Prelesion", verbose = FALSE)
PVAT_Cd8posCcl5posTeff <- FindMarkers(PVAT, only.pos = F, ident.1 = "Cd8+ Ccl5+ Teffs_Late_disease", ident.2 = "Cd8+ Ccl5+ Teffs_Prelesion", verbose = FALSE)
PVAT_Cd4Foposxp3posTreg <- FindMarkers(PVAT, only.pos = F, ident.1 = "Cd4+ Foxp3+ Tregs_Late_disease", ident.2 = "Cd4+ Foxp3+ Tregs_Prelesion", verbose = FALSE)
#MemBcells <- FindMarkers(PVAT, only.pos = F, ident.1 = "Memory B cells ?_Late_disease", ident.2 = "Memory B cells ?_Prelesion", verbose = FALSE)
PVAT_Bcells_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "B cells_Late_disease", ident.2 = "B cells_Prelesion", verbose = FALSE)
PVAT_Macrophagesact_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Macrophages activated_Late_disease", ident.2 = "Macrophages activated_Prelesion", verbose = FALSE)
PVAT_Plasmacells_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Plasma cells ?_Late_disease", ident.2 = "Plasma cells ?_Prelesion", verbose = FALSE)
PVAT_Dividing <- FindMarkers(PVAT, only.pos = F, ident.1 = "Dividing cells_Late_disease", ident.2 = "Dividing cells_Prelesion", verbose = FALSE)
PVAT_Fibro_Pi16 <- FindMarkers(PVAT, only.pos = F, ident.1 = "Pi16+ Fibroblasts_Late_disease", ident.2 = "Pi16+ Fibroblasts_Prelesion", verbose = FALSE)
PVAT_Fibroblasts_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Fibroblasts_Late_disease", ident.2 = "Fibroblasts_Prelesion", verbose = FALSE)
PVAT_Fibro_Mgp_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Mgp+ Fibroblasts_Late_disease", ident.2 = "Mgp+ Fibroblasts_Prelesion", verbose = FALSE)
PVAT_Fibroact_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Fibroblasts activated_Late_disease", ident.2 = "Fibroblasts activated_Prelesion", verbose = FALSE)
PVAT_Macrophages_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Macrophages_Late_disease", ident.2 = "Macrophages_Prelesion", verbose = FALSE)
PVAT_Macrophages_Lgals3_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Lgals3+ Macrophages_Late_disease", ident.2 = "Lgals3+ Macrophages_Prelesion", verbose = FALSE)
PVAT_Int_monocytes <- FindMarkers(PVAT, only.pos = F, ident.1 = "Intermediate monocytes_Late_disease", ident.2 = "Intermediate monocytes_Prelesion", verbose = FALSE)
PVAT_ILC <- FindMarkers(PVAT, only.pos = F, ident.1 = "ILC's_Late_disease", ident.2 = "ILC's_Prelesion", verbose = FALSE)
#NKcellsandorILC2 <- FindMarkers(PVAT, only.pos = F, ident.1 = "NKcells and ILC2's?_Late_disease", ident.2 = "NKcells and ILC2's?_Prelesion", verbose = FALSE)
PVAT_Conv_DC2 <- FindMarkers(PVAT, only.pos = F, ident.1 = "Conv DC2_Late_disease", ident.2 = "Conv DC2_Prelesion", verbose = FALSE)
PVAT_Class_monocytes <- FindMarkers(PVAT, only.pos = F, ident.1 = "Classical and non-classical Monocytes_Late_disease", ident.2 = "Classical and non-classical Monocytes_Prelesion", verbose = FALSE)
PVAT_Granulocytes_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Granulocytes_Late_disease", ident.2 = "Granulocytes_Prelesion", verbose = FALSE)
PVAT_Conv_DC1 <- FindMarkers(PVAT, only.pos = F, ident.1 = "Conv DC1_Late_disease", ident.2 = "Conv DC1_Prelesion", verbose = FALSE)
#Nonclass_monocytes <- FindMarkers(PVAT, only.pos = F, ident.1 = "Non-classical monocytes_Late_disease", ident.2 = "Non-classical monocytes_Prelesion", verbose = FALSE)
PVAT_Endothelialcells_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Endothelial cells_Late_disease", ident.2 = "Endothelial cells_Prelesion", verbose = FALSE)
PVAT_EC_rgs5_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Rgs5+ Endothelial cells_Late_disease", ident.2 = "Rgs5+ Endothelial cells_Prelesion", verbose = FALSE)
PVAT_VSMCs_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "VSMC's_Late_disease", ident.2 = "VSMC's_Prelesion", verbose = FALSE)
PVAT_MAST_reannot <- FindMarkers(PVAT, only.pos = F, ident.1 = "Hdc+ Cpa3+ MAST cells_Late_disease", ident.2 = "Hdc+ Cpa3+ MAST cells_Prelesion", verbose = FALSE)


sctlist <- c("PVAT_Cd4posT", "PVAT_Cd8posT", "PVAT_NK", "PVAT_Cd8posCcl5posTeff", "PVAT_Cd4Foposxp3posTreg", "PVAT_Bcells_reannot", "PVAT_Macrophagesact_reannot", "PVAT_Plasmacells_reannot", "PVAT_Dividing", "PVAT_Fibro_Pi16", 
             "PVAT_Fibroblasts_reannot", "PVAT_Fibroact_reannot","PVAT_Fibro_Mgp_reannot", "PVAT_Macrophages_reannot", "PVAT_Macrophages_Lgals3_reannot", "PVAT_Int_monocytes", "PVAT_ILC", "PVAT_Conv_DC2", "PVAT_Class_monocytes", 
             "PVAT_Granulocytes_reannot", "PVAT_Conv_DC1", "PVAT_Endothelialcells_reannot", "PVAT_EC_rgs5_reannot", "PVAT_VSMCs_reannot", "PVAT_MAST_reannot")


for(i in 1:length(sctlist)) {                
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}

Idents(PVAT) <- "celltype.group"
PVAT_Tcellresponse <- FindMarkers(PVAT, only.pos = T, ident.1 = "T cells_Late_disease", ident.2 = "T cells_Prelesion", verbose = FALSE)

PVAT_Bcellresponse <- FindMarkers(PVAT, only.pos = T, ident.1 = "B cells_Late_disease", ident.2 = "B cells_Prelesion", verbose = FALSE)

PVAT_NKcresponse <- FindMarkers(PVAT, only.pos = T, ident.1 = "NK cells_Late_disease", ident.2 = "NK cells_Prelesion", verbose = FALSE)

PVAT_Macrophageresponse <- FindMarkers(PVAT, only.pos = T, ident.1 = "Macrophages activated_Late_disease", ident.2 = "Macrophages activated_Prelesion", verbose = FALSE)
PVAT_Macrophageresponse2 <- FindMarkers(PVAT, only.pos = T, ident.1 = "Macrophages_Late_disease", ident.2 = "Macrophages_Prelesion", verbose = FALSE)

PVAT_Monocyteresponse <- FindMarkers(PVAT, only.pos = T, ident.1 = "Monocytes_Late_disease", ident.2 = "Monocytes_Prelesion", verbose = FALSE)

PVAT_Granulocyteresponse <- FindMarkers(PVAT, only.pos = T, ident.1 = "Granulocytes_Late_disease", ident.2 = "Granulocytes_Prelesion", verbose = FALSE)

Idents(PVAT) <- "mouse.fine"
DimPlot(PVAT, reduction = "umap", label = T , repel = T, label.size = 3) + NoLegend()
FeaturePlot(PVAT, features = c("Ptprc", "Tcf21"))

topgenes <- c("Ms4a4c", "Tmsb10", "Ifi213", "Fos", "Ly6a", "Ly6e", "Irf7", "Ifi27l2a", "Cd74", "Retnla", "Lyz2", "H2-Ab1", "H2-Eb1", "Dcn", "S100a10")
topgenes <- unique(topgenes)

Idents(PVAT) <- "celltype"
D3 <- DotPlot(PVAT, features = topgenes, idents = cells) + 
  RotatedAxis() +
  ggtitle("PVAT")


PVAT_Bcells_reannot$gene <- rownames(PVAT_Bcells_reannot)
PVAT_Bcells_reannot_sub <- filter(PVAT_Bcells_reannot, !str_detect(gene, "^mt")) %>%
  filter(!str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "^Gm")) %>%
  filter(!str_detect(gene, "^Rp"))

PVAT_Cd4posT$gene <- rownames(PVAT_Cd4posT)
PVAT_Cd4posT_sub <- filter(PVAT_Cd4posT, !str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "^Rp")) %>%
  filter(!str_detect(gene, "^Gm"))

PVAT_Cd8posT$gene <- rownames(PVAT_Cd8posT)
PVAT_Cd8posT_sub <- filter(PVAT_Cd8posT, !str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "^Rp")) %>%
  filter(!str_detect(gene, "^Gm"))

PVAT_Macrophagesact_reannot$gene <- rownames(PVAT_Macrophagesact_reannot)
PVAT_Macrophagesact_reannot_sub <- filter(PVAT_Macrophagesact_reannot, !str_detect(gene, "Gm8797"))

PVAT_Fibroact_reannot$gene <- rownames(PVAT_Fibroact_reannot)
PVAT_Fibroact_reannot_sub <- filter(PVAT_Fibroact_reannot, !str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "^Gm")) %>%
  filter(!str_detect(gene, "^Rp"))

PVAT_VSMCs_reannot$gene <- rownames(PVAT_VSMCs_reannot)
PVAT_VSMCs_reannot_sub <- filter(PVAT_VSMCs_reannot, !str_detect(gene, "AY036118")) %>%
  filter(!str_detect(gene, "^Gm"))

V1<-EnhancedVolcano(PVAT_Bcells_reannot_sub,
                    lab = rownames(PVAT_Bcells_reannot_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in PVAT B cells',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V2<-EnhancedVolcano(PVAT_Cd4posT_sub,
                    lab = rownames(PVAT_Cd4posT_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in PVAT Cd4+ T cells',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V3<-EnhancedVolcano(PVAT_Cd8posT_sub,
                    lab = rownames(PVAT_Cd8posT_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in PVAT Cd8+ T cells',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V4<-EnhancedVolcano(PVAT_Macrophagesact_reannot_sub,
                    lab = rownames(PVAT_Macrophagesact_reannot_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in PVAT activated macrophages',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V5<-EnhancedVolcano(PVAT_Fibroact_reannot_sub,
                    lab = rownames(PVAT_Fibroact_reannot_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in PVAT activated fibroblasts',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V6<-EnhancedVolcano(PVAT_VSMCs_reannot_sub,
                    lab = rownames(PVAT_VSMCs_reannot_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in PVAT VSMCs',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

grid.arrange(V1,V2,V3,V4,V5,V6, ncol=3)


#### Spleen ####

Idents(Spleen) <- "celltype.group"
DimPlot(Spleen, reduction = "umap", group.by = "group_id", label = T , repel = T, label.size = 3) + NoLegend()
p4 <-DimPlot(Spleen, reduction = "umap", split.by = "group_id", group.by = "seurat_clusters", pt.size = 1.2, label = T , repel = T, label.size = 4) + NoLegend() + ggtitle("Spleen")
l4 <- FeaturePlot(Spleen, features = c("Ly6a", "Ly-6A-E-Sca-1")) + 
  ggtitle("Ly-6A-E-Sca-1", subtitle = "Spleen")


Idents(Spleen) <- "celltype.group"
Spleen_Cd4posT <- FindMarkers(Spleen, only.pos = F, ident.1 = "Cd4+ Tcells_Late_disease", ident.2 = "Cd4+ Tcells_Prelesion", verbose = FALSE)
Spleen_NK <- FindMarkers(Spleen, only.pos = F, ident.1 = "NK's_Late_disease", ident.2 = "NK's_Prelesion", verbose = FALSE)
Spleen_Cd8posT <- FindMarkers(Spleen, only.pos = F, ident.1 = "Cd8+ Tcells_Late_disease", ident.2 = "Cd8+ Tcells_Prelesion", verbose = FALSE)
Spleen_Cd8posCcl5posTeff <- FindMarkers(Spleen, only.pos = F, ident.1 = "Cd8+ Ccl5+ Teffs_Late_disease", ident.2 = "Cd8+ Ccl5+ Teffs_Prelesion", verbose = FALSE)
Spleen_Cd4Foposxp3posTreg <- FindMarkers(Spleen, only.pos = F, ident.1 = "Cd4+ Foxp3+ Tregs_Late_disease", ident.2 = "Cd4+ Foxp3+ Tregs_Prelesion", verbose = FALSE)
#MemBcells <- FindMarkers(Spleen, only.pos = F, ident.1 = "Memory B cells ?_Late_disease", ident.2 = "Memory B cells ?_Prelesion", verbose = FALSE)
Spleen_Bcells_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "B cells_Late_disease", ident.2 = "B cells_Prelesion", verbose = FALSE)
Spleen_Macrophagesact_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Macrophages activated_Late_disease", ident.2 = "Macrophages activated_Prelesion", verbose = FALSE)
Spleen_Plasmacells_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Plasma cells ?_Late_disease", ident.2 = "Plasma cells ?_Prelesion", verbose = FALSE)
Spleen_Dividing <- FindMarkers(Spleen, only.pos = F, ident.1 = "Dividing cells_Late_disease", ident.2 = "Dividing cells_Prelesion", verbose = FALSE)
#Spleen_Fibro_Pi16 <- FindMarkers(Spleen, only.pos = F, ident.1 = "Pi16+ Fibroblasts_Late_disease", ident.2 = "Pi16+ Fibroblasts_Prelesion", verbose = FALSE)
#Spleen_Fibroblasts_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Fibroblasts_Late_disease", ident.2 = "Fibroblasts_Prelesion", verbose = FALSE)
#Spleen_Fibro_Mgp_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Mgp+ Fibroblasts_Late_disease", ident.2 = "Mgp+ Fibroblasts_Prelesion", verbose = FALSE)
#Spleen_Fibroact_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Fibroblasts activated_Late_disease", ident.2 = "Fibroblasts activated_Prelesion", verbose = FALSE)
Spleen_Macrophages_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Macrophages_Late_disease", ident.2 = "Macrophages_Prelesion", verbose = FALSE)
Spleen_Macrophages_Lgals3_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Lgals3+ Macrophages_Late_disease", ident.2 = "Lgals3+ Macrophages_Prelesion", verbose = FALSE)
Spleen_Int_monocytes <- FindMarkers(Spleen, only.pos = F, ident.1 = "Intermediate monocytes_Late_disease", ident.2 = "Intermediate monocytes_Prelesion", verbose = FALSE)
Spleen_ILC <- FindMarkers(Spleen, only.pos = F, ident.1 = "ILC's_Late_disease", ident.2 = "ILC's_Prelesion", verbose = FALSE)
#NKcellsandorILC2 <- FindMarkers(Spleen, only.pos = F, ident.1 = "NKcells and ILC2's?_Late_disease", ident.2 = "NKcells and ILC2's?_Prelesion", verbose = FALSE)
Spleen_Conv_DC2 <- FindMarkers(Spleen, only.pos = F, ident.1 = "Conv DC2_Late_disease", ident.2 = "Conv DC2_Prelesion", verbose = FALSE)
Spleen_Class_monocytes <- FindMarkers(Spleen, only.pos = F, ident.1 = "Classical and non-classical Monocytes_Late_disease", ident.2 = "Classical and non-classical Monocytes_Prelesion", verbose = FALSE)
Spleen_Granulocytes_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Granulocytes_Late_disease", ident.2 = "Granulocytes_Prelesion", verbose = FALSE)
Spleen_Conv_DC1 <- FindMarkers(Spleen, only.pos = F, ident.1 = "Conv DC1_Late_disease", ident.2 = "Conv DC1_Prelesion", verbose = FALSE)
#Nonclass_monocytes <- FindMarkers(Spleen, only.pos = F, ident.1 = "Non-classical monocytes_Late_disease", ident.2 = "Non-classical monocytes_Prelesion", verbose = FALSE)
#Spleen_Endothelialcells_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Endothelial cells_Late_disease", ident.2 = "Endothelial cells_Prelesion", verbose = FALSE)
#Spleen_EC_rgs5_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Rgs5+ Endothelial cells_Late_disease", ident.2 = "Rgs5+ Endothelial cells_Prelesion", verbose = FALSE)
#Spleen_VSMCs_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "VSMC's_Late_disease", ident.2 = "VSMC's_Prelesion", verbose = FALSE)
Spleen_MAST_reannot <- FindMarkers(Spleen, only.pos = F, ident.1 = "Hdc+ Cpa3+ MAST cells_Late_disease", ident.2 = "Hdc+ Cpa3+ MAST cells_Prelesion", verbose = FALSE)


sctlist <- c("Spleen_Cd4posT", "Spleen_Cd8posT", "Spleen_NK", "Spleen_Cd8posCcl5posTeff", "Spleen_Cd4Foposxp3posTreg", "Spleen_Bcells_reannot", "Spleen_Macrophagesact_reannot", "Spleen_Plasmacells_reannot", 
             "Spleen_Dividing", "Spleen_Macrophages_reannot", "Spleen_Macrophages_Lgals3_reannot", "Spleen_Int_monocytes", "Spleen_ILC", "Spleen_Conv_DC2", "Spleen_Class_monocytes", 
             "Spleen_Granulocytes_reannot", "Spleen_Conv_DC1", "Spleen_MAST_reannot")

for(i in 1:length(sctlist)) {            #This writes the result files into separate csv-files in your working directory       
  write.csv2(get(sctlist[i]),
             paste("/scratch/project_2005050/Rstats/",
                   sctlist[i],
                   ".csv"),
             row.names = TRUE)
}

Spleen_Granulocytes_reannot$gene <- rownames(Spleen_Granulocytes_reannot)
Spleen_Granulocytes_reannot_sub <- filter(Spleen_Granulocytes_reannot, !str_detect(gene, "^Rp")) %>%
  filter(!str_detect(gene, "^Gm")) %>%
  filter(!str_detect(gene, "^Ig"))

V1<-EnhancedVolcano(Spleen_Granulocytes_reannot_sub,
                    lab = rownames(Spleen_Granulocytes_reannot_sub),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    title = 'Late disease vs Prelesion in Spleen granulocytes',
                    pCutoff = 5e-2,
                    FCcutoff = 0.4,
                    pointSize = 3.0,
                    col=c('grey', 'slateblue2', 'cyan2', 'cyan4'),
                    labSize = 6.0,
                    titleLabSize = 15.0)

V1

topgenes <- c("Cd74", "Tmsb10", "Lcp1", "Ftl1")

Idents(Spleen) <- "celltype"
D4 <- DotPlot(Spleen, features = topgenes, idents = cells) + 
  RotatedAxis() +
  ggtitle("Spleen")



#Plot tissues together

grid.arrange(D1, D2, D3, D4, ncol=2)

p1 + p2 | p3 + p4

l1
l2
l3
l4

#### B cells ####

#Memory B cells
FeaturePlot(TERVA2_harmony, features = c("Ptprc", "Cd80", "Nt5e", "Cd38", "Cd84", "Cd86", "Pax5", "Spib"))

#B1 cells
FeaturePlot(TERVA2_harmony, features = c("Sdc1", "Cd5"))

#Follicular B cells
FeaturePlot(TERVA2_harmony, features = c("Ptprc", "Ighm", "Cd38", "Cr2", "Cd22", "Cd19", "Pax5"))

#Marginal zone B cells
FeaturePlot(TERVA2_harmony, features = c("Ptprc", "Ighm", "Ighd", "R3", "Cd9", "Cd22", "Cr1", "Pax5", "Ebf1", "Tcf3", "Slc22a2"))

Idents(TERVA2_harmony) <- "celltype"
Bcells <- subset(TERVA2_harmony, idents="B cells")

Bcells <- RunPCA(Bcells, verbose = FALSE)
ElbowPlot(Bcells)

Bcells <- RunUMAP(Bcells, dims = 1:7) %>%
  FindNeighbors(dims = 1:7)



Bcells <- FindClusters(Bcells, resolution = 0)
Bcells <- FindClusters(Bcells, resolution = 0.1) #This seems to be optimal resolution
Bcells <- FindClusters(Bcells, resolution = 0.2) 
Bcells <- FindClusters(Bcells, resolution = 0.3)
Bcells <- FindClusters(Bcells, resolution = 0.4)
Bcells <- FindClusters(Bcells, resolution = 0.5) 
Bcells <- FindClusters(Bcells, resolution = 0.6) 
Bcells <- FindClusters(Bcells, resolution = 0.7)
Bcells <- FindClusters(Bcells, resolution = 0.8)
Bcells <- FindClusters(Bcells, resolution = 0.9)
Bcells <- FindClusters(Bcells, resolution = 1)
b_cls1 <- clustree(Bcells)
b_cls1


Bcells <- SetIdent(Bcells, value = "RNA_snn_res.0.2")
B_cell_subtypes1 <- DimPlot(Bcells, group.by = "Sample")
B_cell_subtypes2 <- DimPlot(Bcells, group.by = "tissue_id")
B_cell_subtypes3 <- DimPlot(Bcells, group.by = "RNA_snn_res.0.2")
B_cell_subtypes1 + B_cell_subtypes2 + B_cell_subtypes3

#Mem B cells
DefaultAssay(Bcells) <- "ADT"
badt_mem <- FeaturePlot(Bcells, reduction = "umap", features = c("CD45", "CD80", "CD73", "CD38","CD86"), cols = c("lightgrey", "darkgreen"))
DefaultAssay(Bcells) <- "RNA"
brna_mem <- FeaturePlot(Bcells, features = c("Ptprc", "Cd80", "Nt5e", "Cd38", "Cd84", "Cd86", "Pax5", "Spib"))
badt_mem
brna_mem

#plasma cells
DefaultAssay(Bcells) <- "ADT"
badt_pl <- FeaturePlot(Bcells, reduction = "umap", features = c("CD45", "CD27", "Ly-6A-E-Sca-1", "CD138-Syndecan-1"), cols = c("lightgrey", "darkgreen"))
DefaultAssay(Bcells) <- "RNA"
brna_pl <- FeaturePlot(Bcells, features = c("Ptprc", "Cd27", "Ly6a", "Sdc1"))
badt_pl
brna_pl

#follicular B cells
DefaultAssay(Bcells) <- "ADT"
badt_f <- FeaturePlot(Bcells, reduction = "umap", features = c("CD45", "IgD", "CD185-CXCR5", "CD23"), cols = c("lightgrey", "darkgreen"))
DefaultAssay(Bcells) <- "RNA"
brna_f <- FeaturePlot(Bcells, features = c("Ptprc", "Ighd", "Cxcr5", "Cd23"))
badt_f
brna_f

Bcell_subtype_markers <- FindAllMarkers(Bcells, logfc.threshold = 0.15)
Top10_Bcellmarkersbycl <- Bcell_subtype_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(Bcells, features = Top10_Bcellmarkersbycl$gene) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

write.csv(Bcell_subtype_markers, "/scratch/project_2005050/Rstats/Bcell_subtype_markers.csv")
write.csv(Top10_Bcellmarkersbycl, "/scratch/project_2005050/Rstats/Top10_Bcellmarkersbycl.csv")


####


Bcells <- NormalizeData(Bcells, assay = "ADT", normalization.method = "CLR")
Bcells <- ScaleData(Bcells, assay = "ADT")
Bcells <- RunPCA(Bcells, assay = "ADT", reduction.name = "apca", features = rownames(Bcells@assays$ADT))
Bcells_harmony <- RunUMAP(Bcells_harmony, assay = "ADT", reduction.name = "apca", features = rownames(Bcells@assays$ADT))


Bcells_harmony <- Bcells %>% RunHarmony("Sample", plot_convergence = T)

Bcells_harmony <- Bcells_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:10, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:10)

Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0)
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.1) #This seems to be optimal resolution
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.2) 
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.3)
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.4)
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.5) 
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.6) 
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.7)
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.8)
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 0.9)
Bcells_harmony <- FindClusters(Bcells_harmony, resolution = 1)
b_cls2 <- clustree(Bcells_harmony)
b_cls2


Bcells_harmony <- RunUMAP(Bcells_harmony, assay = "ADT", reduction.name = "apca", features = rownames(Bcells@assays$ADT))


DefaultAssay(Bcells_harmony) <- "ADT"
badt <- FeaturePlot(Bcells_harmony, reduction = "umap", split.by = "tissue_id", features = c("CD40", "CD86", "CD27", "CD5"), cols = c("lightgrey", "darkgreen"))
DefaultAssay(Bcells_harmony) <- "RNA"
brna <- FeaturePlot(Bcells_harmony, split.by = "tissue_id", features = c("Cd40", "Cd86", "Cd27", "Cd5"))
badt
brna


Bcells_harmony <- SetIdent(Bcells_harmony, value = "RNA_snn_res.0.1")
B_cell_subtypes1_harm <- DimPlot(Bcells_harmony, group.by = "Sample")
B_cell_subtypes2_harm <- DimPlot(Bcells_harmony, group.by = "tissue_id")
B_cell_subtypes3_harm <- DimPlot(Bcells_harmony)
B_cell_subtypes1_harm + B_cell_subtypes2_harm + B_cell_subtypes3_harm

Bcell_subtype_markers_harmony <- FindAllMarkers(Bcells_harmony, logfc.threshold = 0.15)
Top10_Bcellmarkersbycl_harmony <- Bcell_subtype_markers_harmony %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(Bcells_harmony, features = Top10_Bcellmarkersbycl_harmony$gene) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

write.csv(Bcell_subtype_markers_harmony, "/scratch/project_2005050/Rstats/Bcell_subtype_markers_harmony.csv")
write.csv(Top10_Bcellmarkersbycl_harmony, "/scratch/project_2005050/Rstats/Top10_Bcellmarkersbycl_harmony.csv")



### Save whole data as H5 object ####

SaveH5Seurat(object = TERVA2_harmony, overwrite = T, verbose = T)


#### TCR - BCR analysis ####
#There were quite large differences in the successfulness of the TCR/BCR libraries between samples. This should be accounted for when drawing any inferences from the data.

#LD eWAT: two errors for TCR/BCR: Low fraction of reads mapped to any VDJ gene (18.80 % in TCR, 28.45 % in BCR, ideal should be > 50 %), and low fraction of valid barcodes in TCR (74.83 %, ideal > 85 %)
#PL eWAT: low fraction valid barcodes in TCR (exactly 85.00 %), and low fraction of reads mapped to any vdj gene (48.32 %, ideal > 50 %)


#BiocManager::install("scRepertoire")
#devtools::install_github("rnabioco/djvdj")

library(djvdj)
library(scRepertoire)

Idents(TERVA2_harmony) <- "celltype"
TERVA2_harmony_small <- subset(TERVA2_harmony, idents = c("29","30","31", "32", "33",
                                                          "34", "35"), invert = TRUE)

TERVA2_harmony_small <- subset(TERVA2_harmony, idents = c("29_Late_disease", "29_Prelesion","30_Late_disease", "30_Prelesion","31_Late_disease", 
                                                          "31_Prelesion", "32_Late_disease", "32_Prelesion", "33_Late_disease",
                                                          "34_Late_disease", "35_Late_disease"), invert = TRUE)

clonotype_id <- TERVA2_harmony_small@meta.data$TCR_clonotype_id
TERVA2_harmony_small@meta.data$clonotype_id <- clonotype_id

chains <- TERVA2_harmony_small@meta.data$TCR_chains
TERVA2_harmony_small@meta.data$chains <- chains


TERVA2_harmony_small <- calc_diversity(
  input       = TERVA2_harmony_small,                       # Seurat object
  cluster_col = "celltype.group",
  clonotype_col = "clonotype_id"
)

TCRusage_simpson <- plot_diversity(
  input       = TERVA2_harmony_small,                       # Seurat object
  cluster_col = "celltype.group",
  clonotype_col = "clonotype_id"  
) +
  ggtitle("TCR diversity (Simpson)") +
  theme(axis.text.x=element_text(angle=45, hjust=1))



TERVA2_harmony_small <- calc_abundance(
  input = TERVA2_harmony_small,
  cluster_col = "celltype.group",
  clonotype_col = "TCR_clonotype_id",
  prefix = "TCR_"
)



TCRsimilarity <- calc_similarity(TERVA2_harmony_small,
                                 cluster_col = "Sample",
                                 method = abdiv::jaccard,
                                 clonotype_col = "clonotype_id",
                                 prefix = "TCR_",
                                 return_mat = TRUE)

plot_similarity(
  TERVA2_harmony_small,
  cluster_col = "Sample"
  ) 



chains <- TERVA2_harmony_small@meta.data$BCR_chains
TERVA2_harmony_small@meta.data$chains <- chains

calc_gene_usage(
  input       = TERVA2_harmony_small,
  gene_cols   = c("BCR_v_gene","BCR_j_gene"),                   # Column containing genes
  cluster_col = "Sample",                    # Column containing cell clusters to compare
  chain       = "IGH",                      # Chain to calculate gene usage for
  chain_col   = "chains"                    # Column containing chains
)


plot_gene_usage(input       = TERVA2_harmony_small,
                gene_cols   = c("BCR_v_gene","BCR_j_gene"),                   # Column(s) containing genes
                type        = "bar",
                cluster_col = "celltype.group",                   # Type of plot, 'heatmap' or 'bar'
                chain       = "IGK",                      # Chain to plot
                n_genes     = 20,                         # The number of top genes to plot
                plot_colors = "cyan4"
)


ggs <- plot_gene_usage(
  input       = TERVA2_harmony_small,
  gene_cols   = c("BCR_v_gene","BCR_j_gene"),      # Column(s) containing genes
  cluster_col = "Sample",                    # Column containing cell clusters to compare
  chain       = "IGH",                      # Chain to plot
  plot_colors = "cyan4",
  n_genes     = 20
)

plot_grid(plotlist = ggs)



