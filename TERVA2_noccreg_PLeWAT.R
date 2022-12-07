#### TERVA2 scRNA-seq experiment 2021 ####
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
sc <- load10X("/Users/limikk/Desktop/TERVA/TERVA2DATA/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/count", includeFeatures = "Gene Expression")
sc

#Channel with 32285 genes and 7292 cells


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

DR = sc$metaData[,sc$DR]
DR = as.data.frame(DR)
genes<- as.vector(c("Igkc","Ighg1","Ighm","Iglc2","Iglc3","Iglc1","Ighe","Ighg3","Ighd")) #With this gene set the global cont. fraction is 6.45 %
#allgenes <- as.vector(sc[["toc"]]@Dimnames[[1]])
#plotMarkerMap(sc, geneSet = genes)
#plotMarkerMap(sc, "Igkc")
#plotMarkerMap(sc, "Ighm")
#plotMarkerMap(sc, "Ighg1")
#plotMarkerMap(sc, "Iglc2")
#plotMarkerMap(sc, "Iglc3")
#plotMarkerMap(sc, "Iglc1")
#plotMarkerMap(sc, "Ighe")
#plotMarkerMap(sc, "Ighg3")
#plotMarkerMap(sc, "Ighd")


sc <- autoEstCont(sc)

#1346 genes passed tf-idf cut-off and 470 soup quantile filter.  Taking the top 100.
#Using 439 independent estimates of rho.
#Estimated global rho of 0.02
#Channel with 32285 genes and 7292 cells


#head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#plotMarkerDistribution(sc)

#useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = genes)) #Use if AutoEstCont is not stringent enough (< 2%)
#plotMarkerMap(sc, geneSet = genes, useToEst = useToEst)
#sc <- calculateContaminationFraction(sc, list(IG = genes), useToEst = useToEst) #Estimated global contamination fragtion of 1.46%
#head(sc$metaData)

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
tenx <- Read10X("/Users/limikk/Desktop/TERVA/TERVA2DATA/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/count/filtered_feature_bc_matrix/")
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
srat <- FindClusters(srat, resolution = 0.5, verbose = F)
DimPlot(srat, reduction = "umap", pt.size=1, group.by = "RNA_snn_res.0.5", label= T) 
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

nExp_poi <- round(0.045*nrow(srat@meta.data)) #0.045 is based on the 10X estimate of multiplet formation rate dependent on cell loading
srat <- doubletFinder_v3(srat, PCs = 1:20, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
#srat_reuse <- doubletFinder_v3(srat, PCs = 1:20, pN = 0.25, pK = 0.22, nExp = nExp_poi, reuse.pANN = "pANN_0.25_0.22_120", sct = TRUE)

DF.name = colnames(srat@meta.data)[grepl("DF.classification", colnames(srat@meta.data))]
VlnPlot(srat, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

plot1 <- DimPlot(srat, group.by="DF.classifications_0.25_0.06_303", reduction="umap", pt.size=1, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB"))
plot2 <- DimPlot(srat, reduction = "umap", pt.size=1, group.by = "RNA_snn_res.0.5")

plot1 + plot2 
#plot1 + plot1_1

srat_DF <- subset(srat, subset = DF.classifications_0.25_0.06_303 == "Singlet") #Creates a new Seurat Object from the singlets.

plot3 <- DimPlot(srat_DF, group.by="DF.classifications_0.25_0.02_303", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB"))
plot4 <- DimPlot(srat_DF, reduction = "umap", pt.size=0.5, label = T, group.by = "RNA_snn_res.0.5")

plot3 + plot4 
plot4 # looks same as before so ok



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
srat_N <- FindClusters(srat_N, resolution = 0.5, verbose = F)
plot7 <- DimPlot(srat_N, reduction = "umap", pt.size=0.5, label = T, group.by = "RNA_snn_res.0.5")
FeaturePlot(srat_N, features = c("Myh11", "Cnn1", "Tcf21", "C1qa", "Lyz2", "H2-Ab1", "Vwf", "Flt1", "Kdr", "Ltb", "Pdgfra", "Lum", "Pecam1", "Cdh5", "Cd79a", "Cd3d", "Cnp"), pt.size = 0.2,
            ncol = 4)
FeaturePlot(srat_N, features = c("Il7r", "Icos", "Kit", "Rag1"), pt.size = 0.2,
            ncol = 3) #Il7r = CD127, Icos = CD278, Kit = CD117
VlnPlot(srat_N, features = c("Il7r", "Icos", "Kit"), pt.size = 0.2, ncol = 3)
DotPlot(srat_N, features = c("Il7r", "Icos", "Kit")) + RotatedAxis()

VlnPlot(srat_N, features = c("Acta2", "Cd3d", "Cd79a", "Spp1", "Ccl7", "Il1b", "Ccl22", "Jchain", "Dcn"), pt.size = 0.2, ncol = 3)
plot9 <- DotPlot(srat_N, features = c("Acta2", "Cd3d", "Cd79a", "Spp1", "Ccl7", "Il1b", "Ccl22", "Jchain", "Dcn")) + RotatedAxis()
plot10 <- DoHeatmap(subset(srat_N, downsample = 100), features = c("Acta2", "Cd3d", "Cd79a", "Spp1", "Ccl7", "Il1b", "Ccl22", "Jchain", "Dcn"), size = 3)

plot7 + plot9 + plot10

#### Add clonotypes ####

#devtools::install_github("rnabioco/djvdj")
library(djvdj)

srat_TCR <- import_vdj(input = srat_N, vdj_dir = "~/Desktop/TERVA/TERVA2DATA/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/vdj_t/", filter_paired = F)
srat_BCR <- import_vdj(input = srat_N, vdj_dir = "~/Desktop/TERVA/TERVA2DATA/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/vdj_b/", filter_paired = F)

srat_VDJ <- import_vdj(input = srat_N, vdj_dir = "~/Desktop/TERVA/TERVA2DATA/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/vdj_t/", prefix = "TCR_", filter_paired = F)
srat_VDJ <- import_vdj(input = srat_VDJ, vdj_dir = "~/Desktop/TERVA/TERVA2DATA/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/vdj_b/", prefix = "BCR_", filter_paired = F)


#### VDJ analysis ####

table(!is.na(srat_VDJ$TCR_clonotype_id), !is.na(srat_VDJ$BCR_clonotype_id))


#BiocManager::install("dittoSeq")
library(dittoSeq)

v1 <- dittoDimPlot(srat_VDJ, var = "TCR_chains", size = 2,
                   cells.use = grepl("TRA", srat_VDJ$TCR_chains, fixed = TRUE), 
                   show.others = TRUE)

#unique(srat_VDJ@meta.data$BCR_isotype)
# NA      "IGHM"  "IGHG"  "None"  "IGHD"  "Multi" "IGHA" 


v2 <- dittoDimPlot(srat_VDJ, var = "BCR_isotype", size = 2,
                   cells.use = grepl("IGHM", srat_VDJ$BCR_isotype, fixed = TRUE), 
                   show.others = TRUE)

v3 <- FeaturePlot(srat_VDJ, features = c("Myh11", "Cnn1", "C1qa", "Lyz2", "Flt1", "Pdgfra", "Lum", "Pecam1", "Cdh5", "Cd79a", "Cd3d"), pt.size = 1,
                  ncol = 4)

plot7 + v1 + v2

#Let's look at what cellranger estimated as T/B cells
library(rjson)
Tcells <- as.data.frame(fromJSON(file = "~/Desktop/TERVA/TERVA2DATA/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/vdj_t/cell_barcodes.json"))
names(Tcells)[1] <- "Tcells" 
Bcells <- as.data.frame(fromJSON(file = "~/Desktop/TERVA/TERVA2DATA/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/outs/per_sample_outs/PL_eWAT_GEX_Ab_VDJ_noIntr_SC5PPE/vdj_b/cell_barcodes.json"))
names(Bcells)[1] <- "Bcells" 

Tcells<- as.vector(Tcells)
Bcells<- as.vector(Bcells)
Tc <- DimPlot(object = srat_VDJ, cells.highlight = Tcells, cols.highlight = "red", cols = "gray", order = TRUE, group.by = "RNA_snn_res.0.7", label = T)
Bc <- DimPlot(object = srat_VDJ, cells.highlight = Bcells, cols.highlight = "blue", cols = "gray", order = TRUE, group.by = "RNA_snn_res.0.7", label = T)
Tc + Bc

#### SAVE H5-object to use in subsequent analyses ####

library(SeuratDisk)
library(SeuratData)

SaveH5Seurat(object = srat_VDJ, overwrite = T, verbose = T)

h5<-LoadH5Seurat("PLeWATafterQC.h5Seurat") #works, srat_VDJ and h5 look similar (small changes in the order of items though)

remove(list = c("srat_BCR","srat_TCR","h5","sub","out"))




#### DE ####
arrange(srat18markers,
        desc(avg_log2FC)) %>%
  rownames()

srat0markers <- FindMarkers(srat_N, ident.1 = 0, min.pct = 0.1) #B cells, top markers Cd79a (B cell marker), Ly6d (May act as a specification marker at earliest stage specification of lymphocytes between B- and T-cell development. Marks the earliest stage of B-cell specification.), H2-DMb2 (MHC), Fcmr (Fc fragment of IgM receptor), Cd74 (B cell marker)
top0 <- head(srat0markers, 10)

tcf21 <- FeaturePlot(srat_N, features = "Tcf21")
ptprc <- FeaturePlot(srat_N, features = "Ptprc")

tcf21 + ptprc


srat0markers <- FindMarkers(srat_N, ident.1 = 0, min.pct = 0.1) # "C1qa"           "C1qc"           "Pf4"            "C1qb"           "F13a1"          "Folr2"          "Cbr2"           "Mrc1"           "Retnla"         "Apoe"
srat1markers <- FindMarkers(srat_N, ident.1 = 1, min.pct = 0.1) # "Dcn"            "Gpx3"           "Gsn"            "Penk"           "Sult1e1"        "Col3a1"         "Dpep1"          "Clec3b"         "Serping1"       "Ccl11"
srat2markers <- FindMarkers(srat_N, ident.1 = 2, min.pct = 0.1) # "Aqp1"          "Gpihbp1"       "Cldn5"         "Fabp4"         "Mgll"          "Rflnb"         "Kdr"           "Cdh5"          "Cav1"          "Tspan13" 
srat3markers <- FindMarkers(srat_N, ident.1 = 3, min.pct = 0.1) # "Lyz1"           "Il1b"           "Tnip3"          "Bcl2a1d"        "Cd52"           "Ear2"           "Bcl2a1b"        "Fcrls"          "Areg"           "Lsp1" 
srat4markers <- FindMarkers(srat_N, ident.1 = 4, min.pct = 0.1) # "Rgs5"          "Kcnj8"         "Ndufa4l2"      "Abcc9"         "Mfge8"         "Cox4i2"        "Pdgfrb"        "Art3"          "Myl9"          "Steap4"
srat5markers <- FindMarkers(srat_N, ident.1 = 5, min.pct = 0.1) # "Fbn1"          "Cd55"          "Mfap5"         "Pi16"          "Tnfaip6"       "Has1"          "Ugdh"          "Fn1"           "Cd248"         "Fstl1" 
srat6markers <- FindMarkers(srat_N, ident.1 = 6, min.pct = 0.1) # "Rbp7"          "Tm4sf1"        "Flt1"          "Ptprb"         "Btnl9"         "Sema3g"        "Cdh5"          "Efnb2"         "Hey1"          "Igfbp3"
srat7markers <- FindMarkers(srat_N, ident.1 = 7, min.pct = 0.1) # "Plac8"         "Gm26917"       "Lgals3"        "Il1b"          "Cybb"          "Gngt2"         "Lst1"          "Cd300c2"       "Tyrobp"        "Cebpb"
srat8markers <- FindMarkers(srat_N, ident.1 = 8, min.pct = 0.1) # "Ccl5"          "Gzma"          "Nkg7"          "Ms4a4b"        "Ly6c2"         "Vps37b"        "Cd3g"          "Il2rb"         "Cd3d"          "Cd3e" 
srat9markers <- FindMarkers(srat_N, ident.1 = 9, min.pct = 0.1) # "Plbd1"          "Cst3"           "Naaa"           "Wdfy4"          "H2afz"          "Irf8"           "Ckb"            "Gm2a"           "Naga"           "Hist1h2ap"
srat10markers <- FindMarkers(srat_N, ident.1 = 10, min.pct = 0.1) # "Tagln"         "Fabp4"         "Rbp7"          "Gng11"         "Cfd"           "Car3"          "Rflnb"         "Tcf15"         "Gpihbp1"       "Sncg"
srat11markers <- FindMarkers(srat_N, ident.1 = 11, min.pct = 0.1) # "Sfrp4"         "Dpt"           "Cilp"          "Gpx3"          "S100a6"        "Timp2"         "Gsn"           "Hmcn2"         "Igfbp6"        "Pi16"
srat12markers <- FindMarkers(srat_N, ident.1 = 12, min.pct = 0.1) # "Pf4"           "Mafb"          "Lgmn"          "Retnla"        "Wfdc17"        "Fxyd2"         "Tgfbi"         "F13a1"         "Clec10a"       "Tfrc"
srat13markers <- FindMarkers(srat_N, ident.1 = 13, min.pct = 0.1) # "Slpi"           "Upk3b"          "Clu"            "Aebp1"          "Nkain4"         "Igfbp6"         "Fmod"           "Igfbp5"         "C2"             "Msln"
srat14markers <- FindMarkers(srat_N, ident.1 = 14, min.pct = 0.1) # "Csf2"          "Hilpda"        "Itk"           "Ccr8"          "Tnfaip3"       "Rora"          "Vps37b"        "Rgs2"          "Areg"          "Tnfrsf18" 
srat15markers <- FindMarkers(srat_N, ident.1 = 15, min.pct = 0.1) # "Cd79a"          "Ly6d"           "Cd79b"          "Ms4a1"          "Cd37"           "Fcmr"           "H2-DMb2"        "Iglc2"          "Iglc3"          "Ccr7" 
srat16markers <- FindMarkers(srat_N, ident.1 = 16, min.pct = 0.1) # "Igkv14-126"     "Igkv1-117"      "Igkv4-55"       "Jchain"         "Igkc"           "Igkv3-2"        "Ighm"           "Ighv11-2"       "Ighv1-80"       "Ighv9-3"
srat17markers <- FindMarkers(srat_N, ident.1 = 17, min.pct = 0.1) # "S100a9"        "S100a8"        "Retnlg"        "Il1b"          "Msrb1"         "Il1r2"         "Csf3r"         "Pglyrp1"       "Slpi"          "Il1rn" 
srat18markers <- FindMarkers(srat_N, ident.1 = 18, min.pct = 0.1) # "Acta2"         "Tagln"         "Myh11"         "Tpm2"          "Myl9"          "Mustn1"        "Tpm1"          "Mylk"          "Des"           "Gm13889"

adip_genes <- c("Adipoq", "Pdgfra", "Ces1f", "Btc", "Apoe", "Cacna1a", "Prune2", "Mt2", "Tcf21")
sctadipmarkers <- FindAllMarkers(srat_N, return.thres= 0.1, min.pct = 0.25, features=intersect(rownames(srat_N), adip_genes)) # The same as on previous row but this targets only your genes of interest.

angio_genes <- c("Cdh5", "Sdpr", "Egfl7", "Ptprb", "Ecscr", "Cldn5", "Icam2", "Slc9a3r2", "Myh11", "Cnn1")
sctangiomarkers <- FindAllMarkers(srat_N, return.thres= 0.1, min.pct = 0.25, features=intersect(rownames(srat_N), angio_genes)) # The same as on previous row but this targets only your genes of interest.



####

#Subset the suspicious cluster 1

srat_NO1 <- subset(srat_N, subset = seurat_clusters %in% c(0,2:9,12,14,16)) #Creates a new Seurat Object from the singlets.
srat_NO1 <- NormalizeData(srat_NO1, assay = "RNA", verbose = F)
srat_NO1 <- FindVariableFeatures(srat_NO1, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(srat_NO1), 10)
plot11 <- VariableFeaturePlot(srat_NO1)
plot12 <- LabelPoints(plot = plot11, points = top10, repel = TRUE)
plot11 + plot12

srat_NO1 <- ScaleData(srat_NO1, assay = "RNA", verbose = T, features = rownames(srat_NO1))
srat_NO1 <- RunPCA(srat_NO1, verbose = F)
ElbowPlot(srat_NO1)

srat_NO1 <- RunUMAP(srat_NO1, dims = 1:18, verbose = F)
srat_NO1 <- FindNeighbors(srat_NO1, dims = 1:18, verbose = F)
srat_NO1 <- FindClusters(srat_NO1, resolution = 0.5, verbose = F)

plot13 <- DimPlot(srat_NO1, reduction = "umap", pt.size=0.5, label = T, group.by = "RNA_snn_res.0.5")
FeaturePlot(srat_NO1, features = c("Acta2", "Cd3d", "Cd79a", "Spp1", "Ccl7", "Il1b", "Ccl22", "Jchain", "Dcn"), pt.size = 0.2,
            ncol = 3)

plot14 <- VlnPlot(srat_NO1, features = c("Cd3d", "Cd79a"), pt.size = 0.2, ncol = 2)
plot15 <- DotPlot(srat_NO1, features = c("Cd3d", "Cd79a")) + RotatedAxis()
plot16 <- DoHeatmap(subset(srat_NO1, downsample = 100), features = c("Cd3d", "Cd79a", "Acta2", "Spp1", "Jchain", "Dcn"), size = 3)

plot13 + plot14 + plot15 + plot16

#### SAVE H5-object to use in subsequent analyses ####

library(SeuratDisk)
library(SeuratData)

SaveH5Seurat(object = srat_VDJ, overwrite = T, verbose = T)

h5<-LoadH5Seurat("PLeWATafterQC.h5Seurat") #works, srat_VDJ and h5 look similar (small changes in the order of items though)

remove(list = c("srat_BCR","srat_TCR","h5","sub","out"))
