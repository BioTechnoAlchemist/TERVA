#Subset the suspicious cluster 1

srat_NO1 <- subset(TERVA2_harmony, subset = seurat_clusters == "10") 
srat_NO1 <- NormalizeData(srat_NO1, assay = "RNA", verbose = F)
srat_NO1 <- FindVariableFeatures(srat_NO1, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(srat_NO1), 10)
VariableFeaturePlot(srat_NO1)
LabelPoints(plot = plot11, points = top10, repel = TRUE)


srat_NO1 <- ScaleData(srat_NO1, assay = "RNA", verbose = T, features = rownames(srat_NO1))
srat_NO1 <- RunPCA(srat_NO1, verbose = F)
ElbowPlot(srat_NO1)

srat_NO1 <- RunUMAP(srat_NO1, dims = 1:20, verbose = F)
srat_NO1 <- FindNeighbors(srat_NO1, dims = 1:20, verbose = F)
srat_NO1 <- FindClusters(srat_NO1, resolution = 0.5, verbose = F)

DimPlot(srat_NO1, reduction = "umap", pt.size=0.5, label = T, group.by = "RNA_snn_res.0.5")
FeaturePlot(srat_NO1, features = c("Nkg7", "Cd3d", "Cd79a"), pt.size = 0.2,
            ncol = 3)

VlnPlot(srat_NO1, features = c("Cd3d", "Cd79a", "Nkg7"), pt.size = 0.2, ncol = 2)
DotPlot(srat_NO1, features = c("Cd3d", "Cd79a", "Nkg7")) + RotatedAxis()
DoHeatmap(subset(srat_NO1, downsample = 100), features = c("Cd3d", "Cd79a", "Nkg7"), size = 3)

VlnPlot(srat_NO1, features = "nFeature_RNA")
