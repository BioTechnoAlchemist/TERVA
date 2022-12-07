# Clustree

setwd("~/Desktop/TERVA/TERVA2DATA/")

#install.packages('clustree')
library(clustree)

#First run FindClusters on the seurat object with different resolutions and the check below that there are at least two columns with the prefix RNA_snn_res.
head(srat_N[[]])

#Run clustree

plot_cls1 <- clustree(srat_N)
plot_cls2 <- clustree(srat, node_colour= "Cd79a", node_colour_aggr= "median")
plot_cls3 <- clustree(srat, node_colour= "Cd3d", node_colour_aggr= "median")


plot_cls1 + plot_cls2 + plot_cls3
