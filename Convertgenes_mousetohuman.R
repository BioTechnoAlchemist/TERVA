setwd("/scratch/project_2005050/Rstats")

.libPaths(c("/projappl/project_2005050/project_rpackages_4.1.1", .libPaths()))

library(SeuratDisk)

LDAO <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDAOafterQC2.h5Seurat") #ReQC
PLAO <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLAOafterQC2.h5Seurat") #ReQC
LDeWAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/LD_eWATafterQC.h5Seurat")
PLeWAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLeWATafterQC2.h5Seurat") #ReQC
LDPVAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDPVATafterQC.h5Seurat")
PLPVAT <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLPVATafterQC2.h5Seurat") #ReQC
LDSpleen <- LoadH5Seurat("/scratch/project_2005050/Rstats/LDSpleenafterQC.h5Seurat")
PLSpleen <- LoadH5Seurat("/scratch/project_2005050/Rstats/PLSpleenafterQC.h5Seurat")


# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

allgenes <- as.data.frame(rownames(LDAO)) #set data
TERVA2_genes_converted <- convertMouseGeneList(allgenes)
TERVA2_genes_converted <-as.data.frame(TERVA2_genes_converted)

write.table(allgenes, "allgenes.txt", row.names=F, col.names = F, quote=F, sep="\t")
write.table(TERVA2_genes_converted, "humanorthologs.txt", row.names=F, col.names = F, quote=F, sep="\t")


