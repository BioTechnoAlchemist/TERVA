#### DESeq2 ####

library(Matrix.utils)
library(edgeR)
library(magrittr)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(dplyr)

counts <- TERVA2_harmony@assays$RNA@counts 

metadata <- TERVA2_harmony@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(TERVA2_harmony@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

#make sample_id and group_id, which are later expected
sample_id <- as.factor(TERVA2_harmony@meta.data$Sample)
sce$sample_id <- sample_id

coldata <- as.data.frame(sce@colData)

coldata <- mutate(coldata,
                  Status = case_when(
                    startsWith(Sample, "LD") ~ "Late_disease",
                    startsWith(Sample, "PL") ~ "Prelesion"
                  ))

#group id to case/control

group_id <- coldata$Status
sce$group_id <- as.factor(group_id)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids


# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))


# Total number of samples 
ns <- length(sids)

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reorder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

# Explore the different components of list
str(pb)

pb$"NA" = NULL

# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)



gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 


metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id) 

metadata$cluster_id <- factor(metadata$cluster_id)

head(metadata, n = 10)

# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
clusters

# Subset the metadata to only the B cells
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[3]), ]
head(cluster_metadata)

# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
cluster_metadata

# Subset the counts to only the B cells
counts <- pb[[clusters[3]]]

cluster_counts <- as.data.frame(as.matrix(counts[, which(colnames(counts) %in% rownames(cluster_metadata))]))


# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))      

# Create DESeq2 object        
dds <- DESeqDataSetFromMatrix(round(cluster_counts), 
                              colData = cluster_metadata, 
                              design = ~ group_id)

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
DESeq2::plotPCA(rld, intgroup = "group_id")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
plotDispEsts(dds)

# Check the coefficients for the comparison
resultsNames(dds)
#[1] "Intercept"                          "group_id_Prelesion_vs_Late_disease"

# Generate results object
res <- results(dds, 
               name = "group_id_Prelesion_vs_Late_disease",
               alpha = 0.05)

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res <- lfcShrink(dds, 
                 coef = "group_id_Prelesion_vs_Late_disease",
                 res=res,
                 type = "apeglm")

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl 

# Function to run DESeq2 and get results for all clusters
## x is index of cluster in clusters vector on which to run function
## A is the sample group to compare
## B is the sample group to compare against (base level)

get_dds_resultsAvsB <- function(x, A, B){
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  counts <- pb[[clusters[x]]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  
  #all(rownames(cluster_metadata) == colnames(cluster_counts))        
  
  dds <- DESeqDataSetFromMatrix(round(cluster_counts), 
                                colData = cluster_metadata, 
                                design = ~ 1)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
  
  # Plot PCA
  
  DESeq2::plotPCA(rld, intgroup = "group_id")
  ggsave(paste0("/scratch/project_2005050/Rstats", clusters[x], "_specific_PCAplot.png"))
  
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Plot heatmap
  png(paste0("/scratch/project_2005050/Rstats", clusters[x], "_specific_heatmap.png"))
  pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
  dev.off()
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  # Plot dispersion estimates
  png(paste0("/scratch/project_2005050/Rstats", clusters[x], "_dispersion_plot.png"))
  plotDispEsts(dds)
  dev.off()
  
  # Output results of Wald test for contrast for A vs B
  contrast <- c("group_id", levels(cluster_metadata$group_id)[A], levels(cluster_metadata$group_id)[B])
  
  resultsNames(dds)
  res <- results(dds, 
                 contrast = contrast,
                 alpha = 0.05)
  
  res <- lfcShrink(dds, 
                   contrast =  contrast,
                   res=res,
                   type = "ashr")
  # Set thresholds
  padj_cutoff <- 0.05
  
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            paste0("/scratch/project_2005050/Rstats", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            paste0("/scratch/project_2005050/Rstats", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  ## ggplot of top genes
  normalized_counts <- counts(dds, 
                              normalized = TRUE)
  
  ## Order results by padj values
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n=20)
  
  
  top20_sig_norm <- data.frame(normalized_counts) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% top20_sig_genes)
  
  gathered_top20_sig <- top20_sig_norm %>%
    gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
  
  gathered_top20_sig <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))
  
  ## plot using ggplot2
  ggplot(gathered_top20_sig) +
    geom_point(aes(x = gene, 
                   y = normalized_counts, 
                   color = group_id), 
               position=position_jitter(w=0.1,h=0)) +
    scale_y_log10() +
    xlab("Genes") +
    ylab("log10 Normalized Counts") +
    ggtitle("Top 20 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("/scratch/project_2005050/Rstats", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_top20_DE_genes.png"))
  
}

map(1:length(clusters), get_dds_resultsAvsB, A = 2, B = 1)