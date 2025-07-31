setwd("~/Desktop/TERVA/TERVA2DATA/DATA and result from CSC server/GSEA/eWAT/")
library(enrichplot)
library(tidyverse)
library(multienrichjam)
library(dplyr)

# Function to capitalize the first letter of each word in a comma-separated string
capitalize_words <- function(words) {
  word_list <- unlist(strsplit(words, ","))
  modified_string <- str_to_sentence(word_list, locale = "en")
  modified_string <- paste(modified_string, collapse = ", ")
  return(modified_string)
}

data_up <- read.csv("upregulated results/eWAT_Fibro_Pi16_gProfiler_mmusculus_26-10-2023_14-00-42__intersections.csv", sep = ";") #This originates from a web interface run of g:profiler at https://biit.cs.ut.ee/gprofiler/gost.
data_up <- data_up %>%
  mutate(term_name_mod = str_to_sentence(data_up$term_name), gene_name = sapply(intersections, capitalize_words)) # change gene names to mouse format and capitalize the first letter of the term sentences.
data_down <- read.csv("downregulated results/eWAT_Pi16Fibro_gProfiler_mmusculus_16-04-2024_21-32-54__intersections.csv", sep = ";")
data_down <- data_down %>%
  mutate(term_name_mod = str_to_sentence(data_down$term_name), gene_name = sapply(intersections, capitalize_words))

gp_mod_up <- data_up
gp_mod_up$GeneRatio = paste0(gp_mod_up$intersection_size,  "/", gp_mod_up$query_size)
gp_mod_up$BgRatio = paste0(gp_mod_up$term_size, "/", gp_mod_up$effective_domain_size)

gp_mod_down <- data_down
gp_mod_down$GeneRatio = paste0(gp_mod_down$intersection_size,  "/", gp_mod_down$query_size)
gp_mod_down$BgRatio = paste0(gp_mod_down$term_size, "/", gp_mod_down$effective_domain_size)

# Visualization of the GSEA results ####

#Only upregulated enrichResults to be used in cnetplot etc.

gp_mod_enrich_up <- enrichDF2enrichResult(enrichDF = gp_mod_up, keyColname = "term_id", geneColname = "gene_name", 
                                          pvalueColname = "p_value", descriptionColname = "term_name_mod", pvalueCutoff = 0.05)

up1 <- barplot(gp_mod_enrich_up, showCategory = 40, x = "GeneRatio", font.size = 14, title = "Up-regulated terms and pathways") +
  scale_fill_viridis_c(direction = -1) +
  theme(text = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(vjust = 1)) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=70))

up2 <- cnetplot(gp_mod_enrich_up, showCategory = 12, layout = "circle", colorEdge = TRUE, color.params = list(edge = "#005C5C")) + 
  theme(text = element_text(size = 16))


#Only downregulated enrichResults to be used in cnetplot etc.

gp_mod_enrich_down <- enrichDF2enrichResult(enrichDF = gp_mod_down, keyColname = "term_id", geneColname = "gene_name", 
                                            pvalueColname = "p_value", descriptionColname = "term_name_mod", pvalueCutoff = 0.05)

down1 <- barplot(gp_mod_enrich_down, showCategory = 40, x = "GeneRatio", font.size = 14, title = "Down-regulated terms and pathways") +
  scale_fill_viridis_c(direction = -1) +
  theme(text = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(vjust = 1)) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=70))

down2 <- cnetplot(gp_mod_enrich_down, showCategory = 12, layout = "circle", colorEdge = TRUE, color.params = list(edge = "#005C5C")) + 
  theme(text = element_text(size = 16))


eWAT <- (up1 + down1) /
  (up2 + down2)

# Text analysis of the results ####

#Upregulated terms ####
terms <- data_up$term_name
source <- data_up$source



datasets <- list()
for (source_name in unique(source)) {
  subset_terms <- terms[source == source_name]
  term_list <- unlist(strsplit(subset_terms, "\n"))
  # Filter terms based on the keywords
  immune_keywords <- c("immune", "immunity", "inflammation", "inflammatory", "mhc", "complement", "c1", "c2", "c3", "c4", "hypersensitivity", 
                       "autoimmune", "leukocyte", "defense", "cytokine", "interferon", "interleukin", "antigen", "neutrophil", "granulocyte", "platelet", "macrophage",
                       "th1", "th2", "th17", "dendritic cell", "t cell", "b cell", "rheuma", "arthritis", "atherosclerosis", "lupus", "diabetes", "viral", "bacterial", "graft")
  #Please note that including "t cell" in the list, some none immune related phrases might appear if the CONTAIN it anywhere: eg. "fa(t cell) differentiation". Then again using " t cell" end up dropping the immune related terms which START with t cell, e.g. "T cell activation". Therefore be careful with this. 
  immune_related_terms <- term_list[grep(paste(immune_keywords, collapse = "|"), tolower(term_list))]
  immune_related_terms <- sort(immune_related_terms)
  datasets[[source_name]] <- immune_related_terms
}

for (source_name in names(datasets)) {
  cat("Data Source:", source_name, "\n")
  print(datasets[[source_name]])
  cat("\n")
}

#Downregulated terms ####
terms <- data_down$term_name
source <- data_down$source

datasets <- list()
for (source_name in unique(source)) {
  subset_terms <- terms[source == source_name]
  term_list <- unlist(strsplit(subset_terms, "\n"))
  immune_keywords <- c("immune", "immunity", "inflammation", "inflammatory", "mhc", "complement", "c1", "c2", "c3", "c4", "hypersensitivity", 
                       "autoimmune", "leukocyte", "defense", "cytokine", "interferon", "interleukin", "antigen", "neutrophil", "granulocyte", "platelet", "macrophage",
                       "th1", "th2", "th17", "dendritic cell", "t cell", "b cell", "rheuma", "arthritis", "atherosclerosis", "lupus", "diabetes", "viral", "bacterial", "graft")
  immune_related_terms <- term_list[grep(paste(immune_keywords, collapse = "|"), tolower(term_list))]
  immune_related_terms <- sort(immune_related_terms)
    datasets[[source_name]] <- immune_related_terms
}

for (source_name in names(datasets)) {
  cat("Data Source:", source_name, "\n")
  print(datasets[[source_name]])
  cat("\n")
}

