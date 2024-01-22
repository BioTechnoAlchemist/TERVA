
setwd("/scratch/project_2005050/Rstats/NicheNet/")

library(nichenetr) 
library(tidyverse)
library(circlize)
library(Seurat)

.libPaths(c("/projappl/project_2005050/project_rpackages_4.1.1/", .libPaths()))
library(colorBlindness)

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "tissue_id") 
Spleen <- subset(TERVA2_harmony, idents = "Spleen")

ligand_target_matrix =  readRDS(file = "ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


#Rename the human genes to mouse orthologs
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

#Define sender and receiver cell populations


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

## receiver
Spleen <- SetIdent(Spleen, value = "celltype")
receiver = "s100a9+/a8+ Granulocytes"  #Just change this to your cell type of interest
expressed_genes_receiver = get_expressed_genes(receiver, Spleen, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("B cells", "Lef1+ Tcf7+ Cd4+ T cells", "Lef1+ Tcf7+ Cd8a+ Cd8b1+ T cells",  "Cd8+ Ccl5+ Nkg7+ Cytotoxic T cells","Foxp3+ Regulatory T cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Spleen, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# Gene set of interest for the receiver population

seurat_obj_receiver= subset(Spleen, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group_id"]])

condition_oi = "Late_disease"
condition_reference = "Prelesion" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10, logfc.threshold = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.10) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# Define possible ligands expressed by the senders
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# Ligand activity analysis

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

DotPlot(subset(Spleen, idents = sender_celltypes), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
DotPlot(subset(Spleen, idents = sender_celltypes), features = best_upstream_ligands %>% rev(), cols = "RdYlBu", split.by = "group_id") + RotatedAxis()

# Active target gene inference
#DOESN'T WORK, INPUT NULL

# Receptors of top target ligands

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

# Receptors of top target ligands (bona fide evidence)

lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

# Define log fold change for ligands in senders
#DOESN'T WORK, INPUT NULL

# Circos plots ####

# DOESN'T WORK AS ONLY ONE UPSTREAM LIGAND AVAILABLE DOWNSTREAM