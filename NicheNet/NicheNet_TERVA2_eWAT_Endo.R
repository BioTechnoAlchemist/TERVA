
setwd("/scratch/project_2005050/Rstats/NicheNet/")

library(nichenetr) 
library(tidyverse)
library(circlize)
library(Seurat)

.libPaths(c("/projappl/project_2005050/project_rpackages_4.1.1/", .libPaths()))
library(colorBlindness)

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "tissue_id") 
eWAT <- subset(TERVA2_harmony, idents = "eWAT")

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
eWAT <- SetIdent(eWAT, value = "celltype")
receiver = "Gpihbp1+ Fabp4+ Endothelial cells"  #Just change this to your cell type of interest
expressed_genes_receiver = get_expressed_genes(receiver, eWAT, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Cd248+ Pi16+ Fibroblasts", "Ccl11+ Fibroblasts", "Mgp+ Aebp1+ Activated fibroblasts", "Folr2+ Lyve1+ M2 Macrophages", "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages", "Pf4+ Retnla+ Macrophages", "Trem2+ Lgals3+ Macrophages",
                     "Innate lymphoid cells", "Conventional Dendritic cells DC1", "Mesothelial cells", "Rgs5+ Endothelial cells", "Vascular smooth muscle cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, eWAT, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# Gene set of interest for the receiver population

seurat_obj_receiver= subset(eWAT, idents = receiver)
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

DotPlot(subset(eWAT, idents = sender_celltypes), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
DotPlot(subset(eWAT, idents = sender_celltypes), features = best_upstream_ligands %>% rev(), cols = "RdYlBu", split.by = "group_id") + RotatedAxis()

# Active target gene inference

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

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

# DE analysis for each sender cell type (This will take some time if many senders)
DE_table_all = Idents(eWAT) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = eWAT, condition_colname = "group_id", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]
colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

DotPlot(eWAT, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

# Circos plots ####

# Calculate average ligand expr in sender cells
avg_expression_ligands = AverageExpression(eWAT, features = best_upstream_ligands, assays = "RNA")

# Assign ligand to senders
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

# Determine which ligands are expressed by specific cells
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

#Check the senders from the beginning of the workflow

Cd248Pi16Fibro_specific_ligands = sender_ligand_assignment$`"Cd248+ Pi16+ Fibroblasts"` %>% names() %>% setdiff(general_ligands)
Ccl11Fibro_specific_ligands = sender_ligand_assignment$`"Ccl11+ Fibroblasts"` %>% names() %>% setdiff(general_ligands)
MgpAebp1Actfibro_specific_ligands = sender_ligand_assignment$`"Mgp+ Aebp1+ Activated fibroblasts"` %>% names() %>% setdiff(general_ligands)
Folr2Lyve1Macro_specific_ligands = sender_ligand_assignment$`Folr2+ Lyve1+ M2 Macrophages` %>% names() %>% setdiff(general_ligands)
Infl_Macro_specific_ligands = sender_ligand_assignment$`Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages` %>% names() %>% setdiff(general_ligands)
Pf4RetnlaMacro_specific_ligands = sender_ligand_assignment$`Pf4+ Retnla+ Macrophages` %>% names() %>% setdiff(general_ligands)
Trem2Lgals3Macro_specific_ligands = sender_ligand_assignment$`Trem2+ Lgals3+ Macrophages` %>% names() %>% setdiff(general_ligands)
ILC_specific_ligands = sender_ligand_assignment$`Innate lymphoid cells` %>% names() %>% setdiff(general_ligands)
ConvDC1_specific_ligands = sender_ligand_assignment$`Conventional Dendritic cells DC1` %>% names() %>% setdiff(general_ligands)
Mesothelial_specific_ligands = sender_ligand_assignment$`Mesothelial cells` %>% names() %>% setdiff(general_ligands)
Rgs5EC_specific_ligands = sender_ligand_assignment$`Rgs5+ Endothelial cells` %>% names() %>% setdiff(general_ligands)
VSMC_specific_ligands = sender_ligand_assignment$`Vascular smooth muscle cells` %>% names() %>% setdiff(general_ligands)


ligand_type_indication_df = tibble(
  ligand_type = c(rep("Cd248Pi16Fibro_specific_ligands", times = Cd248Pi16Fibro_specific_ligands %>% length()),
                  rep("Ccl11Fibro_specific_ligands", times = Ccl11Fibro_specific_ligands %>% length()),
                  rep("MgpAebp1Actfibro_specific_ligands", times = MgpAebp1Actfibro_specific_ligands %>% length()),
                  rep("Folr2Lyve1Macro_specific_ligands", times = Folr2Lyve1Macro_specific_ligands %>% length()),
                  rep("Infl_Macro_specific_ligands", times = Infl_Macro_specific_ligands %>% length()),
                  rep("Pf4RetnlaMacro_specific_ligands", times = Pf4RetnlaMacro_specific_ligands %>% length()),
                  rep("Trem2Lgals3Macro_specific_ligands", times = Trem2Lgals3Macro_specific_ligands %>% length()),
                  rep("ILC_specific_ligands", times = ILC_specific_ligands %>% length()),
                  rep("ConvDC1_specific_ligands", times = ConvDC1_specific_ligands %>% length()),
                  rep("Mesothelial_specific_ligands", times = Mesothelial_specific_ligands %>% length()),
                  rep("Rgs5EC_specific_ligands", times = Rgs5EC_specific_ligands %>% length()),
                  rep("VSMC_specific_ligands", times = VSMC_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(Cd248Pi16Fibro_specific_ligands, Ccl11Fibro_specific_ligands, MgpAebp1Actfibro_specific_ligands, Folr2Lyve1Macro_specific_ligands, Infl_Macro_specific_ligands, Pf4RetnlaMacro_specific_ligands, Trem2Lgals3Macro_specific_ligands, ILC_specific_ligands, ConvDC1_specific_ligands, Mesothelial_specific_ligands, Rgs5EC_specific_ligands, VSMC_specific_ligands, general_ligands))

# Define ligand-target links of interest

active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "Late_disease-DE") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)


displayAllColors(mycolors, color = "white")

grid_col_ligand =c("General" = mycolors[[1]],
                   "MgpAebp1Actfibro_specific_ligands" = mycolors[[2]],
                   "Folr2Lyve1Macro_specific_ligands" = mycolors[[3]],
                   "Infl_Macro_specific_ligands" = mycolors[[4]],
                   "Pf4RetnlaMacro_specific_ligands" = mycolors[[5]],
                   "Trem2Lgals3Macro_specific_ligands" = mycolors[[6]],
                   "ILC_specific_ligands" = mycolors[[7]],
                   "ConvDC1_specific_ligands" = mycolors[[8]],
                   "Mesothelial_specific_ligands" = mycolors[[9]],
                   "Rgs5EC_specific_ligands" = mycolors[[10]],
                   "VSMC_specific_ligands" = mycolors[[11]],
                   "Cd248Pi16Fibro_specific_ligands" = mycolors[[12]],
                   "Ccl11Fibro_specific_ligands" = mycolors[[13]]
)

grid_col_target =c(
  "Late_disease-DE" = "goldenrod1")
grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% dplyr::select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
# Prepare circos visualization

target_order = circos_links$target %>% unique() 
#Check that the ligand order is correct, if not copy from above where the specific ligand types are defined
ligand_order = c(Cd248Pi16Fibro_specific_ligands, Ccl11Fibro_specific_ligands, MgpAebp1Actfibro_specific_ligands, Folr2Lyve1Macro_specific_ligands, Infl_Macro_specific_ligands, Pf4RetnlaMacro_specific_ligands, Trem2Lgals3Macro_specific_ligands, ILC_specific_ligands, ConvDC1_specific_ligands, Mesothelial_specific_ligands, Rgs5EC_specific_ligands, VSMC_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

# Define gaps between segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

unique(circos_links$ligand_type) #check which specific links exist
gaps = c(rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
         width_ligand_target)

# Render a circos plot. This version uses differing transparency levels determined by the regulatory potential value of each ligand-target interaction

circos.par(gap.degree = 2)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)

circos.clear()

