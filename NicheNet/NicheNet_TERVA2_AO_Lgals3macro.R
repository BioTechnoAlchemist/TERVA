setwd("/scratch/project_2005050/Rstats/NicheNet")

.libPaths(c("/projappl/project_2005050/", .libPaths()))

#devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(circlize)

TERVA2_harmony <- SetIdent(TERVA2_harmony,value = "tissue_id") 
AO <- subset(TERVA2_harmony, idents = "Aorta")


ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

#Rename the human genes to mouse orthologs
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

#Define sender and receiver cell populations

## receiver
AO <- SetIdent(AO, value = "celltype")
receiver = "Lgals3+ Macrophages" #Just change this to your cell type of interest
expressed_genes_receiver = get_expressed_genes(receiver, AO, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Cd8+ Tcells", "Cd4+ Tcells", "Cd8+ Ccl5+ Teffs", "Cd4+ Foxp3+ Tregs", "Intermediate Monocytes", "Classical and non-classical Monocytes", "Innate lymphoid cells", "B cells", "Macrophages", "Macrophages activated", "Vascular smooth muscle cells", "Conventional Dendritic cells 2")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, AO, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# Gene set of interest for the receiver population

seurat_obj_receiver= subset(AO, idents = receiver)
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

DotPlot(subset(AO, idents = sender_celltypes), features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
DotPlot(subset(AO, idents = sender_celltypes), features = best_upstream_ligands %>% rev(), cols = "RdYlBu", split.by = "group_id") + RotatedAxis()

# Active target gene inference

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
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
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
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
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
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
DE_table_all = Idents(AO) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = AO, condition_colname = "group_id", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]
colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

# Circos plots ####

# Calculate average ligand expr in sender cells
avg_expression_ligands = AverageExpression(AO, features = best_upstream_ligands, assays = "RNA")

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

Cd8T_specific_ligands = sender_ligand_assignment$`Cd8+ Tcells` %>% names() %>% setdiff(general_ligands)
Cd4T_specific_ligands = sender_ligand_assignment$`Cd4+ Tcells` %>% names() %>% setdiff(general_ligands)
Teffs_specific_ligands = sender_ligand_assignment$`Cd8+ Ccl5+ Teffs` %>% names() %>% setdiff(general_ligands)
Tregs_specific_ligands = sender_ligand_assignment$`Cd4+ Foxp3+ Tregs` %>% names() %>% setdiff(general_ligands)
IntMono_specific_ligands = sender_ligand_assignment$`Intermediate Monocytes` %>% names() %>% setdiff(general_ligands)
Mono_specific_ligands = sender_ligand_assignment$`Classical and non-classical Monocytes` %>% names() %>% setdiff(general_ligands)
ILC_specific_ligands = sender_ligand_assignment$`Innate lymphoid cells` %>% names() %>% setdiff(general_ligands)
VSMC_specific_ligands = sender_ligand_assignment$`Vascular smooth muscle cells` %>% names() %>% setdiff(general_ligands)
Macro_specific_ligands = sender_ligand_assignment$`Macrophages` %>% names() %>% setdiff(general_ligands)
Macroact_specific_ligands = sender_ligand_assignment$`Macrophages activated` %>% names() %>% setdiff(general_ligands)
Bcells_specific_ligands = sender_ligand_assignment$`B cells` %>% names() %>% setdiff(general_ligands)
DC2_specific_ligands = sender_ligand_assignment$`Conventional Dendritic cells 2` %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("Cd8T-specific", times = Cd8T_specific_ligands %>% length()),
                  rep("Cd4T-specific", times = Cd4T_specific_ligands %>% length()),
                  rep("Teffs-specific", times = Teffs_specific_ligands %>% length()),
                  rep("Tregs-specific", times = Tregs_specific_ligands %>% length()),
                  rep("IntMono-specific", times = IntMono_specific_ligands %>% length()),
                  rep("Mono-specific", times = Mono_specific_ligands %>% length()),
                  rep("ILC-specific", times = ILC_specific_ligands %>% length()),
                  rep("VSMC-specific", times = VSMC_specific_ligands %>% length()),
                  rep("Macro-specific", times = Macro_specific_ligands %>% length()),
                  rep("Macroact-specific", times = Macroact_specific_ligands %>% length()),
                  rep("Bcells-specific", times = Bcells_specific_ligands %>% length()),
                  rep("DC2-specific", times = DC2_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(Cd8T_specific_ligands, Cd4T_specific_ligands, Teffs_specific_ligands, Tregs_specific_ligands, IntMono_specific_ligands, Mono_specific_ligands, ILC_specific_ligands, VSMC_specific_ligands, Macro_specific_ligands, Macroact_specific_ligands, Bcells_specific_ligands, DC2_specific_ligands, general_ligands))

# Define ligand-target links of interest

active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "Late_disease-DE") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = mycols[[1]],
                   "Cd8T-specific" = mycols[[11]],
                   "Cd4T-specific" = mycols[[3]],
                   "Teffs-specific" = mycols[[4]],
                   "Tregs-specific" = mycols[[5]],
                   "IntMono-specific" = mycols[[6]],
                   "Mono-specific" = mycols[[2]],
                   "ILC-specific" = mycols[[8]],
                   "VSMC-specific" = mycols[[12]],
                   "Macro-specific" = mycols[[10]],
                   "Macroact-specific" = mycols[[7]],
                   "Bcells-specific" = mycols[[9]],
                   "DC2-specific" = mycols[[14]]
)
grid_col_target =c(
  "Late_disease-DE" = mycols[[13]])


grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

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
ligand_order = c(Cd8T_specific_ligands, Cd4T_specific_ligands, Teffs_specific_ligands, Tregs_specific_ligands, IntMono_specific_ligands, Mono_specific_ligands, ILC_specific_ligands, VSMC_specific_ligands, Macro_specific_ligands, Macroact_specific_ligands, Bcells_specific_ligands, DC2_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

# Define gaps between segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

unique(circos_links$ligand_type) #check which specific links exist
gaps = c(rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target)

# Render the circos plot (the first option will draw the links with same transparency)

circos.par(gap.degree = 5)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

circos.clear() #This empties the plot

# This version uses differing transparency levels determined by the regulatory potential value of each ligand-target interaction

circos.par(gap.degree = 3)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
