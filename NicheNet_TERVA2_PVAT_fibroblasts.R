#This script was created by Aydin Bölük and Lea Mikkola.

setwd("/scratch/project_2005050/Rstats/NicheNet/")

library(nichenetr)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(circlize)


#### PVAT - Fibroblasts ####
Idents(TERVA2_harmony) <- "tissue_id"
pvat <- subset(TERVA2_harmony, tissue_id == "PVAT")

# Read in NicheNet's ligand-target prior model, ligand-receptor network and weighted integrated networks.

lr_network = readRDS("lr_network.rds")
ligand_target_matrix =  readRDS(file = "ligand_target_matrix.rds")
weighted_networks = readRDS("weighted_networks.rds")


lr_network <- lr_network |>
  distinct(from, to, database)

weighted_networks_lr <- weighted_networks$lr_sig |>
  inner_join(lr_network, by = c("from", "to"))

table(pvat$celltype)

pvat <- SetIdent(pvat, value = "celltype")
pvat <- UpdateSeuratObject(pvat)

celltypes_non_obese <- unique(pvat$celltype.group)[grepl("Non-Obese", x = unique(pvat$celltype.group))]

celltypes_obese <- unique(pvat$celltype.group)[!(grepl("Non-Obese", x = unique(pvat$celltype.group)))]

celltypes_obese <- gsub("_Obese", "", x = celltypes_obese)
celltypes_non_obese <- gsub("_Non-Obese", "", x = celltypes_non_obese)

as.data.frame(table(pvat$celltype))
intersect(celltypes_non_obese, celltypes_obese)

sender_cells <- intersect(celltypes_non_obese, celltypes_obese)[c(-27, -12)]

sender_celltypes = c("B cells", "Mfap4+ Fibroblasts", "Ccl11+ Fibroblasts", "Mgp+ Aebp1+ Activated fibroblasts", "Lef1+ Tcf7+ Cd4+ T cells", "Gpihbp1+ Fabp4+ Endothelial cells",
                     "Classical and non-classical monocytes", "Conventional Dendritic cells DC1", "Folr2+ Lyve1+ M2 Macrophages", "Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages",
                     "Pf4+ Retnla+ Macrophages", "Mesothelial cells", "Fscn1+ Apol7c+ Dendritic cells", "Intermediate monocytes", "Natural killer cells", "Plasma cells", "s100a9+/a8+ Granulocytes",
                     "Vascular smooth muscle cells")

#### Run NicheNet ####

nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = pvat,
  receiver = "Cd248+ Pi16+ Fibroblasts",
  condition_colname = "group_id", condition_oi = "Obese", condition_reference = "Non-Obese",
  sender = sender_celltypes,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks,
  assay_oi = "RNA",
  lfc_cutoff = 0.10,
  geneset = "DE",
  top_n_ligands = 30,
  organism = "mouse", 
  expression_pct = 0.10)

#Make visualizations

vis_ligand_aupr <- nichenet_output$ligand_activities %>% filter(test_ligand %in% nichenet_output$top_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity",
                                     color = "darkorange", legend_title = "AUPR") + 
  theme(axis.text.x.top = element_blank(),
        text = element_text(size = 14),
        legend.key.width = unit(1, "cm"),
        legend.position = "right")
  

p_ligand_expression_sender = nichenet_output$ligand_expression_dotplot + 
  ggtitle("Expression in Sender") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(1, "cm"),
        legend.position = "right")

p_lfc_sender = nichenet_output$ligand_differential_expression_heatmap + 
  theme(text = element_text(size = 14),
        axis.text.x.top = element_text(angle = 60),
        legend.key.width = unit(1, "cm"),
        legend.position = "right")

p_regulatorypotential = nichenet_output$ligand_target_heatmap +
  theme(text = element_text(size = 14),
        axis.text.x.top = element_text(angle = 60),
        legend.key.width = unit(1, "cm"),
        legend.position = "top")

p_interactionpotential = nichenet_output$ligand_receptor_heatmap + 
  theme(text = element_text(size = 14),
  axis.text.x.top = element_text(angle = 60),
  legend.key.width = unit(1, "cm"),
  legend.position = "top",
  axis.title.x = element_blank())

(p_ligand_aupr + p_lfc_sender + p_ligand_expression_sender +
  patchwork::plot_layout(widths = c(1, 1.5, 2), ncol = 3))

p_regulatorypotential + p_interactionpotential
 

#DotPlot(pvat, features = nichenet_output$top_ligands %>% rev(), split.by = "group_id", cols = "RdYlBu") + RotatedAxis()


# Circos Plot ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#### Calculate average ligand expression in sender cells ####

avg_expression_ligands <- AggregateExpression(
  pvat,
  assays = "RNA",
  features = nichenet_output$top_ligands,
)

#### Assign ligands to sender cells ####
sender_ligand_assignment <- avg_expression_ligands$RNA |>
  apply(MARGIN = 1, FUN = function(ligand_expression){
    ligand_expression > (ligand_expression |> mean() + ligand_expression |> sd())
  }) |>
  t()

sender_ligand_assignment <- sender_ligand_assignment |>
  apply(MARGIN = 2, FUN = function(x){
    x[x == TRUE]
  }) |>
  keep(function(x){length(x) > 0})


sender_ligand_assignment <- sender_ligand_assignment

all_assigned_ligands <- sender_ligand_assignment |>
  lapply(FUN = function(x){names(x)}) |>
  unlist()

unique_ligands <- all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands <- nichenet_output$top_ligands |> setdiff(unique_ligands)

Bcell_specific_ligands = sender_ligand_assignment$`B cells` |> names() |> setdiff(general_ligands)
Mfap4Fibroblasts_specific_ligands = sender_ligand_assignment$`Mfap4+ Fibroblasts` |> names() |> setdiff(general_ligands)
Ccl11Fibroblasts_specific_ligands = sender_ligand_assignment$`Ccl11+ Fibroblasts` |> names() |> setdiff(general_ligands)
MgpFibroblasts_specific_ligands = sender_ligand_assignment$`Mgp+ Aebp1+ Activated fibroblasts` |> names() |> setdiff(general_ligands)
Cd4Tcells_specific_ligands = sender_ligand_assignment$`Lef1+ Tcf7+ Cd4+ T cells` |> names() |> setdiff(general_ligands)
EC_specific_ligands = sender_ligand_assignment$`Gpihbp1+ Fabp4+ Endothelial cells` |> names() |> setdiff(general_ligands)
Clandnonclmonocytes_specific_ligands = sender_ligand_assignment$`"Classical and non-classical monocytes` |> names() |> setdiff(general_ligands)
ConvDC1_specific_ligands = sender_ligand_assignment$`Conventional Dendritic cells DC1` |> names() |> setdiff(general_ligands)
M2Macrophages_specific_ligands = sender_ligand_assignment$`Folr2+ Lyve1+ M2 Macrophages` |> names() |> setdiff(general_ligands)
InflMacrophages_specific_ligands = sender_ligand_assignment$`Ccl4+ Cxcl2+ Ccl3+ Inflammatory Macrophages` |> names() |> setdiff(general_ligands)
Pf4Macro_specific_ligands = sender_ligand_assignment$`Pf4+ Retnla+ Macrophages` |> names() |> setdiff(general_ligands)
Mesothelialcells_specific_ligands = sender_ligand_assignment$`Mesothelial cells` |> names() |> setdiff(general_ligands)
Fscn1DC_specific_ligands = sender_ligand_assignment$`Fscn1+ Apol7c+ Dendritic cells` |> names() |> setdiff(general_ligands)
Intermonocytes_specific_ligands = sender_ligand_assignment$`Intermediate monocytes` |> names() |> setdiff(general_ligands)
NK_specific_ligands = sender_ligand_assignment$`Natural killer cells` |> names() |> setdiff(general_ligands)
plasmacells_specific_ligands = sender_ligand_assignment$`Plasma cells` |> names() |> setdiff(general_ligands)
Granulocytes_specific_ligands = sender_ligand_assignment$`s100a9+/a8+ Granulocytes` |> names() |> setdiff(general_ligands)
VSMC_specific_ligands = sender_ligand_assignment$`Vascular smooth muscle cells` |> names() |> setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("Bcell_specific_ligands", times = Bcell_specific_ligands |> length()),
                  rep("Mfap4Fibroblasts_specific_ligands", times = Mfap4Fibroblasts_specific_ligands |> length()),
                  rep("Ccl11Fibroblasts_specific_ligands", times = Ccl11Fibroblasts_specific_ligands |> length()),
                  rep("MgpFibroblasts_specific_ligands", times = MgpFibroblasts_specific_ligands |> length()),
                  rep("Cd4Tcells_specific_ligands", times = Cd4Tcells_specific_ligands |> length()),
                  rep("EC_specific_ligands", times = EC_specific_ligands |> length()),
                  rep("Clandnonclmonocytes_specific_ligands", times = Clandnonclmonocytes_specific_ligands |> length()),
                  rep("ConvDC1_specific_ligands", times = ConvDC1_specific_ligands |> length()),
                  rep("M2Macrophages_specific_ligands", times = M2Macrophages_specific_ligands |> length()),
                  rep("InflMacrophages_specific_ligands", times = InflMacrophages_specific_ligands |> length()),
                  rep("Pf4Macro_specific_ligands", times = Pf4Macro_specific_ligands |> length()),
                  rep("Mesothelialcells_specific_ligands", times = Mesothelialcells_specific_ligands |> length()),
                  rep("Fscn1DC_specific_ligands", times = Fscn1DC_specific_ligands |> length()),
                  rep("Intermonocytes_specific_ligands", times = Intermonocytes_specific_ligands |> length()),
                  rep("NK_specific_ligands", times = NK_specific_ligands |> length()),
                  rep("plasmacells_specific_ligands", times = plasmacells_specific_ligands |> length()),
                  rep("Granulocytes_specific_ligands", times = Granulocytes_specific_ligands |> length()),
                  rep("VSMC_specific_ligands", times = VSMC_specific_ligands |> length()),
                  rep("General", times = general_ligands |> length())),
  ligand = c(Bcell_specific_ligands, 
             Mfap4Fibroblasts_specific_ligands, 
             Ccl11Fibroblasts_specific_ligands, 
             MgpFibroblasts_specific_ligands, 
             Cd4Tcells_specific_ligands, 
             EC_specific_ligands, 
             Clandnonclmonocytes_specific_ligands, 
             ConvDC1_specific_ligands, 
             M2Macrophages_specific_ligands, 
             InflMacrophages_specific_ligands, 
             Pf4Macro_specific_ligands, 
             Mesothelialcells_specific_ligands, 
             Fscn1DC_specific_ligands, 
             Intermonocytes_specific_ligands, 
             NK_specific_ligands, 
             plasmacells_specific_ligands, 
             Granulocytes_specific_ligands, 
             VSMC_specific_ligands, 
             general_ligands))

ligand_type_indication_df

## Only one endothelial cell specific ligand

#### Define the ligand-target links of interest ####

active_ligand_target_links_df <- nichenet_output$ligand_target_df |>
  mutate(target_type = "Obese-DE") |>
  inner_join(ligand_type_indication_df)

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand <- c(
  "General" = "#b73779",
  "EC_specific_ligands" = "#5ec962"
)

grid_col_target =c(
  "Obese-DE" = "#51127c")

grid_col_tbl_ligand <- tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target <- tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links <- circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links <- circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle <- circos_links %>% select(ligand,target, weight)

ligand_color <- circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color <- ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color <- circos_links %>% distinct(target,color_target_type)
grid_target_color <- target_color$color_target_type %>% set_names(target_color$target)

grid_col <- c(grid_ligand_color,grid_target_color)

transparency <- circos_links %>%
  mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>%
  mutate(transparency = 1-weight) %>% .$transparency



#### Prepare the circos visualization: order ligands and targets ####

target_order <- circos_links$target |> unique()
ligand_order = c(general_ligands, EC_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

#### Prepare the circos visualization: define the gaps between the different segments ####

width_same_cell_same_ligand_type <- 0.5
width_different_cell <- 6
width_ligand_target <- 15
width_same_cell_same_target_type <- 1


gaps <- c(rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "EC_specific_ligands") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Obese-DE") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)


#### Render the circos plot (all links same transparancy) ####
# Only the widths of the blocks that indicate each target gene is proportional the ligand-target regulatory potential (\~prior knowledge supporting the regulatory interaction).

circos.par(gap.degree = 1)
chordDiagram(
  links_circle,
  directional = 1,
  order=order,
  link.sort = TRUE,
  link.decreasing = FALSE,
  grid.col = grid_col,
  transparency = transparency,
  diffHeight = 0.005,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.075))

# we go back to the first track and customize sector labels
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter,
                CELL_META$ylim[1],
                CELL_META$sector.index,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.55),
                cex = 0.7)
  },
  bg.border = NA)

circos.clear()
