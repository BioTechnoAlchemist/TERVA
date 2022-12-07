# Proportion of cell types in tissues by state
setwd("/scratch/project_2005050/Rstats")

Idents(TERVA2_harmony) <- "celltype"
count_table <- as.data.frame(table(TERVA2_harmony@meta.data$celltype, TERVA2_harmony@meta.data$Sample, TERVA2_harmony@meta.data$group_id))

library(tidyverse)

#### All in one go ####
rm <- c("29":"35")
celltable <- count_table %>%
  group_by(Var2) %>%
#  filter(!Var1 %in% rm)  %>%
  mutate(sum_cells = sum(Freq),
    Prop_cells = Freq/sum_cells*100,
    sum_prop = sum(Prop_cells))
  

celltableplot <- lapply(split(celltable, celltable$Var2), 
                        function(x) 
                          ggplot(x, aes(x=Var1, y=Prop_cells)) +
                          geom_bar(stat="identity", fill= "blue") +
                          ylim(0, 30) +
                          theme(axis.text.x=element_text(angle=45, hjust=1)) +
                          ggtitle(x$Var2[1]))
  
#### Sample by sample ####

celltable$status <- celltable$Var3

p1 <- celltable %>%
  filter(Var2 %in% c("LDAO", "PLAO"))
p1 <- ggplot(p1, aes(x=Var1, y=Prop_cells, fill= status)) +
  geom_bar(stat="identity", position = "dodge") +
  ylim(0, 30) +
  scale_fill_manual(values = c("darkgoldenrod3", "dodgerblue3")) +
  labs(title = "Proportion of cell types in late disease and prelesion aorta", x = "Cell types"
       , y = "Proportion, max 30 %") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 14),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))


p2 <- celltable %>%
  filter(Var2 %in% c("LDeWAT", "PLeWAT"))
p2 <- ggplot(p2, aes(x=Var1, y=Prop_cells, fill= status)) +
  geom_bar(stat="identity", position = "dodge") +
  ylim(0, 30) +
  scale_fill_manual(values = c("darkgoldenrod3", "dodgerblue3")) +
  labs(title = "Proportion of cell types in late disease and prelesion eWAT", x = "Cell types"
       , y = "Proportion, max 30 %") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 14),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  NoLegend()

p3 <- celltable %>%
  filter(Var2 %in% c("LDPVAT", "PLPVAT"))
p3 <- ggplot(p3, aes(x=Var1, y=Prop_cells, fill= status)) +
  geom_bar(stat="identity", position = "dodge") +
  ylim(0, 30) +
  scale_fill_manual(values = c("darkgoldenrod3", "dodgerblue3")) +
  labs(title = "Proportion of cell types in late disease and prelesion PVAT", x = "Cell types"
       , y = "Proportion, max 30 %") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 14),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  NoLegend()

p4 <- celltable %>%
  filter(Var2 %in% c("LDSpleen", "PLSpleen"))
p4 <- ggplot(p4, aes(x=Var1, y=Prop_cells, fill= status)) +
  geom_bar(stat="identity", position = "dodge") +
  ylim(0, 30) +
  scale_fill_manual(values = c("darkgoldenrod3", "dodgerblue3")) +
  labs(title = "Proportion of cell types in late disease and prelesion Spleen", x = "Cell types"
       , y = "Proportion, max 30 %") +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 14),
        plot.title = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  NoLegend()

grid.arrange(p1,p2,p3,p4, ncol=2)

library(patchwork)
p1 | p2
p1 
p2
p3
p4

#### Tests hide ####
element_expression_plot <- function(index,elements=TSSs4YCNS_permissive){
  plot <- data.frame(TPM = assay(elements, “TPM”)[index,]) # in the TPM normalized expression assay, select row of the enhancer candidate of interest and save as a column of a dataframe
  plot$term <- elements$term # next to each count add the term of the sample
  plot$organ_system <- elements$organ_system # next to each count add the organ system of the sample
  plot$animal_id<- elements$animal_id
  plot$sample_origin<- elements$sample_origin
  plot$Library_id<- elements$Library_id
  plot <- arrange(plot, organ_system, term) # sort the dataframe by organ system, so that terms of the same organ will be grouped together
  plot$term <- factor(plot$term, levels=unique(plot$term)) # make terms a factor so that above order is kept in plot
  ggplotly(ggplot(plot, aes(x=term, y=TPM), guide_legend(title=“term”)) + # basic plot of expression values per term with title
             geom_jitter(aes(color=sample_origin), width=0.25) +
             labs(title=“TPM normalized expression”, subtitle=index, x=element_blank())+
             scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
             theme(axis.text.x = element_text(size=12,angle = 30))
  )
}

cellorg <- mutate(celltable,
                     org = case_when(
                       contains(celltable$Var1 == "B cells") ~ "B cells",
                       contains(celltable$Var1 == "Plasma cells") ~ "B cells",
                       contains(celltable$Var1 == "Fibro") ~ "Fibroblasts",
                       contains(celltable$Var1 == "Tcells" ) ~ "T cells",
                       contains(celltable$Var1 == "Tregs" ) ~ "T cells",
                       contains(celltable$Var1 == "Teffs" ) ~ "T cells",
                       contains(celltable$Var1 == "EC's") ~ "Endothelial cells",
                       contains(celltable$Var1 == "Macrophages") ~ "Macrophages",
                       contains(celltable$Var1 == "onocytes") ~ "Monocytes",
                       contains(celltable$Var1 == "VSMC") ~ "VSMC's",
                       contains(celltable$Var1 == "DC") ~ "DC's",
                       contains(celltable$Var1 == "ILC") ~ "ILC's and NK cells",
                       contains(celltable$Var1 == "NK") ~ "ILC's and NK cells",
                       contains(celltable$Var1 == "Granulocytes" ) ~ "Granulocytes",
                       contains(celltable$Var1 == "MAST" ) ~ "Granulocytes",
                       contains(celltable$Var1 == "Divi") ~ "Dividing cells"
                     ))
tissue_id <- tissuedata$Tissue
TERVA2_harmony@meta.data$tissue_id <- as.factor(tissue_id)
