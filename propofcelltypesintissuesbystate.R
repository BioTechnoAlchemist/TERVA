# Proportion of cell types in tissues by state
setwd("/scratch/project_2005050/Rstats")

Idents(TERVA2_harmony) <- "celltype"
count_table <- as.data.frame(table(TERVA2_harmony@meta.data$celltype, TERVA2_harmony@meta.data$Sample, TERVA2_harmony@meta.data$group_id))

library(tidyverse)
library(patchwork)

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

p1 | p2
p1 
p2
p3
p4

