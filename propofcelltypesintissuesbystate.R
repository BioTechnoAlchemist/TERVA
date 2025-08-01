# Proportion of cell types in tissues by state
# Author: Lea Mikkola

setwd("/scratch/project_2005050/Rstats")

Idents(TERVA2_harmony) <- "celltype"
count_table <- as.data.frame(table(TERVA2_harmony@meta.data$celltype, TERVA2_harmony@meta.data$Sample, TERVA2_harmony@meta.data$group_id))

library(tidyverse)
library(patchwork)


celltable <- count_table %>%
  group_by(Var2) %>%
  mutate(sum_cells = sum(Freq),
    Prop_cells = Freq/sum_cells*100,
    sum_prop = sum(Prop_cells))

write.csv(celltable, file = "cellproportions.csv", row.names = F, na = "na")
  
#### By tissue in disease states ####

celltable$status <- celltable$Var3

t1 <- celltable %>%
  filter(Var2 %in% c("OBAO", "NOBAO"))
p1 <- ggplot(t1, aes(x=Var1, y=Prop_cells, fill= status)) +
  geom_bar(stat="identity", position = "dodge") +
  ylim(0, 30) +
  scale_fill_manual(values = c("darkgoldenrod3", "dodgerblue3")) +
  labs(title = "Aorta", x = "Cell types"
       , y = "Proportion, max 30 %") +
  theme(plot.title = element_text(size = 16, margin = margin(t=40,b=-30)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size=14)) +
  NoLegend()


t2 <- celltable %>%
  filter(Var2 %in% c("OBeWAT", "NOBeWAT"))
p2 <- ggplot(t2, aes(x=Var1, y=Prop_cells, fill= status)) +
  geom_bar(stat="identity", position = "dodge") +
  ylim(0, 30) +
  scale_fill_manual(values = c("darkgoldenrod3", "dodgerblue3")) +
  labs(title = "eWAT", x = "Cell types"
       , y = "Proportion, max 30 %") +
  theme(plot.title = element_text(size = 16, margin = margin(t=40,b=-30)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size=14)) +
  NoLegend()

t3 <- celltable %>%
  filter(Var2 %in% c("OBPVAT", "NOBPVAT"))
p3 <- ggplot(t3, aes(x=Var1, y=Prop_cells, fill= status)) +
  geom_bar(stat="identity", position = "dodge") +
  ylim(0, 30) +
  scale_fill_manual(values = c("darkgoldenrod3", "dodgerblue3")) +
  labs(title = "PVAT", x = "Cell types"
       , y = "Proportion, max 30 %") +
  theme(plot.title = element_text(size = 16, margin = margin(t=40,b=-30)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size=14)) +
  NoLegend()

t4 <- celltable %>%
  filter(Var2 %in% c("OBSpleen", "NOBSpleen"))
p4 <- ggplot(t4, aes(x=Var1, y=Prop_cells, fill= status)) +
  geom_bar(stat="identity", position = "dodge") +
  ylim(0, 30) +
  scale_fill_manual(values = c("darkgoldenrod3", "dodgerblue3")) +
  labs(title = "Spleen", x = "Cell types"
       , y = "Proportion, max 30 %") +
  theme(axis.text.x=element_text(angle=70, hjust=1, size = 20),
        plot.title = element_text(size = 16, margin = margin(t=40,b=-30)),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size=14)) +
  NoLegend()

p1 / p2 / p3 / p4



