setwd("/scratch/project_2005050/Rstats")
.libPaths(c("/projappl/project_2005050/project_rpackages_4.1.1", .libPaths()))

library(tidyverse)
library(immunarch)

# 1.2) Replace with the path to your processed 10x data or to the clonotypes file
file_path = "/scratch/project_2005050/Rstats/filtered_contigs/LDAO/TCR"

# 1.3) Load 10x data with repLoad
immdata_10x <- repLoad(file_path)
