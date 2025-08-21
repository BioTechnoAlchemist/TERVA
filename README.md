Session information of the single-cell and spatial analyses in R for the TERVA project (Lea Mikkola).

Single-cell analyses:

R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 8.6 (Green Obsidian)

Matrix products: default
BLAS/LAPACK: /opt/intel/oneapi/mkl/2022.1.0/lib/intel64/libmkl_gf_lp64.so.2

locale:
[1] C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] sp_1.5-0               SeuratObject_4.1.0     Seurat_4.1.1           SeuratDisk_0.0.0.9020  EnhancedVolcano_1.14.0 ggrepel_0.9.1          colorBlindness_0.1.9   harmony_0.1.0          Rcpp_1.0.9             clustree_0.5.0        
[11] ggraph_2.0.5           factoextra_1.0.7       RColorBrewer_1.1-3     gridExtra_2.3          viridis_0.6.2          viridisLite_0.4.0      cowplot_1.1.1          patchwork_1.1.1        forcats_0.5.1          stringr_1.4.0         
[21] dplyr_1.0.9            purrr_0.3.4            readr_2.1.2            tidyr_1.2.0            tibble_3.1.8           ggplot2_3.3.6          tidyverse_1.3.1       

loaded via a namespace (and not attached):
  [1] pacman_0.5.1                  utf8_1.2.2                    reticulate_1.25               tidyselect_1.1.2              RSQLite_2.2.14                AnnotationDbi_1.60.2          htmlwidgets_1.5.4             grid_4.2.1                   
  [9] Rtsne_0.16                    munsell_0.5.0                 codetools_0.2-18              ica_1.0-2                     future_1.26.1                 miniUI_0.1.1.1                withr_2.5.0                   spatstat.random_2.2-0        
 [17] colorspace_2.0-3              progressr_0.10.1              Biobase_2.58.0                filelock_1.0.2                rstudioapi_0.13               stats4_4.2.1                  ROCR_1.0-11                   tensor_1.5                   
 [25] listenv_0.8.0                 MatrixGenerics_1.8.1          GenomeInfoDbData_1.2.9        polyclip_1.10-0               farver_2.1.0                  bit64_4.0.5                   parallelly_1.32.0             vctrs_0.4.1                  
 [33] generics_0.1.2                BiocFileCache_2.4.0           R6_2.5.1                      GenomeInfoDb_1.34.9           graphlayouts_0.8.0            hdf5r_1.3.5                   gridGraphics_0.5-1            bitops_1.0-7                 
 [41] spatstat.utils_2.3-1          cachem_1.0.6                  DelayedArray_0.22.0           assertthat_0.2.1              promises_1.2.0.1              scales_1.2.0                  rgeos_0.5-9                   gtable_0.3.0                 
 [49] globals_0.15.1                goftest_1.2-3                 tidygraph_1.2.1               rlang_1.0.4                   splines_4.2.1                 lazyeval_0.2.2                spatstat.geom_2.4-0           broom_0.8.0                  
 [57] BiocManager_1.30.18           yaml_2.3.5                    reshape2_1.4.4                abind_1.4-5                   modelr_0.1.8                  backports_1.4.1               httpuv_1.6.5                  tools_4.2.1                  
 [65] ellipsis_0.3.2                spatstat.core_2.4-4           BiocGenerics_0.44.0           ggridges_0.5.3                plyr_1.8.7                    sparseMatrixStats_1.8.0       zlibbioc_1.44.0               RCurl_1.98-1.7               
 [73] rpart_4.1.16                  deldir_1.0-6                  pbapply_1.5-0                 S4Vectors_0.36.2              zoo_1.8-10                    SummarizedExperiment_1.26.1   haven_2.5.0                   cluster_2.1.3                
 [81] fs_1.5.2                      magrittr_2.0.3                data.table_1.14.2             scattermore_0.8               lmtest_0.9-40                 reprex_2.0.1                  RANN_2.6.1                    fitdistrplus_1.1-8           
 [89] matrixStats_0.62.0            hms_1.1.1                     mime_0.12                     xtable_1.8-4                  readxl_1.4.0                  IRanges_2.32.0                compiler_4.2.1                KernSmooth_2.23-20           
 [97] crayon_1.5.1                  htmltools_0.5.3               mgcv_1.8-40                   later_1.3.0                   tzdb_0.3.0                    lubridate_1.8.0               DBI_1.1.3                     tweenr_1.0.2                 
[105] ExperimentHub_2.4.0           dbplyr_2.2.1                  MASS_7.3-57                   rappdirs_0.3.3                Matrix_1.4-1                  cli_3.3.0                     parallel_4.2.1                igraph_1.3.2                 
[113] GenomicRanges_1.48.0          pkgconfig_2.0.3               plotly_4.10.0                 spatstat.sparse_2.1-1         xml2_1.3.3                    XVector_0.38.0                rvest_1.0.2                   digest_0.6.29                
[121] sctransform_0.3.3             RcppAnnoy_0.0.19              spatstat.data_2.2-0           Biostrings_2.66.0             cellranger_1.1.0              leiden_0.4.2                  uwot_0.1.11                   DelayedMatrixStats_1.18.0    
[129] curl_4.3.2                    shiny_1.7.1                   lifecycle_1.0.1               nlme_3.1-157                  jsonlite_1.8.0                fansi_1.0.3                   pillar_1.8.0                  lattice_0.20-45              
[137] KEGGREST_1.38.0               fastmap_1.1.0                 httr_1.4.3                    survival_3.3-1                interactiveDisplayBase_1.34.0 glue_1.6.2                    remotes_2.4.2                 celldex_1.6.0                
[145] png_0.1-7                     BiocVersion_3.15.2            bit_4.0.4                     ggforce_0.3.3                 stringi_1.7.8                 blob_1.2.3                    AnnotationHub_3.4.0           memoise_2.0.1                
[153] irlba_2.3.5                   future.apply_1.9.0     
 
Spatial analyses:

R version 4.3.0 (2023-04-21)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.6.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Europe/Helsinki
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] clustree_0.5.0     ggraph_2.1.0       viridis_0.6.4      viridisLite_0.4.2  hdf5r_1.3.8        patchwork_1.1.3    ggplot2_3.4.3     
 [8] Seurat_4.4.0       dplyr_1.1.3        Matrix_1.6-1.1     SeuratObject_4.1.4 sp_2.1-0          

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3     jsonlite_1.8.7         magrittr_2.0.3         spatstat.utils_3.0-3   farver_2.1.1           rmarkdown_2.25        
  [7] vctrs_0.6.3            ROCR_1.0-11            spatstat.explore_3.2-3 htmltools_0.5.6        sctransform_0.4.0      parallelly_1.36.0     
 [13] KernSmooth_2.23-20     htmlwidgets_1.6.2      ica_1.0-3              plyr_1.8.9             plotly_4.10.2          zoo_1.8-12            
 [19] igraph_1.5.1           mime_0.12              lifecycle_1.0.3        pkgconfig_2.0.3        R6_2.5.1               fastmap_1.1.1         
 [25] fitdistrplus_1.1-11    future_1.33.0          shiny_1.7.5            digest_0.6.33          colorspace_2.1-0       tensor_1.5            
 [31] irlba_2.3.5.1          progressr_0.14.0       fansi_1.0.4            spatstat.sparse_3.0-2  httr_1.4.7             polyclip_1.10-6       
 [37] abind_1.4-5            compiler_4.3.0         bit64_4.0.5            withr_2.5.1            ggforce_0.4.1          MASS_7.3-58.4         
 [43] tools_4.3.0            lmtest_0.9-40          httpuv_1.6.11          future.apply_1.11.0    goftest_1.2-3          glue_1.6.2            
 [49] nlme_3.1-162           promises_1.2.1         grid_4.3.0             Rtsne_0.16             cluster_2.1.4          reshape2_1.4.4        
 [55] generics_0.1.3         gtable_0.3.4           spatstat.data_3.0-1    tidyr_1.3.0            data.table_1.14.8      tidygraph_1.2.3       
 [61] utf8_1.2.3             spatstat.geom_3.2-5    RcppAnnoy_0.0.21       ggrepel_0.9.3          RANN_2.6.1             pillar_1.9.0          
 [67] stringr_1.5.0          later_1.3.1            splines_4.3.0          tweenr_2.0.2           lattice_0.21-8         survival_3.5-5        
 [73] bit_4.0.5              deldir_1.0-9           tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2          knitr_1.44            
 [79] gridExtra_2.3          scattermore_1.2        xfun_0.40              graphlayouts_1.0.1     matrixStats_1.0.0      stringi_1.7.12        
 [85] lazyeval_0.2.2         yaml_2.3.7             pacman_0.5.1           evaluate_0.22          codetools_0.2-19       tibble_3.2.1          
 [91] cli_3.6.1              uwot_0.1.16            xtable_1.8-4           reticulate_1.32.0      munsell_0.5.0          Rcpp_1.0.11           
 [97] globals_0.16.2         spatstat.random_3.1-6  png_0.1-8              parallel_4.3.0         ellipsis_0.3.2         listenv_0.9.0         
[103] scales_1.2.1           ggridges_0.5.4         leiden_0.4.3           purrr_1.0.2            rlang_1.1.1            cowplot_1.1.1         
