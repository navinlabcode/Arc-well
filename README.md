# Arc-well
This repository contains the scripts, meta data metrics files used in the manuscript: Archival single cell sequencing reveals persistent subclones over years to decades of DCIS progression.


#### scripts 
R scripts in this directory includes all of the codes that could reproduce figures of Arc-well paper.
- _snakemake_files_: includes the workflow of downsampling coverage calculations
- _CNA_pipeline_: includes the workflow of initial QC and variable binning, prepares the inputs for downstream analysis, can be found [here](https://github.com/navinlabcode/CNV_pipeline).

#### pre_load_data 
data required for the R scripts under scripts folder

#### metrics  
includes filtering files, QC metrics, clinical, computational meta data from cell lines, FFPE samples.

### _Dependencies_
------------
R scripts are dependent on [*CopyKit*](https://github.com/navinlabcode/copykit)(*v0.1.0*), which can be installed by:
``` r
devtools::install_github(repo = "navinlabcode/copykit",ref="f709a48")
```
Session info:
```
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggtree_3.2.1                            ape_5.6-2                               dendextend_1.15.2                      
 [4] dbscan_1.1-10                           copykit_0.1.0                           DNAcopy_1.68.0                         
 [7] Rsubread_2.8.1                          SingleCellExperiment_1.16.0             SummarizedExperiment_1.24.0            
[10] MatrixGenerics_1.6.0                    matrixStats_0.61.0                      ggalt_0.4.0                            
[13] DEGreport_1.33.1                        ggpubr_0.4.0                            Homo.sapiens_1.3.1                     
[16] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 org.Hs.eg.db_3.14.0                     GO.db_3.14.0                           
[19] OrganismDbi_1.36.0                      GenomicFeatures_1.46.5                  GenomicRanges_1.46.1                   
[22] GenomeInfoDb_1.30.1                     AnnotationDbi_1.56.2                    IRanges_2.28.0                         
[25] S4Vectors_0.32.3                        Biobase_2.54.0                          BiocGenerics_0.40.0                    
[28] RColorBrewer_1.1-2                      ComplexHeatmap_2.10.0                   forcats_0.5.1                          
[31] stringr_1.4.0                           purrr_0.3.4                             readr_2.1.2                            
[34] tidyverse_1.3.1                         tibble_3.1.6                            useful_1.2.6                           
[37] ggplot2_3.3.5                           cowplot_1.1.1                           tidyr_1.2.0                            
[40] dplyr_1.0.8                            

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3              rtracklayer_1.54.0          ggthemes_4.2.4              prabclus_2.3-2             
  [5] bit64_4.0.5                 knitr_1.37                  DelayedArray_0.20.0         data.table_1.14.2          
  [9] KEGGREST_1.34.0             RCurl_1.98-1.6              doParallel_1.0.17           generics_0.1.2             
 [13] RSQLite_2.2.10              bit_4.0.4                   tzdb_0.2.0                  httpuv_1.6.5               
 [17] xml2_1.3.3                  lubridate_1.8.0             assertthat_0.2.1            viridis_0.6.2              
 [21] amap_0.8-18                 xfun_0.30                   hms_1.1.1                   promises_1.2.0.1           
 [25] DEoptimR_1.0-10             fansi_1.0.2                 restfulr_0.0.13             progress_1.2.2             
 [29] dbplyr_2.1.1                readxl_1.3.1                igraph_1.2.11               DBI_1.1.2                  
 [33] geneplotter_1.72.0          tmvnsim_1.0-2               reshape_0.8.8               ellipsis_0.3.2             
 [37] ggnewscale_0.4.6            backports_1.4.1             annotate_1.72.0             biomaRt_2.50.3             
 [41] ggalluvial_0.12.3           vctrs_0.3.8                 abind_1.4-5                 cachem_1.0.6               
 [45] withr_2.4.3                 robustbase_0.93-9           lasso2_1.2-22               GenomicAlignments_1.30.0   
 [49] treeio_1.18.1               prettyunits_1.1.1           mclust_5.4.9                mnormt_2.0.2               
 [53] cluster_2.1.2               segmented_1.4-0             lazyeval_0.2.2              crayon_1.5.0               
 [57] genefilter_1.76.0           edgeR_3.36.0                pkgconfig_2.0.3             nlme_3.1-155               
 [61] vipor_0.4.5                 nnet_7.3-17                 pals_1.7                    diptest_0.76-0             
 [65] rlang_1.0.1                 miniUI_0.1.1.1              lifecycle_1.0.1             filelock_1.0.2             
 [69] extrafontdb_1.0             BiocFileCache_2.2.1         modelr_0.1.8                dichromat_2.0-0            
 [73] ggrastr_1.0.1               cellranger_1.1.0            graph_1.72.0                Matrix_1.4-0               
 [77] aplot_0.1.2                 carData_3.0-5               reprex_2.0.1                beeswarm_0.4.0             
 [81] GlobalOptions_0.1.2         png_0.1-7                   viridisLite_0.4.0           rjson_0.2.21               
 [85] bitops_1.0-7                ConsensusClusterPlus_1.58.0 KernSmooth_2.23-20          Biostrings_2.62.0          
 [89] blob_1.2.2                  shape_1.4.6                 rstatix_0.7.0               gridGraphics_0.5-1         
 [93] ggsignif_0.6.3              scales_1.1.1                memoise_2.0.1               magrittr_2.0.2             
 [97] plyr_1.8.6                  zlibbioc_1.40.0             compiler_4.1.2              BiocIO_1.4.0               
[101] ash_1.0-15                  clue_0.3-60                 DESeq2_1.34.0               Rsamtools_2.10.0           
[105] snakecase_0.11.0            cli_3.2.0                   XVector_0.34.0              patchwork_1.1.1            
[109] MASS_7.3-55                 tidyselect_1.1.2            stringi_1.7.6               proj4_1.0-11               
[113] yaml_2.3.5                  locfit_1.5-9.4              ggrepel_0.9.1               tools_4.1.2                
[117] parallel_4.1.2              rio_0.5.29                  circlize_0.4.14             rstudioapi_0.13            
[121] bluster_1.4.0               foreach_1.5.2               foreign_0.8-82              logging_0.10-108           
[125] janitor_2.1.0               gridExtra_2.3               digest_0.6.29               BiocManager_1.30.16        
[129] shiny_1.7.1                 fpc_2.2-9                   Rcpp_1.0.8.1                car_3.0-11                 
[133] broom_0.7.12                later_1.3.0                 httr_1.4.2                  ggdendro_0.1.23            
[137] psych_2.1.9                 kernlab_0.9-29              colorspace_2.0-3            rvest_1.0.2                
[141] XML_3.99-0.9                fs_1.5.2                    splines_4.1.2               uwot_0.1.11                
[145] yulab.utils_0.0.4           RBGL_1.70.0                 tidytree_0.3.8              flexmix_2.3-17             
[149] mapproj_1.2.8               ggplotify_0.1.0             xtable_1.8-4                jsonlite_1.8.0             
[153] modeltools_0.2-23           ggfun_0.0.5                 R6_2.5.1                    mime_0.12                  
[157] htmltools_0.5.2             pillar_1.7.0                glue_1.6.2                  fastmap_1.1.0              
[161] BiocParallel_1.28.3         BiocNeighbors_1.12.0        class_7.3-20                codetools_0.2-18           
[165] maps_3.4.0                  utf8_1.2.2                  lattice_0.20-45             mixtools_1.2.0             
[169] curl_4.3.2                  ggbeeswarm_0.6.0            gtools_3.9.2                zip_2.2.0                  
[173] openxlsx_4.2.5              Rttf2pt1_1.3.10             survival_3.3-0              limma_3.50.1               
[177] munsell_0.5.0               GetoptLong_1.0.5            fastcluster_1.2.3           GenomeInfoDbData_1.2.7     
[181] iterators_1.0.14            haven_2.4.3                 reshape2_1.4.4              gtable_0.3.0               
[185] extrafont_0.17                     
```

### _Data source_
------------
The original sequencing data from this study has been deposited to the Sequence Read Archive (SRA): PRJNA799605.

### _Contact_
------------
For any additional information, please [email](mailto:nnavin@mdanderson.org) corresponding author.

