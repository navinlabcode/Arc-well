# Arc-well
This repository contains the scripts, meta data metrics files used in the manuscript: Archival single cell sequencing reveals persistent subclones over years to decades of DCIS progression.


#### scripts 
R scripts in this directory includes all of the codes that could reproduce figures of Arc-well paper.
- _snakemake_files_: includes the workflow of downsampling coverage calculations
- _CNA_pipeline_: includes the workflow of initial QC and variable binning, prepares the inputs for downstream analysis, can be found [here](https://github.com/navinlabcode/CNV_pipeline).

#### pre_load_data 
files that required during excuting the R scripts under scripts folder
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
```

### _Data source_
------------
The original sequencing data from this study has been deposited to the Sequence Read Archive (SRA): PRJNA799605.

### _Contact_
------------
For more detailed information, please [email](mailto:nnavin@mdanderson.org) corresponding author.

