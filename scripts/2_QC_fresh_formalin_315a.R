#####-----------Process-------######
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
library(copykit)
library(stringr)
library(ComplexHeatmap)

source("./scripts/0_arcwell_functions.R")
load("./pre_load_data/pre_load_data.rda")
chr_name <- unlist(read.table("./pre_load_data/chr_name.txt"), use.names = F)
chr_color <- read.table("./pre_load_data/chr_color.txt")
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", 
            "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")

setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload")

#-----315A_fresh----
pro_name <- c("fresh_315a")
#-----315A_formalin----
pro_name <- c("formalin_315a")

raw_path <- paste0("./map_seg_output/", pro_name)
filter_cells = read.table(paste0("./metrics/", pro_name, "_filtered_bincounts_newnormal.txt"), header = T)

#----------continue procedure------
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx # 2164/1660

filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] %>% janitor::make_clean_names()
length(filt_cells_names) 
varbin_mtx_tumor <- varbin_mtx[,filt_cells_names]
varbin_mtx_tumor@colData  # 1932/1272

#####--NEED TO MODIFY---add meta info------
# ---fresh_315a---
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         cellstate = stringr::str_extract(my_name, "nofix_(.*?)_w"), 
         cellrow = str_extract(my_name, "r\\d+"), cellcol = str_extract(my_name, "c\\d+")) %>%
  mutate(cellstate = str_sub(cellstate, 7, -3), cellrow = str_sub(cellrow, 2), cellcol = str_sub(cellcol, 2)) %>% dplyr::select(-1)

#---formalin_315a---
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         cellstate = stringr::str_extract(my_name, "yes_(.*?)_s"), 
         cellrow = str_extract(my_name, "r\\d+"), cellcol = str_extract(my_name, "c\\d+")) %>% 
  mutate(cellstate = str_sub(cellstate, 8, -3), cellrow = str_sub(cellrow, 2), cellcol = str_sub(cellcol, 2)) %>% dplyr::select(-1)

############------processing-------##############
varbin_mtx_tumor@colData <- cbind(varbin_mtx_tumor@colData, name_meta)
varbin_mtx_tumor_log <- logNorm(varbin_mtx_tumor, transform = "log2")

saveRDS(varbin_mtx_tumor_log, file = paste0("./objects/", pro_name, c("_filtered_copykit.rds")))
# varbin_mtx_tumor_log <- readRDS(paste0("./objects/", pro_name, c("_filtered_copykit.rds")))

#----merge fresh and formalin 315A ---to plot heatmap together-------
varbin_mtx_fresh <- readRDS("./objects/fresh_315a_filtered_copykit.rds")
varbin_mtx_formalin <- readRDS("./objects/formalin_315a_filtered_copykit.rds")

#---header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

my_seg <- log2(rbind(t(varbin_mtx_fresh@assays@data$segment_ratios), t(varbin_mtx_formalin@assays@data$segment_ratios)))
my_seg <- my_seg[sample(rownames(my_seg)),]
name_meta <- rownames(my_seg) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%
  mutate(sample = ifelse(stringr::str_detect(my_name, "_nofix_"), "fresh_315a", "formalin_315a")) %>% dplyr::select("sample")

sample_col <- c("#fb5607", "#1D4E89")
names(sample_col) <- c("formalin_315a", "fresh_315a")
ha_row=rowAnnotation(df = name_meta, col = list(sample=sample_col), show_annotation_name = F)

breaks = c(-2,0,2)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue4", "white", "firebrick4"))
pdf(paste0("./figures/fresh_formalin_315a_combined_complexHeatmap_unclustered.pdf"), height = 8, width = 8)
Heatmap(my_seg, cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, name = "scheatmap", 
        show_row_names = F, show_column_names = F, row_title = paste0("single cells:", nrow(my_seg)), 
        column_title = paste0("315A_scHeatmap"), use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row)
dev.off()

meta_merged <- as.data.frame(rbind(varbin_mtx_fresh@colData, varbin_mtx_formalin@colData))
rownames(meta_merged) <- NULL
meta_merged$condition = ifelse(stringr::str_detect(meta_merged$sample, "_nofix_"), "fresh_315a", "formalin_315a")
table(meta_merged$condition)

#---merge copykit object---#######
varbin_mtx_fresh@colData$condition <- "fresh_315a"
varbin_mtx_formalin@colData$condition <- "formalin_315a"
varbin_mtx_merged <- cbind(varbin_mtx_fresh, varbin_mtx_formalin)
bin_count <- varbin_mtx_merged@assays@data$bin_counts

#----calculate overdispersion----
bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% dplyr::rename(over_disp = V1) %>% dplyr::select("over_disp")
varbin_mtx_merged@colData <- cbind(varbin_mtx_merged@colData, bin_count_overdisp2)

saveRDS(varbin_mtx_merged, file = "./objects/fresh_formalin_315a_merged_final_filtered_overdisp_copykit.rds")
# varbin_mtx_merged <- readRDS("./objects/fresh_formalin_315a_merged_final_filtered_overdisp_copykit.rds")

#-----QC matrics-----
merged_meta <- as.data.frame(varbin_mtx_merged@colData)
meta_fro <- merged_meta %>% filter(condition == "fresh_315a")
meta_form <- merged_meta %>% filter(condition == "formalin_315a")

mean(meta_fro$reads_total)
mean(meta_form$reads_total)
mean(meta_fro$reads_assigned_bins)
mean(meta_form$reads_assigned_bins)
mean(meta_fro$percentage_duplicates)
mean(meta_form$percentage_duplicates)
mean(meta_fro$median_bin_count)
mean(meta_form$median_bin_count)

