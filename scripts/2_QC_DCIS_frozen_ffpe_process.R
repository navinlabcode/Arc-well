#####-----------Process-------######
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
library(copykit)
library(BuenColors)
library(ggpubr)

source("./scripts/0_arcwell_functions.R")
load("./pre_load_data/pre_load_data.rda")
chr_name <- unlist(read.table("./pre_load_data/chr_name.txt"), use.names = F)
chr_color <- read.table("./pre_load_data/chr_color.txt")
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", 
            "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")

setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload")

#-----BCIS28_frozen----
pro_name <- c("bcis28_frozen_qc")
#-----BCIS28_ffpe----
pro_name <- c("bcis28_ffpe_qc")

raw_path <- paste0("./map_seg_output/", pro_name)
filter_cells = read.table(paste0("./metrics/", pro_name, "_filtered_bincounts_newnormal.txt"), header = T)

#----------continue procedure------
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx # 1478/1583

dim(filter_cells)
topleft(filter_cells)
filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] %>% janitor::make_clean_names()
length(filt_cells_names) 
varbin_mtx_tumor <- varbin_mtx[,filt_cells_names]
varbin_mtx_tumor@colData  # 1119/831

#####--NEED TO MODIFY---add meta info------
# ---bcis28_frozen---
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_28t_"), "primary", "recurrence")) %>% dplyr::select(-1)

#---bcis28_ffpe---
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%  
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_28t_"), "primary", "recurrence")) %>% dplyr::select(-1)

#####------processing-------##############
varbin_mtx_tumor@colData <- cbind(varbin_mtx_tumor@colData, name_meta)
table(varbin_mtx_tumor@colData$peak)
table(varbin_mtx_tumor@colData$dispense)
table(varbin_mtx_tumor@colData$timepoint)

varbin_mtx_tumor_log <- logNorm(varbin_mtx_tumor, transform = "log2")
saveRDS(varbin_mtx_tumor_log, file = paste0("./objects/", pro_name, c("_filtered_copykit.rds")))
# varbin_mtx_tumor_log <- readRDS(paste0("./objects/", pro_name, c("_filtered_copykit.rds")))

#----merge BCIS28 frozen and ffpe---to plot heatmap together-------
pro_name <- "bcis28_frozen_ffpe_merged"
varbin_mtx_frozen <- readRDS("./objects/bcis28_frozen_qc_filtered_copykit.rds")
varbin_mtx_ffpe <- readRDS("./objects/bcis28_ffpe_qc_filtered_copykit.rds")
varbin_mtx_frozen@colData$condition <- "Frozen"
varbin_mtx_ffpe@colData$condition <- "FFPE"
varbin_mtx_merged <- cbind(varbin_mtx_frozen, varbin_mtx_ffpe)


#----calculate UMAP-----
set.seed(31)
near_nb <- 20
umap_data <-  data.frame(uwot::umap(log2(t(varbin_mtx_merged@assays@data$segment_ratios)), 
                                    metric = "manhattan", min_dist = 0.1, spread = 3, n_neighbors = near_nb))

umap_data2 <- umap_data %>% dplyr::rename("UMAP_1m" = "X1", "UMAP_2m" = "X2")
varbin_mtx_merged@int_colData$reducedDims <- umap_data2
varbin_mtx_merged@colData <- cbind(varbin_mtx_merged@colData, umap_data2)

ggplot(umap_data2) + geom_point(aes(x = UMAP_1m, y = UMAP_2m))

#---find clusters----
set.seed(17)
hdb_data2 <- dbscan::hdbscan(umap_data2, minPts = nrow(umap_data2)*0.02)
subclones <- paste0("c", as.character(hdb_data2$cluster))

varbin_mtx_merged@colData <- cbind(varbin_mtx_merged@colData, subclones)
varbin_mtx_merged@colData$subclones <- factor(varbin_mtx_merged@colData$subclones, 
                                              levels = paste0("c", 0:(length(table(subclones)))))
varbin_mtx_merged@colData$tp_clstm <- paste(varbin_mtx_merged@colData$condition, varbin_mtx_merged@colData$subclones, sep="_")

if(names(table(varbin_mtx_merged@colData$subclones)[1]) == 0) {my_col = new_pal} else {my_col = c("grey",new_pal)}

p1 <- ggplot(as.data.frame(varbin_mtx_merged@colData),aes(x = UMAP_1m, y = UMAP_2m, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = my_col) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name, "_passInitialQC_cells_umap.pdf"), p1, width = 5, height = 4)

#----remove outlier cluster c0 and cluster that have less 6 cells.
#----paired samples
clone_num <- table(varbin_mtx_merged@colData$subclones, varbin_mtx_merged@colData$condition)
clone_num_less6 <- as.data.frame(clone_num) %>% filter(Var1 != "c0") %>% filter(Freq > 0 & Freq <6) %>% 
  mutate(comb = paste(Var2, Var1, sep = "_"))

varbin_mtx_merged2 <- varbin_mtx_merged[, !(subclones == "c0" | (varbin_mtx_merged@colData$tp_clstm %in% clone_num_less6$comb))]
varbin_mtx_merged2@colData

#----Plot UMAPs-----
p1 <- ggplot(as.data.frame(varbin_mtx_merged2@colData),aes(x = UMAP_1m, y = UMAP_2m, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = new_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap2.pdf"), p1, width = 3.6, height = 2.5)

tp_col <- c("#fb5607", "#1D4E89")
p1 <- ggplot(shuf(as.data.frame(varbin_mtx_merged2@colData)),aes(x = UMAP_1m, y = UMAP_2m, fill = condition)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = tp_col) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap_by_conditions_shuffle2.pdf"), p1, 
                 width = 3.6, height = 2.5)

# pro_name <- "BCIS28_frozen_ffpe_merged"
saveRDS(varbin_mtx_merged2, file = paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
# varbin_mtx_merged2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
#----Complexheatmap----
library(ComplexHeatmap)
mtx_srt <- as.data.frame(varbin_mtx_merged2@colData) %>% 
  arrange(factor(subclones, levels = paste0("c", 1:length(unique(varbin_mtx_merged2@colData$subclones)))))
shuf_name <- mtx_srt %>% group_by(subclones) %>% mutate(sample = sample(sample)) %>% pull(sample)
mtx_srt <- mtx_srt[shuf_name,]
ht_mtx <- log2(t(varbin_mtx_merged2@assays@data$segment_ratios))[mtx_srt$sample,]

#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("peak", "subclones", "condition")) 
rownames(anno_mtx) <- NULL

peak_col <- c("#219ebc","#f4a261")
names(peak_col) <- c("d", "a")
clst_col <- new_pal[1:length(unique(anno_mtx$subclones))]
names(clst_col) <- paste0("c", 1:length(unique(anno_mtx$subclones)))
tp_col <- c("#fb5607", "#1D4E89")
names(tp_col) <- c("FFPE", "Frozen")

ha_row=rowAnnotation(df = anno_mtx, col = list(peak=peak_col, subclones= clst_col, condition=tp_col), 
                     show_annotation_name = F)
#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-2,0,2)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue4", "white", "firebrick4"))

pdf(paste0("./figures/", pro_name, "_complexHeatmap.pdf"), height = 8, width = 8)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        row_split = anno_mtx$subclones,
        name = "scheatmap", show_row_names = F, show_column_names = F, 
        row_title = paste0("single cells: ", nrow(ht_mtx), " (",names(table(anno_mtx$condition)[1]), ": ", 
                           table(anno_mtx$condition)[1], "; ",
                           names(table(anno_mtx$condition)[2]), ": ", table(anno_mtx$condition)[2], ")"),
        column_title = paste0(pro_name, "_scHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                              labels_gp = gpar(fontsize = 12)), top_annotation = ha_col, left_annotation = ha_row)
dev.off()

#---calculate basic matrics------
pro_name <- "bcis28_frozen_ffpe_merged"
varbin_mtx_merged2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
bin_count <- varbin_mtx_merged2@assays@data$bin_counts

#----calculate overdispersion---
bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% rename(over_disp = V1) %>% dplyr::select("over_disp")
varbin_mtx_merged2@colData <- cbind(varbin_mtx_merged2@colData, bin_count_overdisp2) 

saveRDS(varbin_mtx_merged2, file = paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))
# varbin_mtx_merged2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))

merged_meta <- as.data.frame(varbin_mtx_merged2@colData) 
meta_fro <- merged_meta %>% filter(condition == "Frozen") 
meta_form <- merged_meta %>% filter(condition == "FFPE") 
dim(meta_fro)
dim(meta_form)

mean(meta_fro$reads_total)
mean(meta_form$reads_total)
mean(meta_fro$reads_assigned_bins)
mean(meta_form$reads_assigned_bins)
mean(meta_fro$percentage_duplicates)
mean(meta_form$percentage_duplicates)
mean(meta_fro$median_bin_count)
mean(meta_form$median_bin_count)

