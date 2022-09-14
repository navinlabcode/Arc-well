#####-----------Process-------######
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
library(copykit)

source("./scripts/0_arcwell_functions.R")
load("./pre_load_data/pre_load_data.rda")
chr_name <- unlist(read.table("./pre_load_data/chr_name.txt"), use.names = F)
chr_color <- read.table("./pre_load_data/chr_color.txt")
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", 
            "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")

setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload")

#-----lung 30----
pro_name <- c("lung30_p1")
#-----lung 31----
pro_name <- c("lung31_p2")
#-----prostate ax4_bl----
#---NOTE:: This sample has multiple cells from ax62, need to remove those first!!!
pro_name <- c("prostate_p1_ax4bl")
raw_path <- paste0("./map_seg_output/", pro_name)
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx <- varbin_mtx[,stringr::str_detect(varbin_mtx@colData$sample, "ax4")]
varbin_mtx@colData # 554
#-----prostate pcf169----
pro_name <- c("prostate_p2_pcf169")

#----------processing------
raw_path <- paste0("./map_seg_output/", pro_name)
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx@colData # lung30:1202;  

filter_cells = read.table(paste0("./metrics/", pro_name, "_filtered_bincounts_newnormal.txt"), header = T)
filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] %>% janitor::make_clean_names()
length(filt_cells_names) 
varbin_mtx_tumor <- varbin_mtx[,filt_cells_names]
varbin_mtx_tumor@colData  # lung30: 946;   

#####--NEED TO MODIFY---add meta info------
#---lung 30---
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_s30_"), "primary", "recurrence"))

#---lung 31---
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_s31_"), "primary", "recurrence")) 

#---ax4_bl---
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%  
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_ax4_"), "primary", "recurrence")) 

#---pcf169---
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%
  mutate(peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_pcf"), "primary", "recurrence")) %>% dplyr::select(c("peak", "timepoint"))


############------processing-------##############
varbin_mtx_tumor@colData <- cbind(varbin_mtx_tumor@colData, name_meta)
table(varbin_mtx_tumor@colData$peak)
table(varbin_mtx_tumor@colData$dispense)
table(varbin_mtx_tumor@colData$timepoint)

varbin_mtx_tumor_log <- logNorm(varbin_mtx_tumor, transform = "log2")
saveRDS(varbin_mtx_tumor_log, file = paste0("./objects/", pro_name, c("_filtered_copykit.rds")))
# varbin_mtx_tumor_log <- readRDS(paste0("./objects/", pro_name, c("_filtered_copykit.rds")))

#----UMAP-----
set.seed(31)
near_nb <- 20
umap_data <-  data.frame(uwot::umap(log2(t(varbin_mtx_tumor_log@assays@data$segment_ratios)), 
                                    metric = "manhattan", min_dist = 0.1, spread = 3, n_neighbors = near_nb))

umap_data2 <- umap_data %>% dplyr::rename("UMAP_1" = "X1", "UMAP_2" = "X2")
varbin_mtx_tumor_log@int_colData$reducedDims <- umap_data2
varbin_mtx_tumor_log@colData <- cbind(varbin_mtx_tumor_log@colData, umap_data2)

ggplot(umap_data2) + geom_point(aes(x = UMAP_1, y = UMAP_2))


#---find clusters----
set.seed(17)
hdb_data <- dbscan::hdbscan(umap_data2, minPts = nrow(umap_data2)*0.015)
subclones <- paste0("c", as.character(hdb_data$cluster))

varbin_mtx_tumor_log@colData <- cbind(varbin_mtx_tumor_log@colData, subclones)
varbin_mtx_tumor_log@colData$subclones <- factor(varbin_mtx_tumor_log@colData$subclones, 
                                                 levels = paste0("c", 0:(length(table(subclones))-1)))

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = c("grey",new_pal)) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name, "_passInitialQC_cells_umap.pdf"), p1, width = 5, height = 4)

#----remove outlier cluster c0 and cluster that have less 3 cells.
#----paired samples
clone_num <- table(varbin_mtx_tumor_log@colData$subclones, varbin_mtx_tumor_log@colData$timepoint)
clone_num_less3 <- as.data.frame(clone_num) %>% filter(Var1 != "c0") %>% filter(Freq > 0 & Freq <4)

varbin_mtx_tumor_log2 <- varbin_mtx_tumor_log[, !(subclones == "c0" | 
                                                    ((varbin_mtx_tumor_log@colData$subclones %in% clone_num_less3$Var1) & 
                                                       (varbin_mtx_tumor_log@colData$timepoint %in% clone_num_less3$Var2)))]
varbin_mtx_tumor_log2@colData

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = new_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap.pdf"), p1, width = 3.05, height = 2)

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
varbin_mtx_tumor_log2@colData

#----Number of normal cells-----
filter_mtx <- read.table(paste0("./metrics/", pro_name, "_metadata.metrics_newnormal.txt"), header = T) 
rownames(filter_mtx) <- filter_mtx$sample %>% janitor::make_clean_names() 
filter_mtx2 <- filter_mtx %>% mutate(dups_percentage = dups_removed/total_reads) %>% rownames_to_column() %>% 
  dplyr::select(c("rowname", "filter_corr_value","reads_kept","median_bin_count","filtered","is_normal"))
table(filter_mtx2$filtered, filter_mtx2$is_normal)

#----single cell heatmap----
library(ComplexHeatmap)
mtx_srt <- as.data.frame(varbin_mtx_tumor_log2@colData) %>% 
  arrange(factor(subclones, levels = paste0("c", 1:length(unique(varbin_mtx_tumor_log2@colData$subclones)))))
head(mtx_srt)
ht_mtx <- log2(t(varbin_mtx_tumor_log2@assays@data$segment_ratios))[mtx_srt$sample,]
topleft(ht_mtx)

#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("peak", "subclones")) 
rownames(anno_mtx) <- NULL

peak_col <- c("#219ebc","#f4a261")
names(peak_col) <- c("d", "a")
clst_col <- new_pal[1:length(unique(anno_mtx$subclones))]
names(clst_col) <- paste0("c", 1:length(unique(anno_mtx$subclones)))

ha_row=rowAnnotation(df = anno_mtx, col = list(peak=peak_col, subclones= clst_col), show_annotation_name = F)

#-----annotated genes---
gene_bin <- read.table("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/hg19_gene_binpos_map.tsv", header = 1)
#---lung--
keep_gene <- c("KRAS","MYC","TP53","ERBB2","EGFR", "FGFR1","FGFR2","ALK","MET","ROS1",
               "NTRK1","RET","DDR2","PDGFRA","NF1","BRAF","MAP2K1","NOTCH1","KMT2D",
               "EZH2","TET2","DNMT3A","SOX2","KEAP1","CDKN2A","NRG1","STK11","PTEN")
#---prostate---
keep_gene <- c("MDM4", "PIK3CA", "PIK3CB", "PIK3R1", "CDK6", "MET", "BRAF", 
               "NKX3-1", "MYC", "PTEN", "ATM", "BRCA2", "RB1", "MSH2", "MSH6")

gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% keep_gene)
ha_bottom_cs = columnAnnotation(clonal_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                         side = "bottom", labels_gp = gpar(fontsize = 14)))

#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue4", "white", "firebrick4"))

pdf(paste0("./figures/", pro_name, "_complexHeatmap2.pdf"), height = 8.5, width = 8)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        row_split = anno_mtx$hdb_clst,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0("single cells:", nrow(ht_mtx)),
        column_title = paste0(pro_name, "_scHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, left_annotation = ha_row, bottom_annotation = ha_bottom_cs)
dev.off()


#----calculate overdispersion---
bin_count <- varbin_mtx_tumor_log2@assays@data$bin_counts
bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% rename(over_disp = V1) %>% dplyr::select("over_disp")
varbin_mtx_tumor_log2@colData <- cbind(varbin_mtx_tumor_log2@colData, bin_count_overdisp2) 
saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))

#---calculate basic matrices------
#---reads, pcr dup, and bin for aneuploid cells
# pro_name <- "lung30_p1"
# pro_name <- "lung31_p2"
# pro_name <- "prostate_p1_ax4bl"
# pro_name <- "prostate_p2_pcf169"
# varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/",pro_name, c("_final_filtered_meta_overdisp_copykit.rds")))

varbin_mtx_tumor_log2@colData
mean(varbin_mtx_tumor_log2@colData$reads_total)
mean(varbin_mtx_tumor_log2@colData$reads_assigned_bins)
mean(varbin_mtx_tumor_log2@colData$percentage_duplicates)
mean(varbin_mtx_tumor_log2@colData$median_bin_count)


#-----Overdipsersion Figure for Supplementary Figure 5a------
all_name <- c("lung30_p1", "lung31_p2", "prostate_p1_ax4bl", "prostate_p2_pcf169")

lung_prostate_merged <- data.frame()
for (pro_name in all_name) {
  varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_meta_overdisp_copykit.rds")))
  mtx_tmp <- as.data.frame(varbin_mtx_tumor_log2@colData) %>% 
    dplyr::select(c("sample","reads_total","reads_assigned_bins","percentage_duplicates",
                    "median_bin_count","peak","subclones","over_disp")) %>% dplyr::mutate(tissue = pro_name)
  lung_prostate_merged <- rbind(lung_prostate_merged, mtx_tmp)
}
write.csv(lung_prostate_merged, "./metrics/lung_prostate_merged_meta.csv")
# lung_prostate_merged <- read.csv("./metrics/lung_prostate_merged_meta.csv", row.names = 1)

library(ggpubr)
tp_col <- c("#4cc9f0","#8cb369","#ffb703","#e76f51")
p1 <- ggplot(lung_prostate_merged, aes(x = tissue, y = over_disp, fill = tissue)) + geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = tp_col) + 
  theme_classic2() + ylim(0, 0.1) + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) 
cowplot::ggsave2("./figures/lung_prostate_merged_over_disp_boxplot.pdf", p1, width = 2.2, height = 3)



