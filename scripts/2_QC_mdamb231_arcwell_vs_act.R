#####-----------Process-------######
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
library(copykit)
library(stringr)
library(ComplexHeatmap)
library(tibble)
library(BuenColors)
library(ggpubr)

source("./scripts/arcwell_functions.R")
load("./pre_load_data/pre_load_data.rda")
chr_name <- unlist(read.table("./pre_load_data/chr_name.txt"), use.names = F)
chr_color <- read.table("./pre_load_data/chr_color.txt")
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", 
            "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")

setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload")

pro_name <- c("arcwell_act")
#-----Arc-well----
pro_name1 <- c("arc_well_mda231")
raw_path <- paste0("./map_seg_output/", pro_name1)
filter_cells = read.table(paste0("./metrics/", pro_name1, "_filtered_bincounts_newnormal.txt"), header = T)
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx # 1205

filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] %>% janitor::make_clean_names()
length(filt_cells_names) 
varbin_mtx_tumor <- varbin_mtx[,filt_cells_names]
varbin_mtx_tumor@colData  # 1054/

#-----ACT----
pro_name2 <- c("act_mda231")
raw_path2 <- paste0("./map_seg_output/", pro_name2)
filter_cells2 = read.table(paste0("./metrics/", pro_name2, "_filtered_bincounts_newnormal.txt"), header = T)
varbin_mtx2 <- readVarbinCNA(raw_path2, remove_Y = TRUE)
varbin_mtx2 # 1261

filt_cells_names2 <- colnames(filter_cells2)[4:ncol(filter_cells2)] %>% janitor::make_clean_names()
length(filt_cells_names2) 
varbin_mtx_tumor2 <- varbin_mtx2[,filt_cells_names2]
varbin_mtx_tumor2@colData  # 1213
#----subset MDA231 cell only from ACT---
filt_cells_names3 <- varbin_mtx_tumor2@colData %>% as.data.frame() %>% 
  dplyr::filter(stringr::str_detect(sample, "mdamb231")) %>% dplyr::pull(sample)
varbin_mtx_tumor2 <- varbin_mtx_tumor2[,filt_cells_names3]
varbin_mtx_tumor2@colData # 940

#------merge Arc-well and ACT-----
varbin_mtx_tumor@colData$condition <- "arc_well_mda231"
varbin_mtx_tumor2@colData$condition <- "act_mda231"
varbin_mtx_merged <- cbind(varbin_mtx_tumor, varbin_mtx_tumor2)

varbin_mtx_tumor_log <- logNorm(varbin_mtx_merged, transform = "log2")
saveRDS(varbin_mtx_tumor_log, file = paste0("objects/", pro_name, c("_filtered_copykit.rds")))
# varbin_mtx_tumor_log <- readRDS(paste0("objects/", pro_name, c("_filtered_copykit.rds")))

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
hdb_data <- dbscan::hdbscan(umap_data2, minPts = nrow(umap_data2)*0.02)
subclones <- paste0("c", as.character(hdb_data$cluster))

varbin_mtx_tumor_log@colData <- cbind(varbin_mtx_tumor_log@colData, subclones)
varbin_mtx_tumor_log@colData$subclones <- factor(varbin_mtx_tumor_log@colData$subclones, 
                                                 levels = paste0("c", 0:(length(table(subclones))-1)))
varbin_mtx_tumor_log@colData$tp_clst <- paste(varbin_mtx_tumor_log@colData$condition, 
                                              varbin_mtx_tumor_log@colData$subclones, sep="_")

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = c("grey",new_pal)) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name, "_passInitialQC_cells_umap.pdf"), p1, width = 5, height = 4)

#----remove outlier cluster c0 and cluster that have less 3 cells.
#----paired samples
clone_num <- table(varbin_mtx_tumor_log@colData$subclones, varbin_mtx_tumor_log@colData$condition)
clone_num_less3 <- as.data.frame(clone_num) %>% filter(Var1 != "c0") %>% filter(Freq > 0 & Freq <4) %>% 
  mutate(comb = paste(Var2, Var1, sep = "_"))

varbin_mtx_tumor_log2 <- varbin_mtx_tumor_log[, !(subclones == "c0" | (varbin_mtx_tumor_log@colData$tp_clst %in% clone_num_less3$comb))]
varbin_mtx_tumor_log2@colData

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = new_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap.pdf"), p1, width = 5, height = 4)

tp_col <- c("#C27739", "#3A85C2")
p2 <- ggplot(shuf(as.data.frame(varbin_mtx_tumor_log2@colData)),aes(x = UMAP_1, y = UMAP_2, fill = condition)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = tp_col) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(), 
        legend.text = element_text(size = 7))
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap2.pdf"), p2, width = 5, height = 4)

p3 <- plot_grid(plotlist=list(p1,p2), ncol=2, align='h')
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap_merge.pdf"), p3, width = 7.8, height = 2.8)

saveRDS(varbin_mtx_tumor_log2, file = paste0("objects/", pro_name, c("_final_filtered_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("objects/", pro_name, c("_final_filtered_copykit.rds")))
#----plot heatmap ----
library(ComplexHeatmap)
mtx_srt <- as.data.frame(varbin_mtx_tumor_log2@colData) %>% 
  arrange(factor(subclones, levels = paste0("c", 1:length(unique(varbin_mtx_tumor_log2@colData$subclones)))))
shuf_name <- mtx_srt %>% group_by(subclones) %>% mutate(sample = sample(sample)) %>% pull(sample)
mtx_srt <- mtx_srt[shuf_name,]
ht_mtx <- log2(t(varbin_mtx_tumor_log2@assays@data$segment_ratios))[mtx_srt$sample,]

#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("subclones", "condition")) 
rownames(anno_mtx) <- NULL

clst_col <- new_pal[1:length(unique(anno_mtx$subclones))]
names(clst_col) <- paste0("c", 1:length(unique(anno_mtx$subclones)))

tp_col <- c("#C27739", "#3A85C2")
names(tp_col) <- c("act_mda231", "arc_well_mda231")
ha_row=rowAnnotation(df = anno_mtx, col = list(subclones= clst_col, condition=tp_col), show_annotation_name = F)

#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")
breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue4", "white", "firebrick4"))

pdf(paste0("./figures/", pro_name, "_complexHeatmap3.pdf"), height = 8, width = 8)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        row_split = anno_mtx$subclones, name = "scheatmap", show_row_names = F, show_column_names = F, 
        row_title = paste0("single cells: ", nrow(ht_mtx), " (",names(table(anno_mtx$condition)[1]), ": ", 
                           table(anno_mtx$condition)[1], "; ",
                           names(table(anno_mtx$condition)[2]), ": ", table(anno_mtx$condition)[2], ")"),
        column_title = paste0(pro_name, "_scHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                   labels_gp = gpar(fontsize = 12)), top_annotation = ha_col, left_annotation = ha_row)
dev.off()


#------Aneuploid QC-----#####
pro_name <- c("arcwell_act")
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
my_mtx <- varbin_mtx_tumor_log2@colData %>% as.data.frame()
my_mtx %>% dplyr::filter(condition == "arc_well_mda231") %>% pull(reads_total) %>% mean()
my_mtx %>% dplyr::filter(condition == "arc_well_mda231") %>% pull(reads_assigned_bins) %>% mean()
my_mtx %>% dplyr::filter(condition == "arc_well_mda231") %>% pull(percentage_duplicates) %>% mean()
my_mtx %>% dplyr::filter(condition == "arc_well_mda231") %>% pull(median_bin_count) %>% mean()

#-----Arc-well and ACT correlation-----#######
seg_all <- log2(varbin_mtx_tumor_log2@assays@data$segment_ratios)
seg_arc <- apply(seg_all[,str_detect(colnames(seg_all), "mda231_9x")], 1, median)
seg_act <- apply(seg_all[,!str_detect(colnames(seg_all), "mda231_9x")], 1, median)
cor(seg_arc, seg_act)

