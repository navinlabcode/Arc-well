#####-----------Process-------######
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
library(copykit)
library(tidyverse)
library(dbscan)
library(dendextend)
library(ComplexHeatmap)
library(ggplot2)
require(RColorBrewer)
library(ape)
library(ggtree)
library(Homo.sapiens)
options(max.print = 200)

source("./scripts/0_arcwell_functions.R")
load("./pre_load_data/pre_load_data.rda")
chr_name <- unlist(read.table("./pre_load_data/chr_name.txt"), use.names = F)
chr_color <- read.table("./pre_load_data/chr_color.txt")
gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", 
            "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")
new_pal2 <- new_pal
names(new_pal2) <- paste0("c", 1:length(new_pal2))
tp_col <- c("deeppink", "chartreuse1")

setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload")
pro_name <- c("duke249_p1")
raw_path <- paste0("./map_seg_output/", pro_name)

keep_gene <- c("SHC1","ESR1","EGFR", "MYC", "CDKN2A","GATA3", "MDM2", "BRCA2", "RB1", "TP53", "ERBB2", "BRCA1","BCL2", "AURKA")

#---------processing------####
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx # 1082
filter_cells = read.table(paste0("./metrics/", pro_name, "_filtered_bincounts_newnormal.txt"), header = T)
filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] %>% janitor::make_clean_names()
length(filt_cells_names) # 897

varbin_mtx_tumor <- varbin_mtx[,filt_cells_names]
varbin_mtx_tumor@colData
#---add meta info--
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%  
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_duke_26_"), "primary", 
                            ifelse(stringr::str_detect(my_name, "_duke_27_"), "recurrence", "notimepoint"))) %>% 
  dplyr::select(c("dispense", "peak", "timepoint"))

varbin_mtx_tumor@colData <- cbind(varbin_mtx_tumor@colData, name_meta)
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

#---find clusters----
set.seed(17)
hdb_data <- dbscan::hdbscan(umap_data2, minPts = nrow(umap_data2)*0.015)
subclones <- paste0("c", as.character(hdb_data$cluster))

varbin_mtx_tumor_log@colData <- cbind(varbin_mtx_tumor_log@colData, subclones)
varbin_mtx_tumor_log@colData$subclones <- factor(varbin_mtx_tumor_log@colData$subclones, 
                                                levels = paste0("c", 0:(length(table(subclones))-1)))
varbin_mtx_tumor_log@colData$tp_clst <- paste(varbin_mtx_tumor_log@colData$timepoint, 
                                              varbin_mtx_tumor_log@colData$subclones, sep="_")

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = c("grey",new_pal)) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name, "_passInitialQC_cells_umap.pdf"), p1, width = 5, height = 4)

#----remove outlier cluster c0 and cluster that have less 3 cells.
#----paired samples
clone_num <- table(varbin_mtx_tumor_log@colData$subclones, varbin_mtx_tumor_log@colData$timepoint)
clone_num_less3 <- as.data.frame(clone_num) %>% filter(Var1 != "c0") %>% 
  filter(Freq > 0 & Freq <4) %>% mutate(comb = paste(Var2, Var1, sep = "_"))

varbin_mtx_tumor_log2 <- varbin_mtx_tumor_log[, !(subclones == "c0" | (varbin_mtx_tumor_log@colData$tp_clst %in% clone_num_less3$comb))]
varbin_mtx_tumor_log2@colData

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = new_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap.pdf"), p1, width = 5, height = 4)

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
#----plot heatmap ----
#----scheatmap and concensus heatmap were plotted separately, since ComplexHeatmap will treat one row as one cell, 
#--so consensus heatmap will become very narrow if plot together.
mtx_srt <- as.data.frame(varbin_mtx_tumor_log2@colData) %>% 
  arrange(factor(subclones, levels = paste0("c", 1:length(unique(varbin_mtx_tumor_log2@colData$subclones)))))
ht_mtx <- log2(t(varbin_mtx_tumor_log2@assays@data$segment_ratios))[mtx_srt$sample,]

#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("peak", "subclones")) 
rownames(anno_mtx) <- NULL
peak_col <- c("#219ebc","#f4a261")
names(peak_col) <- c("d", "a")
clst_col <- new_pal[1:length(unique(anno_mtx$subclones))]
names(clst_col) <- paste0("c", 1:length(unique(anno_mtx$subclones)))
ha_row=rowAnnotation(df = anno_mtx, col = list(peak=peak_col, subclones= clst_col), show_annotation_name = F)
#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue4", "white", "firebrick4"))

pdf(paste0("./figures/", pro_name, "_complexHeatmap.pdf"), height = 8, width = 8)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        row_split = anno_mtx$subclones,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0("single cells:", nrow(ht_mtx)),
        column_title = paste0(pro_name, "_scHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)), top_annotation = ha_col, left_annotation = ha_row)
dev.off()

#----calculate overdispersion---######
bin_count <- varbin_mtx_tumor_log2@assays@data$bin_counts
bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% dplyr::rename(over_disp = V1) %>% dplyr::select("over_disp")
varbin_mtx_tumor_log2@colData <- cbind(varbin_mtx_tumor_log2@colData, bin_count_overdisp2) 

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))

#-----event matrices and integer copy number----####
ploidy_pri <- 2.804038169
varbin_mtx_tumor_log2@colData$ploidy_used <- ploidy_pri

varbin_mtx_tumor_log2 <- calcInteger(varbin_mtx_tumor_log2, assay = "segment_ratios", method = "fixed", 
                                     ploidy_value = colData(varbin_mtx_tumor_log2)$ploidy_used)
varbin_mtx_tumor_log2 <- calcConsensus(varbin_mtx_tumor_log2, assay = "integer", fun="median", consensus_by = "subclones")

ploidy_trunc <- round(2*mean(colData(varbin_mtx_tumor_log2)$ploidy_used))
eventmat <- getEventMat(varbin_mtx_tumor_log2, bin_adj = 2, ploidy_trunc = ploidy_trunc)

## convert back to segmentation level
popseg_long <- as.data.frame(apply(as.data.frame(t(eventmat %>% dplyr::select(matches("c[0-9]+")))), 1, 
                                   function(m) {rep.int(m, eventmat$n.bins)}))

attr(popseg_long, "consensus_by") <- "subclones"
attr(popseg_long, "consensus_assay") <- "integer"
copykit::consensus(varbin_mtx_tumor_log2) <- popseg_long
varbin_mtx_tumor_log2 <- runConsensusPhylo(varbin_mtx_tumor_log2)

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name, c("_final_filtered_overdisp_integerCN_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_overdisp_integerCN_copykit.rds")))

#----plot consensus heatmap ----
cs_mtx <- t(varbin_mtx_tumor_log2@consensus)
clst_col_cs <- new_pal[1:length(unique(varbin_mtx_tumor_log2@colData$subclones))]
subclone_cs <- paste0("c", 1:length(unique(varbin_mtx_tumor_log2@colData$subclones)))
names(clst_col_cs) <- subclone_cs
ha_row_cs=rowAnnotation(subclones = subclone_cs, col = list(subclones = clst_col_cs), show_annotation_name = F)

#----header---
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")
#-----annotated genes---
gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% keep_gene)
ha_bottom_cs = columnAnnotation(clonal_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                      side = "bottom", labels_gp = gpar(fontsize = 14)))

col_vec_cs <- structure(pals::ocean.balance(length(0:ploidy_trunc)), names = 0:ploidy_trunc)

pdf(paste0("./figures/", pro_name, "_complexHeatmap_consensus_integer.pdf"), height = 3, width = 8)
Heatmap(as.matrix(cs_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = F, show_column_names = F, column_title = paste0(pro_name, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec_cs, 
        heatmap_legend_param = list(title = "copy number", 
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs)
dev.off()

#----Number of normal cells-----
filter_mtx <- read.table(paste0("./metrics/", pro_name, "_metadata.metrics_newnormal.txt"), header = T) 
rownames(filter_mtx) <- filter_mtx$sample %>% janitor::make_clean_names() 
filter_mtx2 <- filter_mtx %>% mutate(dups_percentage = dups_removed/total_reads) %>% rownames_to_column() %>% 
  dplyr::select(c("rowname", "timepoint","reads_kept","median_bin_count","filtered","is_normal"))

filter_mtx2 %>% dplyr::filter(filtered == "kept" & is_normal == T) %>% pull(timepoint) %>% table()

#------Aneuploid QC-----#####
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))
meta_mtx <- varbin_mtx_tumor_log2@colData %>% as.data.frame()
table(meta_mtx$timepoint)

meta_mtx %>% dplyr::filter(timepoint == "recurrence") %>% pull(reads_total) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "recurrence") %>% pull(reads_assigned_bins) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "recurrence") %>% pull(percentage_duplicates) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "recurrence") %>% pull(median_bin_count) %>% mean()


#------filtering status heatmap---####
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))

all_mtx <- as.data.frame(varbin_mtx@colData) %>% 
  mutate(filter_state = ifelse(sample %in% rownames(varbin_mtx_tumor_log2@colData),"kept", "removed")) %>% dplyr::select("filter_state")
all_name_meta <- rownames(varbin_mtx@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%  
  mutate(peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_duke_26_"), "primary", 
                            ifelse(stringr::str_detect(my_name, "_duke_27_"), "recurrence", "notimepoint"))) %>% 
  dplyr::select(c("peak", "timepoint"))

varbin_mtx@colData <- cbind(varbin_mtx@colData, all_mtx, all_name_meta)

all_mtx_srt <- as.data.frame(varbin_mtx@colData) %>% arrange(factor(filter_state, levels = c("kept","removed")))
head(all_mtx_srt)
ht_all_mtx <- log2(t(varbin_mtx@assays@data$segment_ratios))[all_mtx_srt$sample,]
topleft(ht_all_mtx)

#----annotation bar--
anno_mtx <- all_mtx_srt %>% dplyr::select(c("filter_state","peak")) 
rownames(anno_mtx) <- NULL
peak_col <- c("#219ebc","#f4a261","#d9d9d9")
names(peak_col) <- c("d", "a","nopeakinfo")
fs_col <- c("#2E8B58", "#BFBEBE")
names(fs_col) <- c("kept", "removed")
ha_row=rowAnnotation(df = anno_mtx, col = list(peak=peak_col, filter_state=fs_col), show_annotation_name = F)
#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue4", "white", "firebrick4"))

pdf(paste0("./figures/", pro_name, "_complexHeatmap_filter_states.pdf"), height = 8, width = 8)
Heatmap(as.matrix(ht_all_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = TRUE, show_row_dend = FALSE, 
        row_split = anno_mtx$filter_state,
        name = "scheatmap", show_row_names = F, show_column_names = F, 
        row_title = paste0("single cells: total: ",nrow(anno_mtx)," removed: ", table(anno_mtx$filter_state)[2], " kept: ", 
                           table(anno_mtx$filter_state)[1]),
        column_title = paste0(pro_name, "_scHeatmap"), use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)), top_annotation = ha_col, left_annotation = ha_row)
dev.off()




