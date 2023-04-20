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
tp_col <- c("#EA3291", "#96C942")

setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload")
pro_name <- c("nki23_p5")
raw_path <- paste0("./map_seg_output/", pro_name)

keep_gene <- c("SHC1","EGFR", "FGFR1","MYC", "CDKN2A", "PTEN", "CDK4", "MDM2", "CCNE1", "AURKA")
keep_gene2 <- c("AKT3","FGFR4","EGFR","FGFR1", "MYC", "CDK4", "CCNE1", "AURKA")

#---------processing------####
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx # 2372
filter_cells = read.table(paste0("./metrics/", pro_name, "_filtered_bincounts_newnormal.txt"), header = T)
filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] %>% janitor::make_clean_names()
length(filt_cells_names) # 1966

varbin_mtx_tumor <- varbin_mtx[,filt_cells_names]
varbin_mtx_tumor@colData

#---add meta info--
name_meta <- rownames(varbin_mtx_tumor@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%  
  mutate(dispense = stringr::str_extract(my_name, "yes|no"), 
         peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_nki23p_"), "primary", 
                            ifelse(stringr::str_detect(my_name, "_nki23re_"), "recurrence", "notimepoint"))) %>% 
  dplyr::select(c("dispense", "peak", "timepoint"))

varbin_mtx_tumor@colData <- cbind(varbin_mtx_tumor@colData, name_meta)
varbin_mtx_tumor_log <- logNorm(varbin_mtx_tumor, transform = "log2")
saveRDS(varbin_mtx_tumor_log, file = paste0("./objects/", pro_name, c("_filtered_copykit.rds")))
# varbin_mtx_tumor_log <- readRDS(paste0("./objects/", pro_name, c("_filtered_copykit.rds")))
#----UMAP-----
set.seed(31)
near_nb <- 20
umap_data <-  data.frame(uwot::umap(log2(t(varbin_mtx_tumor_log@assays@data$segment_ratios)), metric = "manhattan", 
                                    min_dist = 0.1, spread = 3, n_neighbors = near_nb))

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
cowplot::ggsave2(paste0("./figures/", pro_name, "_passDarlanQC_cells_umap.pdf"), p1, width = 5, height = 4)

#----remove outlier cluster c0 and cluster that have less 6 cells.
#----paired samples
clone_num <- table(varbin_mtx_tumor_log@colData$subclones, varbin_mtx_tumor_log@colData$timepoint)
clone_num_less6 <- as.data.frame(clone_num) %>% filter(Var1 != "c0") %>% 
  filter(Freq > 0 & Freq <6) %>% mutate(comb = paste(Var2, Var1, sep = "_"))

varbin_mtx_tumor_log2 <- varbin_mtx_tumor_log[, !(subclones == "c0" | 
                                                    (varbin_mtx_tumor_log@colData$tp_clst %in% clone_num_less6$comb))]
varbin_mtx_tumor_log2@colData

p1 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData),aes(x = UMAP_1, y = UMAP_2, fill = subclones)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = new_pal) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p1
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap.pdf"), p1, width = 5, height = 4)

tp_col <- c("deeppink", "chartreuse1")
p2 <- ggplot(as.data.frame(varbin_mtx_tumor_log2@colData),aes(x = UMAP_1, y = UMAP_2, fill = timepoint)) + 
  geom_point(shape = 21, size=2.5, stroke = 0.03) + 
  scale_fill_manual(values = tp_col) + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank())
p2
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap2.pdf"), p2, width = 5, height = 4)

p3 <- plot_grid(plotlist=list(p1,p2), ncol=2, align='h')
p3
cowplot::ggsave2(paste0("./figures/", pro_name, "_final_filtered_cells_umap_merge.pdf"), p3, width = 7.8, height = 2.8)

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))
#------heatmap with clonal status bar----#####
mtx_srt <- as.data.frame(varbin_mtx_tumor_log2@colData) %>% 
  arrange(factor(subclones, levels = paste0("c", 1:length(unique(varbin_mtx_tumor_log2@colData$subclones)))))
ht_mtx <- log2(t(varbin_mtx_tumor_log2@assays@data$segment_ratios))[mtx_srt$sample,]

#----annotation bar--
anno_mtx <- mtx_srt %>% dplyr::select(c("peak", "subclones", "timepoint")) 
rownames(anno_mtx) <- NULL
peak_col <- c("#219ebc","#f4a261")
names(peak_col) <- c("d", "a")
clst_col <- new_pal[1:length(unique(anno_mtx$subclones))]
names(clst_col) <- paste0("c", 1:length(unique(anno_mtx$subclones)))
tp_col <- c("deeppink", "chartreuse1")
names(tp_col) <- c("primary", "recurrence")
#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue4", "white", "firebrick4"))

#----clonal status---
low_cutoff <- -0.15
up_cutoff <-  0.15
cell_cutoff <- 0.95
neu_cutoff <- 0.90
my_sample <- clonality_log_trinary_neu(log_ratio_df=ht_mtx, lower_cutoff = low_cutoff, upper_cutoff = up_cutoff, 
                                       cell_pct = cell_cutoff, neu_pct = neu_cutoff)
table(my_sample)
gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% keep_gene)

gene_bin_sel_sample <- as.data.frame(my_sample) %>% rownames_to_column() %>% dplyr::mutate(rowname = as.numeric(rowname)) %>% 
  dplyr::right_join(gene_bin_sel, by = c("rowname" = "pos")) %>% 
  dplyr::mutate(my_color = plyr::mapvalues(my_sample, from = c("sCNA", "cCNA","neu"), to = c("grey70","purple2","grey92"))) 

my_col = structure(c("grey70","grey92","purple2"), names = c("sCNA","neu","cCNA"))
my_sample_df <- as.data.frame(t(my_sample))

cs <- as.data.frame(my_sample)
ha_bottom = columnAnnotation(df = cs, col = list(my_sample=c("sCNA"="grey70", "neu" = "grey92","cCNA"="purple2")), 
                             clonal_state = anno_mark(at=gene_bin_sel_sample$rowname, labels = gene_bin_sel_sample$gene, 
                                           side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel_sample$my_color)))

ha_row2=rowAnnotation(df = anno_mtx[,3:2], col = list(subclones= clst_col, timepoint=tp_col), 
                      show_annotation_name = F, simple_anno_size = unit(0.7, "cm"))

pdf(paste0("./figures/", pro_name, "_complexHeatmap_with_clonal_state.pdf"), height = 10, width = 10)
Heatmap(as.matrix(ht_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        row_split = anno_mtx$subclones,
        name = "scheatmap", show_row_names = F, show_column_names = F, row_title = paste0(nrow(ht_mtx), " single cells"),
        column_title = paste0(pro_name, "_low_",low_cutoff, "_up_", up_cutoff, "_cell_pct_", cell_cutoff,"_",neu_cutoff), 
        use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", 
                                    title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row2, bottom_annotation = ha_bottom)
dev.off()

#----calculate overdispersion---######
bin_count <- varbin_mtx_tumor_log2@assays@data$bin_counts
bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% dplyr::rename(over_disp = V1) %>% dplyr::select("over_disp")
varbin_mtx_tumor_log2@colData <- cbind(varbin_mtx_tumor_log2@colData, bin_count_overdisp2) 

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_overdisp_copykit.rds")))

#-----event matrices and integer copy number----####
ploidy_pri <- 3.581553487
ploidy_rec <- 5.23689588
varbin_mtx_tumor_log2@colData$ploidy_used <- as.numeric(plyr::mapvalues(varbin_mtx_tumor_log2@colData$timepoint, 
                                                                        from = c("primary","recurrence"), 
                                                                        to = c(ploidy_pri, ploidy_rec)))

varbin_mtx_tumor_log2 <- calcInteger(varbin_mtx_tumor_log2, assay = "segment_ratios", method = "fixed", 
                                     ploidy_value = colData(varbin_mtx_tumor_log2)$ploidy_used)
varbin_mtx_tumor_log2 <- calcConsensus(varbin_mtx_tumor_log2, assay = "integer", fun="median", consensus_by = "subclones")

ploidy_trunc <- round(2*mean(colData(varbin_mtx_tumor_log2)$ploidy_used))
eventmat <- getEventMat(varbin_mtx_tumor_log2, bin_adj = 2, ploidy_trunc = ploidy_trunc)

## convert back to segmentation level
popseg_long <- as.data.frame(apply(as.data.frame(t(eventmat %>% dplyr::select(matches("c[0-9]+")))), 1, 
                                   function(m) {rep.int(m, eventmat$n.bins)}))

## make medicc input matrix
medicc_input <- cbind(SummarizedExperiment::rowRanges(varbin_mtx_tumor_log2) %>% dplyr::as_tibble() %>% 
                        dplyr::select(seqnames, start, end), popseg_long) %>%
  mutate(diploid=2) %>%          ## here we manually add a diploid cell CN profile to be used as root in the medicc tree
  dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% 
  dplyr::select(sample_id,chrom, everything())

dir.create(paste0("./metrics/medicc_files/", pro_name), recursive = T)
write_tsv(medicc_input, file = paste0("./metrics/medicc_files/", pro_name, "/",pro_name, c("_medicc2_input.tsv")))
# medicc_input <- read_tsv(file = paste0("./metrics/medicc_files/", pro_name, "/",pro_name, c("_medicc2_input.tsv")))

attr(popseg_long, "consensus_by") <- "subclones"
attr(popseg_long, "consensus_assay") <- "integer"
copykit::consensus(varbin_mtx_tumor_log2) <- popseg_long
varbin_mtx_tumor_log2 <- runConsensusPhylo(varbin_mtx_tumor_log2)

saveRDS(varbin_mtx_tumor_log2, file = paste0("./objects/", pro_name, c("_final_filtered_overdisp_integerCN_copykit.rds")))
# varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_overdisp_integerCN_copykit.rds")))

#-----MEDICC2 Tree----######
mywd <- getwd()
#----go to terminal and run medicc2---
# conda activate medicc_env
# medicc2 -a CN --total-copy-numbers -j 40 -vv input_path output_path

tsv_out_path <- paste0(mywd, "/metrics/medicc_files/", pro_name, "/")
medic <- ape::read.tree(paste0(tsv_out_path, pro_name,"_medicc2_input_final_tree.new"))
tree_events <- read_tsv(paste0(tsv_out_path, pro_name,"_medicc2_input_copynumber_events_df.tsv"))
tree_profile <- read_tsv(paste0(tsv_out_path, pro_name,"_medicc2_input_final_cn_profiles.tsv"))
ploidy_trunc <- round(2*mean(colData(varbin_mtx_tumor_log2)$ploidy_used))

meta_con <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% 
  dplyr::select(subclones, timepoint) %>% group_by(subclones, timepoint) %>%
  summarise(n=n()) %>% spread(timepoint, n) 
meta_con[is.na(meta_con)] <- 0
meta_con_p <- meta_con %>% mutate(sum=primary+recurrence) %>% mutate(primary=primary/sum, recurrence=recurrence/sum) %>% 
  dplyr::select(-sum) %>% as.data.frame()

list_samples <- split(meta_con$subclones, meta_con$subclones)
tree <- ggtree::groupOTU(medic, list_samples)
treeplt <- ggtree::ggtree(ape::ladderize(tree), ladderize = FALSE, size = .2) +
  ggtree::geom_tiplab(size=6, aes(color=group),hjust = -0.6, alpha=1)+
  scale_colour_manual(values = new_pal2,breaks = names(new_pal2)) +
  geom_text(aes(x=branch, label=round(branch.length, 0), vjust=-.5), size = 3) +
  theme(legend.position = "none") + ggtree::geom_rootpoint() 

treeplt_nodis <- ggtree::ggtree(ape::ladderize(tree), ladderize = FALSE, size = .2) +
  ggtree::geom_tiplab(size=6, aes(color=group),hjust = -0.6, alpha=1)+
  scale_colour_manual(values = new_pal2,breaks = names(new_pal2)) +
  # geom_text(aes(x=branch, label=round(branch.length, 0), vjust=-.5), size = 3) +
  theme(legend.position = "none") + ggtree::geom_rootpoint() 
#print(treeplt)

pie_mat_m <- meta_con
pies = list()
for (i in 1:length(tree$tip.label)) {
  if(tree$tip.label[i] %in% pie_mat_m$subclones){
    curr_dat = reshape2::melt(pie_mat_m[pie_mat_m$subclones==tree$tip.label[i],]) %>% dplyr::filter(value!=0)
    ## create a ggplot object for each pie chart
    pies[[i]] =  ggplot(curr_dat, aes(y = value, fill = variable, x="")) +
      geom_bar(stat = "identity") +
      coord_polar("y", start=0) +
      theme_void() + scale_fill_manual(values = c("primary" = "#ED2A91", "recurrence" = "#94C93D"), guide = F)
  }else{
    pies[[i]] = ggplot(data.frame(value=100), aes(y = value, x=""), fill="gray")+
      geom_bar(stat = "identity") +
      coord_polar("y", start=0) +
      theme_void() 
  }
}

names(pies) = 1:length(tree$tip.label)
treeplt2 <- ggtree::inset(treeplt, pies, width=0.1, height=0.1)
treeplt2
cowplot::ggsave2(paste0("./figures/", pro_name, "_medicc2_tree.pdf"), treeplt2, width = 4.5, height = 2.25)

lab <-  (treeplt[["data"]] %>% arrange(y))$label
my_order <- rev(lab[grepl("c[0-9]+", lab, perl = T)])

#----plot consensus heatmap ----
cs_mtx <- t(varbin_mtx_tumor_log2@consensus)
cs_mtx_order <- cs_mtx[my_order,]

clst_col_cs <- new_pal[1:length(unique(varbin_mtx_tumor_log2@colData$subclones))]
subclone_cs <- paste0("c", 1:length(unique(varbin_mtx_tumor_log2@colData$subclones)))
names(clst_col_cs) <- subclone_cs
ha_row_cs=rowAnnotation(subclones = my_order, col = list(subclones = clst_col_cs), show_annotation_name = F, show_legend = F)

meta_con_p2 <- meta_con_p %>% column_to_rownames(var = c("subclones"))
ha_barplot <- rowAnnotation(bar= anno_barplot(meta_con_p2[my_order,] %>% dplyr::select(primary, recurrence),
                                              gp = gpar(fill = c("primary" = "#ED2A91", "recurrence" = "#94C93D"), 
                                                        col= c("primary" = "#ED2A91", "recurrence" = "#94C93D"))),
                            show_annotation_name = FALSE)

#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

#-----annotated genes---
gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% keep_gene2)
ha_bottom_cs = columnAnnotation(clonal_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                         side = "bottom", labels_gp = gpar(fontsize = 14)))

col_vec_cs <- structure(pals::ocean.balance(length(0:ploidy_trunc)), names = 0:ploidy_trunc)

pdf(paste0("./figures/", pro_name, "_complexHeatmap_consensus_integer_with_barplot.pdf"), height = 4, width = 10)
Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = paste0(pro_name, "_consensusHeatmap"),
        use_raster = T, raster_quality = 5, col = col_vec_cs, 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)), 
        top_annotation = ha_col, left_annotation = ha_row_cs, bottom_annotation = ha_bottom_cs, right_annotation = ha_barplot)
dev.off()

#-----PCA and RCA-----########
#--check PCA and RCA node---
treeplt_temp <- ggtree::ggtree(ape::ladderize(tree), ladderize = FALSE, size = .2) + ggtree::geom_tiplab()+ggtree::geom_nodelab()
ggtree::inset(treeplt_temp, pies, width=0.07, height=0.07)

pca <- "internal_4"
rca <- "internal_8"

mrca_mat <- tree_profile %>% dplyr::filter(sample_id %in% c(pca,rca)) %>%
  dplyr::select(-starts_with("is_")) %>% spread(sample_id, CN) %>%
  dplyr::select(-chrom, -start, -end) %>% dplyr::select(all_of(pca),all_of(rca))

colnames(mrca_mat) <- c("PCA", "RCA")
gc_mrca <- consensus_genomic_classes(consensus_int = t(mrca_mat))

saveRDS(gc_mrca, file = paste0("./objects/cs_bin_class/", pro_name, c("_cs_bin_class.rds")))
# gc_mrca <- readRDS(file = paste0("./objects/cs_bin_class/", pro_name, c("_cs_bin_class.rds")))

gc_mrca_df <- gc_mrca %>% as.data.frame()
colnames(gc_mrca_df) <- "clonal_states"
gene_bin_sel$col <- plyr::mapvalues(gc_mrca[gene_bin_sel$pos], from = c("sCNA","cCNA"), to = c("#F47F20", "#424451"))
gene_bin_sel2 <- gene_bin_sel %>% dplyr::filter(col == "#F47F20")

ha_bottom_pc = columnAnnotation(df = gc_mrca_df, col = list(clonal_states=c("sCNA"="#F47F20", "cCNA"="#424451")), 
                                clonal_state = anno_mark(at=gene_bin_sel2$pos, labels = gene_bin_sel2$gene, 
                                                         side = "bottom", labels_gp = gpar(fontsize = 14, col = gene_bin_sel2$col)))

col_vec_cs <- structure(pals::ocean.balance(length(0:ploidy_trunc)), names = 0:ploidy_trunc)

pc <- Heatmap(as.matrix(cs_mtx_order), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
              name = "csheatmap", show_row_names = T, row_names_side = "left",
              show_column_names = F, column_title = paste0(pro_name, "_consensusHeatmap"),
              use_raster = T, raster_quality = 5, col = col_vec_cs, 
              heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                          labels_gp = gpar(fontsize = 12)), 
              top_annotation = ha_col, left_annotation = ha_row_cs, right_annotation = ha_barplot)

ppr <- Heatmap(t(mrca_mat), cluster_columns = FALSE, border = TRUE, cluster_rows = FALSE, show_row_dend = FALSE, 
               name = "pcarcaheatmap", show_row_names = T, row_names_side = "left",
               show_column_names = F, column_title = paste0(pro_name, "_consensusHeatmap"),
               use_raster = T, raster_quality = 5, col = col_vec_cs, show_heatmap_legend = F,
               bottom_annotation = ha_bottom_pc)

pdf(paste0("./figures/", pro_name, "_complexHeatmap_consensus_integer_pca_rca_merged.pdf"), height = 5.5, width = 9)
draw(pc %v% ppr, ht_gap = unit(0.5, 'cm'), merge_legend = T,auto_adjust = T)
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

meta_mtx %>% dplyr::filter(timepoint == "primary") %>% pull(reads_total) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "primary") %>% pull(reads_assigned_bins) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "primary") %>% pull(percentage_duplicates) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "primary") %>% pull(median_bin_count) %>% mean()

meta_mtx %>% dplyr::filter(timepoint == "recurrence") %>% pull(reads_total) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "recurrence") %>% pull(reads_assigned_bins) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "recurrence") %>% pull(percentage_duplicates) %>% mean()
meta_mtx %>% dplyr::filter(timepoint == "recurrence") %>% pull(median_bin_count) %>% mean()

#------filtering status heatmap---####
varbin_mtx <- readVarbinCNA(raw_path, remove_Y = TRUE)
varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name, c("_final_filtered_copykit.rds")))

all_mtx <- as.data.frame(varbin_mtx@colData) %>% 
  mutate(filter_state = ifelse(sample %in% rownames(varbin_mtx_tumor_log2@colData),"kept", "removed")) %>% 
  dplyr::select("filter_state")
all_name_meta <- rownames(varbin_mtx@colData) %>% as.data.frame() %>% dplyr::rename(my_name = ".") %>%  
  mutate(peak = ifelse(stringr::str_detect(my_name, "_d_"), "d", ifelse(stringr::str_detect(my_name, "_a_"), "a", "nopeakinfo")),
         timepoint = ifelse(stringr::str_detect(my_name, "_nki23p_"), "primary", 
                            ifelse(stringr::str_detect(my_name, "_nki23re_"), "recurrence", "notimepoint"))) %>% 
  dplyr::select(c("peak", "timepoint"))

varbin_mtx@colData <- cbind(varbin_mtx@colData, all_mtx, all_name_meta)
all_mtx_srt <- as.data.frame(varbin_mtx@colData) %>% arrange(factor(filter_state, levels = c("kept","removed")))
ht_all_mtx <- log2(t(varbin_mtx@assays@data$segment_ratios))[all_mtx_srt$sample,]

#----annotation bar--
anno_mtx <- all_mtx_srt %>% dplyr::select(c("filter_state","timepoint","peak")) 
rownames(anno_mtx) <- NULL
peak_col <- c("#219ebc","#f4a261","#d9d9d9")
names(peak_col) <- c("d", "a","nopeakinfo")
fs_col <- c("#2E8B58", "#BFBEBE")
names(fs_col) <- c("kept", "removed")
tp_col <- c("deeppink", "chartreuse1")
names(tp_col) <- c("primary", "recurrence")
ha_row=rowAnnotation(df = anno_mtx, col = list(peak=peak_col, filter_state=fs_col, timepoint=tp_col), show_annotation_name = F)

#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")

breaks = c(-1,0,1)
col_vec = circlize::colorRamp2(breaks =breaks, c("dodgerblue4", "white", "firebrick4"))
# col_vec = circlize::colorRamp2(breaks =breaks, c("blue", "white", "red"))

pdf(paste0("./figures/", pro_name, "_complexHeatmap_filter_states.pdf"), height = 8, width = 6)
Heatmap(as.matrix(ht_all_mtx), cluster_columns = FALSE, border = TRUE, cluster_rows = TRUE, show_row_dend = FALSE, 
        row_split = anno_mtx$filter_state,
        name = "scheatmap", show_row_names = F, show_column_names = F, 
        row_title = paste0("single cells: total: ",nrow(anno_mtx)," (",names(table(anno_mtx$filter_state)[1]), ": ", 
                           table(anno_mtx$filter_state)[1],"; ", names(table(anno_mtx$filter_state)[2]), ": ", 
                           table(anno_mtx$filter_state)[2], ")"), 
        column_title = paste0(pro_name, "_scHeatmap"), use_raster = T, raster_quality = 5, col = col_vec, 
        heatmap_legend_param = list(title = "Log2 (Ratio)", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)), top_annotation = ha_col, left_annotation = ha_row)
dev.off()

