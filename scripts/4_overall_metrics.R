#####-----------Process-------######
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)
require(RColorBrewer)
library(Homo.sapiens)
library(ggpubr)
library(DEGreport)
library(ggplot2)
library(dplyr)
library(ggalt)

options(max.print = 200)

source("./scripts/arcwell_functions.R")
load("./pre_load_data/pre_load_data.rda")
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", 
            "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")
new_pal2 <- new_pal
names(new_pal2) <- paste0("c", 1:length(new_pal2))
load("./pre_load_data/pre_load_data.rda")
chr_name <- unlist(read.table("./pre_load_data/chr_name.txt"), use.names = F)
chr_color <- read.table("./pre_load_data/chr_color.txt")
gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)

tp_col <- c("#EA3291", "#96C942")
sample_name <- c("duke248_p3","duke254_p4","nki23_p5","nki26_p6","nki28_p7","nki19_p8","nki22_p9","nki15_p10","nki31_p11","nki12_p12")

setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload")

#-----clinical features--(Figure 1a)---######
meta_mtx <- read.csv("./metrics/clinical_metadata_table.csv")
meta_mtx2 <- meta_mtx %>% column_to_rownames(var = "sample")
meta_mtx2t <- t(meta_mtx2[,c("ER","PR","HER2","histology","grade","timepoint")])
my_col = structure(c("#EA3291","#96C942","#fdb863","#BBBBBB","#F1F2F2","#B2ABD2","#AD4599","#672872","#d1e5f0","#4393c3"),
                   names = c("Primary","Recurrence","pos","neg","Unk","DCIS","IBC_DCIS","IBC","2","3"))

pdf("./figures/clinical_meta_data.pdf", width = 7.5, height = 3.5)
Heatmap(meta_mtx2t, name = "clinic_meta", col = my_col, row_names_side = c("left"), rect_gp = gpar(col="black", lwd=1),
        column_split = c(letters[1:2], rep(letters[3:12], each =2)))
dev.off()

#-----correlation with block ages--(Figure S4a)---######
all_meta <- read.table(file = "./metrics/merged_QC_meta_all_DCIS_FFPE.txt", sep = "\t", header = T)
din_age_meta <- read.csv("./metrics/din_block_ages.csv", row.names = 1)
timepoint_col <- c("primary" = "#EA3291", "recurrence" = "#96C942")
pro_name2 <- c("duke248_p3","duke254_p4","nki23_p5","nki26_p6","nki28_p7","nki19_p8","nki22_p9","nki15_p10","nki31_p11","nki12_p12")
time_point <- c("primary","recurrence")
sample_levels <- c("duke249_p1_recurrence","nki17_p2_primary",paste(rep(pro_name2, each = length(time_point)), time_point, sep = "_"))

all_meta2 <- all_meta %>% group_by(patient) %>% summarise(mean_overdis=mean(over_disp), mean_pcr_dup = mean(percentage_duplicates))

all_meta3 <- all_meta2 %>% as.data.frame() %>% left_join(din_age_meta, by = c("patient" = "patient")) %>% 
  separate(patient, c("sample", "patient","time_point")) %>% mutate(sample = paste0(sample, "_",patient))

p1 <- ggplot(all_meta3, aes(x=block_age_years, y=mean_overdis)) + geom_point(shape=21, aes(fill=time_point), size =5)+
  geom_smooth(method=lm, color="black", se = F, size = 1) + geom_cor(method = "pearson") + theme_bw() + 
  scale_fill_manual(values = tp_col) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p2 <- ggplot(all_meta3, aes(x=block_age_years, y=mean_pcr_dup)) + geom_point(shape=21, aes(fill=time_point), size =5)+
  geom_smooth(method=lm, color="black", se = F, size = 1) + geom_cor(method = "pearson") + theme_bw() + 
  scale_fill_manual(values = tp_col) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p3 <- ggplot(all_meta3, aes(x=block_age_years, y=DIN)) + geom_point(shape=21, aes(fill=time_point), size =5)+
  geom_smooth(method=lm, color="black", se = F, size = 1) + geom_cor(method = "pearson") + theme_bw() + 
  scale_fill_manual(values = tp_col)

p4 <- plot_grid(plotlist=list(p1,p2,p3), ncol=1, align='v')
p4
cowplot::ggsave2("./figures/block_age_correlation.pdf", p4, width = 4.5, height = 8)

p5 <- ggplot(all_meta3, aes(x=DIN, y=mean_overdis)) + geom_point(shape=21, aes(fill=time_point), size =5)+
  geom_smooth(method=lm, color="black", se = F, size = 1) + geom_cor(method = "pearson") + theme_bw() + 
  scale_fill_manual(values = c("#EA3291","#96C942")) + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + ggtitle("DIN")
p6 <- ggplot(all_meta3, aes(x=DIN, y=mean_pcr_dup)) + geom_point(shape=21, aes(fill=time_point), size =5)+
  geom_smooth(method=lm, color="black", se = F, size = 1) + geom_cor(method = "pearson") + theme_bw() + 
  scale_fill_manual(values = c("#EA3291","#96C942"))

p7 <- plot_grid(plotlist=list(p5,p6), ncol=1, align='v')
p7
cowplot::ggsave2("./figures/din_correlation.pdf", p7, width = 4.5, height = 6)

#--------Number of subclones and cells (Figure 3a)-------#######
sclone <- read.csv("./metrics/subclone_number.csv", header = T)
all_meta <- read.table(file = "./metrics/merged_QC_meta_all_DCIS_FFPE.txt", sep = "\t", header = T)
cell_num <- as.data.frame(table(all_meta$patient)) %>% mutate(sample = stringr::str_extract(Var1, ".+_p\\d+"), 
                                                              timepoint = stringr::str_extract(Var1, "primary|recurrence")) %>%
  dplyr::filter(!(sample %in% c("duke249_p1", "nki17_p2")))
write.csv(cell_num, "./metrics/cell_number.csv")
# cell_num <- read.csv("./metrics/cell_number.csv", row.names = 1)
cell_num$sample <- factor(cell_num$sample, levels = sample_name)

cell_num %>% filter(timepoint == "primary") %>% pull(Freq) %>% mean()
cell_num %>% filter(timepoint == "primary") %>% pull(Freq) %>% sd()
cell_num %>% filter(timepoint == "recurrence") %>% pull(Freq) %>% mean()
cell_num %>% filter(timepoint == "recurrence") %>% pull(Freq) %>% sd()


all_meta %>% dplyr::select("tp_clst","patient") %>% dplyr::filter(!str_detect(patient, "duke249|nki17")) %>% 
  distinct() %>% pull(patient) %>% table() %>% as.data.frame()
tes <- all_meta %>% dplyr::select("tp_clst","patient") %>% dplyr::filter(!str_detect(patient, "duke249|nki17")) %>% 
  distinct() %>% separate(tp_clst, c("timepoint","subclone")) %>% mutate(patient2 = stringr::str_extract(patient, ".+_p\\d+")) %>% 
  dplyr::select("subclone","patient2") %>% group_by_all() %>% summarise(count = n()) %>% dplyr::filter(count > 1) %>% 
  pull(patient2) %>% table() %>% as.data.frame()
tes
sclone_l <- gather(sclone, samples, subclone_number, sample_name)
sclone_l$timepoints <- factor(sclone_l$timepoints, levels = c("primary","shared","recurrence"))
sclone_l$samples <- factor(sclone_l$samples, levels = sample_name)

my_pale <- c("#EC2E91","#af8dc3","#95C942")
p1 <- ggplot(sclone_l, aes(fill=timepoints, y=subclone_number, x=samples, label= subclone_number)) +
  geom_bar(position = "stack", stat = "identity") + scale_fill_manual(values =  my_pale) + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +theme_classic() + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p2 <- ggplot(cell_num, aes(fill=timepoint, y=Freq, x=sample, label= Freq)) + 
  geom_bar(position = "stack", stat = "identity") + scale_fill_manual(values =  tp_col) + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank())
p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v')

cowplot::ggsave2("./figures/subclone_and_cell_number_in_allsamples.pdf", p3, width = 7, height = 4)

#-----Number of CNA events (Figure 3c)-----#####
cna_num <- data.frame()
for (i in sample_name) {
  print(paste0("Now is runnig: ", i))
  obj <- readRDS(paste0("./objects/", i, c("_final_filtered_overdisp_integerCN_copykit.rds")))
  pri_n = countEvents(obj, "primary", 1)
  rec_n = countEvents(obj, "recurrence", 1)
  
  my_res1 <- cbind(i, pri_n, "Primary")
  my_res2 <- cbind(i, rec_n, "Recurrence")
  my_res <- rbind(my_res1, my_res2)
  cna_num <- rbind(cna_num, my_res)
}  

colnames(cna_num) <- c("sample","num_cna","timepoint")
write.csv(cna_num, "./metrics/number_cna_events.csv")
# cna_num <- read.csv("./metrics/number_cna_events.csv", row.names = 1)

sample_col <- c("#565FAB","#CF3D32","#5DB1DD","#812268","#78C36A","#486A84","#BB6438","#F0E584","#D595A7","#749B59")
cna_num <- cna_num %>% mutate(num_cna = as.numeric(num_cna))%>% mutate(sample = fct_relevel(sample, sample_name))

p1 <- ggplot(cna_num, aes(x=timepoint, y=num_cna)) + geom_line(aes(group=sample), linetype = "dashed", colour = "grey", size = 0.5) + 
  geom_point(aes(color=sample), size=4) + theme_bw() + scale_color_manual(values =  sample_col) + theme(axis.title = element_blank())
p1
cowplot::ggsave2("./figures/CNA_events_number_in_allsamples.pdf", p1, width = 3.75, height = 4)

#---statistic analysis---
cna_num2 <- cna_num[order(cna_num$sample),]
pri_sample <- subset(cna_num2, timepoint == "Primary", num_cna, drop = TRUE)
rec_sample <- subset(cna_num2, timepoint == "Recurrence", num_cna, drop = TRUE)
wilcox.test(pri_sample, rec_sample, paired = TRUE)

#-----single cell MPD (Figure 3d)-----#####
mpd_res <- data.frame()
for (i in sample_name) {
  print(paste0("Now is runnig: ", i))
  obj <- readRDS(paste0("./objects/", i, c("_final_filtered_copykit.rds")))
  pri_cell <- obj@colData %>% as.data.frame() %>% dplyr::filter(timepoint == "primary") %>% pull(sample)
  rec_cell <- obj@colData %>% as.data.frame() %>% dplyr::filter(timepoint == "recurrence") %>% pull(sample)
  
  pri_df <- as_tibble(obj@assays@data$segment_ratios[, pri_cell])
  rec_df <- as_tibble(obj@assays@data$segment_ratios[, rec_cell])
  
  Primary <- mpd_scTree(df = pri_df,n_threads = 10)
  Recurrence <- mpd_scTree(df = rec_df,n_threads = 10)
  
  my_res <- cbind(Primary, Recurrence, i)
  mpd_res <- rbind(mpd_res, my_res)
}  

colnames(mpd_res) <- c("Primary","Recurrence","sample")
write.csv(mpd_res, "./metrics/mpd_scTree.csv")
# mpd_res <- read.csv("./metrics/mpd_scTree.csv", row.names = 1)

mpd2 <- gather(mpd_res, key = "timepoint", value = "values", 1:2)
mpd2 <- mpd2 %>% mutate(values = as.numeric(values)) %>% mutate(sample = fct_relevel(sample, sample_name))

p1 <- ggplot(mpd2, aes(x=timepoint, y=values)) + geom_line(aes(group=sample), linetype = "dashed", colour = "grey", size = 0.5) + 
  geom_point(aes(color=sample), size=4) + theme_bw() + scale_color_manual(values =  sample_col) + theme(axis.title = element_blank())
p1
cowplot::ggsave2("./figures/singlecell_MPD_in_allsamples.pdf", p1, width = 3.75, height = 4)

#---statistic analysis---
mpd22 <- mpd2[order(mpd2$sample),]
pri_sample <- subset(mpd22, timepoint == "Primary", values, drop = TRUE)
rec_sample <- subset(mpd22, timepoint == "Recurrence", values, drop = TRUE)
wilcox.test(pri_sample, rec_sample, paired = TRUE)

#-----MPD downsample---(Figure S6b)---######
BiocParallel::register(BiocParallel::MulticoreParam(workers = 10, progressbar = T), default = T)
BiocParallel::bpparam()
set.seed(1)
# dir.create("./metrics/mpd_downsample")
mpd_downsam_res <- data.frame()
for(i in sample_name){
  print(paste0("Now is runnig: ", i))
  obj <- readRDS(paste0("./objects/", i, c("_final_filtered_overdisp_copykit.rds")))
  pri_cell <- obj@colData %>% as.data.frame() %>% dplyr::filter(timepoint == "primary") %>% pull(sample)
  rec_cell <- obj@colData %>% as.data.frame() %>% dplyr::filter(timepoint == "recurrence") %>% pull(sample)
  pri_df <- as_tibble(obj@assays@data$segment_ratios[, pri_cell])
  rec_df <- as_tibble(obj@assays@data$segment_ratios[, rec_cell])
  
  Primary <- BiocParallel::bplapply(seq(50,950,25), function(n){
    df_list <- sample_select(pri_df, n)
    if(!is.null(df_list)){
      mpds <- BiocParallel::bplapply(df_list, function(df_sub){
        return(mpd_scTree(df = df_sub,n_threads = 1))
      })
      return(c("primary",n,mean(unlist(mpds))))
    }
  })
  
  Recurrence <- BiocParallel::bplapply(seq(50,950,25), function(n){
    df_list <- sample_select(rec_df, n)
    if(!is.null(df_list)){
      mpds <- BiocParallel::bplapply(df_list, function(df_sub){
        return(mpd_scTree(df = df_sub,n_threads = 1))
      })
      return(c("recurrence",n,mean(unlist(mpds))))
    }
  })
  
  my_res <- do.call(rbind, c(Primary, Recurrence)) %>% as.data.frame() %>% cbind(i)
  mpd_downsam_res <- rbind(mpd_downsam_res, my_res)
  
  colnames(my_res) <- c("timepoint","cell_num","mpd","sample")
  write.csv(my_res, paste0("./metrics/mpd_downsample/mpd_downsample_", i, ".csv"))
}
colnames(mpd_downsam_res) <- c("timepoint","cell_num","mpd","sample")
# write.csv(mpd_downsam_res, "./metrics/mpd_downsample_all_samples.csv")
mpd_downsam_res <- read.csv("./metrics/mpd_downsample_all_samples.csv", row.names = 1)

colnames(res) <- c("Timepoint", "nCells", "singlecell_mpd","sample")
sample_col <- c("#565FAB","#CF3D32","#5DB1DD","#812268","#78C36A","#486A84","#BB6438","#F0E584","#D595A7","#749B59")
names(sample_col) <- sample_name

p1 <- ggplot(mpd_downsam_res, aes(x=cell_num, y=mpd)) +
  geom_line(aes(linetype=timepoint, colour=sample), size = 0.5) +
  geom_point(aes(shape=timepoint, colour=sample), size = 1.8) +
  scale_color_manual(values = sample_col) + scale_x_continuous(breaks = c(0, 200, 400, 600, 800)) +
  theme_bw() + ylab("MPD (calculated by single cell)")

cowplot::ggsave2("./figures/single_cell_mpd_downsample.pdf", p1, width = 5.5, height = 4)

#-----FACS ploidy (Figure 3e)-----#####
facs_ploid <- read.csv("./metrics/facs_ploidy.csv", row.names = 1)
facs_ploid2 <- facs_ploid %>% mutate(paired=rep(1:10, each=2))
facs_ploid2$sample <- factor(facs_ploid$sample, levels = rev(sample_name))

p1 <- ggplot(facs_ploid2, aes(x=facs_ploidies, y=sample)) + geom_line(aes(group=paired)) + 
  geom_point(aes(color=timepoint), size=4) + xlim(0,6) +
  theme_bw() + scale_color_manual(values =  tp_col)
p1
cowplot::ggsave2("./figures/facs_ploidy_of_allsamples.pdf", p1, width = 4, height = 4)

#-----# CNA, #subclone, MPD vs clinical features-(Figure S6c)----######
meta_mtx <- read.csv("./metrics/clinical_metadata_table.csv")
din_age_meta <- read.csv("./metrics/din_block_ages.csv", row.names = 1)
sclone <- read.csv("./metrics/subclone_number.csv", header = T)
cna_num <- read.csv("./metrics/number_cna_events.csv", row.names = 1)
mpd_res <- read.csv("./metrics/mpd_scTree.csv", row.names = 1)
facs_ploid <- read.csv("./metrics/facs_ploidy.csv", row.names = 1)
cell_num <- read.csv("./metrics/cell_number.csv", row.names = 1)

sclone2 <- sclone %>% column_to_rownames(var = "timepoints") %>% t() %>% as.data.frame() %>%
  dplyr::mutate(Primary = primary + shared, Recurrence= recurrence+ shared) %>%
  dplyr::select(Primary,Recurrence) %>% rownames_to_column() %>% gather(key = "timepoints", value = "subclone", -rowname) %>%
  dplyr::mutate(sample = paste0(rowname, ifelse(timepoints == "Primary", "_primary", "_recurrence")))
cna_num2 <- cna_num %>% mutate(sample2 = paste0(sample, ifelse(timepoint == "Primary", "_primary", "_recurrence")))
mpd2 <- gather(mpd_res, key = "timepoint", value = "values", 1:2)  %>% 
  dplyr::mutate(sample2 = paste0(sample, ifelse(timepoint == "Primary", "_primary", "_recurrence"))) %>% 
  dplyr::rename(mpd = values)
facs_ploid2 <- facs_ploid %>% dplyr::mutate(sample2 = paste0(sample, ifelse(timepoint == "Primary", "_primary", "_recurrence")))
cell_num2 <- cell_num %>% dplyr::rename(cell_num = Freq)

subclo_cna_meta <- inner_join(meta_mtx, sclone2, by = c("sample" = "sample")) %>%
  inner_join(din_age_meta, by = c("sample" = "patient")) %>% 
  inner_join(cna_num2, by = c("sample" = "sample2")) %>%
  inner_join(mpd2, by = c("sample" = "sample2")) %>%
  inner_join(facs_ploid2, by = c("sample" = "sample2")) %>%
  inner_join(cell_num2, by = c("sample" = "Var1")) %>%
  dplyr::select("sample","timepoint","ER","PR","HER2","histology","grade","subclone","DIN","block_age_years","num_cna",
                "mpd","facs_ploidies","cell_num")

write.csv(subclo_cna_meta, "./metrics/paired_clin_comp_meta.csv")
# subclo_cna_meta <- read.csv("./metrics/paired_clin_comp_meta.csv", row.names = 1)

subclo_cna_meta2 <- subclo_cna_meta %>% mutate(PR = na_if(PR, "Unk"), HER2 = na_if(HER2, "Unk"), 
                                               histology = ifelse(histology == "DCIS", "DCIS", "IBC"))
write.csv(subclo_cna_meta2, "./metrics/paired_clin_comp_meta_NA.csv")
# subclo_cna_meta2 <- read.csv("./metrics/paired_clin_comp_meta_NA.csv", row.names = 1)

my_clin <- c("ER","PR","HER2","histology", "grade")
plist <- list()
j <- 1
for (i in my_clin) {
  if(i %in% c("PR", "HER2")){
    p2 <- subclo_cna_meta2 %>% drop_na(all_of(i)) %>% ggboxplot(x = i, y = "mpd", color = "timepoint", palette = tp_col, 
                    add = "jitter", ylab = "mpd", xlab = i) + stat_compare_means()
  }else{
    p2 <- ggboxplot(subclo_cna_meta2, x = i, y = "mpd", color = "timepoint", palette = tp_col, 
                    add = "jitter", ylab = "mpd", xlab = i) + stat_compare_means()
  }
  plist[[j]] <- p2
  j <- j + 1
}

plist2 <- list()
j <- 1
for (i in my_clin) {
  if(i %in% c("PR", "HER2")){
    p4 <- subclo_cna_meta2 %>% drop_na(all_of(i)) %>% ggboxplot(x = i, y = "num_cna", color = "timepoint", palette = tp_col, 
                    add = "jitter", ylab = "num_cna", xlab = i) + stat_compare_means() 
  }else{
    p4 <- ggboxplot(subclo_cna_meta2, x = i, y = "num_cna", color = "timepoint", palette = tp_col, 
                    add = "jitter", ylab = "num_cna", xlab = i) + stat_compare_means() 
  }
  plist2[[j]] <- p4
  j <- j + 1
}

p_list<-plot_grid(plotlist = c(plist, plist2), ncol = 5)
cowplot::ggsave2("./figures/MPD_NumCNA_tp_er_pr_her2_his_grade2.pdf", p_list, width = 15, height = 6, limitsize = F)


subclo_cna_meta2 %>% filter(timepoint == "primary") %>% pull(subclone) %>% summary()
subclo_cna_meta2 %>% filter(timepoint == "recurrence") %>% pull(subclone) %>% summary()
subclo_cna_meta2 %>% filter(timepoint == "primary") %>% pull(cell_num) %>% summary()
subclo_cna_meta2 %>% filter(timepoint == "recurrence") %>% pull(cell_num) %>% summary()
subclo_cna_meta2 %>% filter(timepoint == "primary") %>% pull(cell_num) %>% sd()
subclo_cna_meta2 %>% filter(timepoint == "recurrence") %>% pull(cell_num) %>% sd()


#-----RCA-PCA--(Figure S10c)-#######
pro_name <- c("duke248_p3","duke254_p4","nki23_p5","nki26_p6","nki28_p7","nki19_p8","nki22_p9")
pca <- c("internal_2", "internal_7", "internal_4", "c11", "internal_1", "internal_18", "internal_2")
rca <- c("c4", "internal_4", "internal_8", "internal_2", "internal_7", "internal_14", "internal_6")
pca_rca <- cbind(pro_name, pca, rca) %>% data.frame()

write.csv(pca_rca, "./metrics/pca_rca_nodes.csv")
# pca_rca <- read.csv("./metrics/pca_rca_nodes.csv", row.names = 1)
#-----read data--
all_rca_pca <- data.frame()
for (i in 1:length(pro_name)) {
  # i <- 1
  tsv_out_path <- paste0("./metrics/medicc_files/", pro_name[i], "/")
  tree_profile <- read_tsv(paste0(tsv_out_path, pro_name[i],"_medicc2_input_final_cn_profiles.tsv"))
  
  mypca <- pca_rca %>% dplyr::filter(pro_name == all_of(pro_name[i])) %>% pull(pca)
  myrca <- pca_rca %>% dplyr::filter(pro_name == all_of(pro_name[i])) %>% pull(rca)
  
  mrca_mat <- tree_profile %>% dplyr::filter(sample_id %in% c(mypca,myrca)) %>%
    dplyr::select(-starts_with("is_")) %>% spread(sample_id, CN) %>% dplyr::select(all_of(mypca),all_of(myrca))
  
  colnames(mrca_mat) <- c("PCA", "RCA")
  myr_p <- mrca_mat %>% dplyr::mutate(rca_pca = RCA - PCA) %>% dplyr::pull(rca_pca)
  
  all_rca_pca <- rbind(all_rca_pca, myr_p)
}

rownames(all_rca_pca) <- pro_name
write.table(all_rca_pca, file = "./metrics/all_rca_pca_mtx.csv", sep = "\t", quote = F)
# all_rca_pca <- read.table(file = "./metrics/all_rca_pca_mtx.csv", sep = "\t", header = T)

#----plot consensus heatmap---
#-----header
ha_col=HeatmapAnnotation(foo=anno_text(chr_name, rot = 0, gp = gpar(fontsize =10)), df =chr_color, 
                         col = list(chr=c("1"="black", "2"="grey")), show_legend = F, annotation_name_side = "left")
#-----annotated genes---
keep_gene <- c("PIK3CA","CCNE2","MYC", "BRCA2", "RB1","ZNF217","AURKA")
gene_bin_sel <- gene_bin %>% dplyr::filter(gene %in% keep_gene)
ha_bottom_cs = columnAnnotation(clonal_state = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                         side = "bottom", labels_gp = gpar(fontsize = 14)))

sample_col <- c("#565FAB","#CF3D32","#5DB1DD","#812268","#78C36A","#486A84","#BB6438")
names(sample_col) <- pro_name
ha_row_cs=rowAnnotation(sample_name = pro_name, col = list(sample_name = sample_col), show_annotation_name = F)

pdiff <- max(all_rca_pca) - min(all_rca_pca)
col_vec_cs <- structure(pals::ocean.balance(length(0:pdiff)), names = min(all_rca_pca): max(all_rca_pca))

anno_mtx_pos <- apply(all_rca_pca, 2, function(x){sum(x>0)})
anno_mtx_neg <- apply(all_rca_pca, 2, function(x){sum(x<0)})
anno_mtx_zero <- apply(all_rca_pca, 2, function(x){sum(x==0)})

anno_mtx <- rbind(anno_mtx_pos, anno_mtx_zero, anno_mtx_neg)
write.table(anno_mtx, file = "./metrics/all_rca_pca_gain_and_loss_mtx.csv", sep = "\t", quote = F)
# anno_mtx <- read.table(file = "./metrics/all_rca_pca_gain_and_loss_mtx.csv", sep = "\t", header = T)

ha_bottom_cs = columnAnnotation(shared_event = anno_barplot(t(anno_mtx), 
                                                            gp = gpar(lwd = 0, fill = c("anno_mtx_pos" = "red", 
                                                                                        "anno_mtx_zero" = "grey",
                                                                                        "anno_mtx_neg" = "blue"), 
                                                                      col= c("anno_mtx_pos" = "red", 
                                                                             "anno_mtx_zero" = "grey", 
                                                                             "anno_mtx_neg" = "blue")), 
                                                            bar_width = 1), border = FALSE,
                                gene_mark = anno_mark(at=gene_bin_sel$pos, labels = gene_bin_sel$gene, 
                                                      side = "bottom", labels_gp = gpar(fontsize = 14)))

pdf(paste0("./figures/all_RCA_PCA_complexHeatmap_consensus_integer_with_barplot.pdf"), height = 4.5, width = 10)
Heatmap(all_rca_pca, cluster_columns = FALSE, border = TRUE, cluster_rows = T, show_row_dend = F,
        name = "csheatmap", show_row_names = T, row_names_side = "left",
        show_column_names = F, column_title = "RCA-PCA",
        use_raster = T, raster_quality = 5, 
        left_annotation = ha_row_cs, col = rev(col_vec_cs), 
        heatmap_legend_param = list(title = "copy number", title_gp = gpar(fontsize = 12, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12)),
        top_annotation = ha_col, bottom_annotation = ha_bottom_cs)
dev.off()



