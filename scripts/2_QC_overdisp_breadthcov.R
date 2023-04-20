#####-----------Process-------######
library(dplyr)
library(tidyr)
library(cowplot)
library(useful)
library(tibble)
# library(copykit)
library(tidyverse)
library(ggplot2)
require(RColorBrewer)
library(Homo.sapiens)
library(ggpubr)
options(max.print = 200)

setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload3")
source("./scripts/0_arcwell_functions.R")
load("./pre_load_data/pre_load_data.rda")
gene_bin <- read.table("./pre_load_data/hg19_gene_binpos_map.tsv", header = 1)
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", 
            "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42","#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE",
            "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0", "#E6E6FA", "#8B7E66")
new_pal2 <- new_pal
names(new_pal2) <- paste0("c", 1:length(new_pal2))
tp_col <- c("#EA3291", "#96C942")
colorpal <- c("#66C5CC","#8BE0A4","#F6CF71","#F89C74","#DCB0F2","#87C55F","#9EB9F3","#FE88B1","#C9DB74","#B497E7")

pro_name <- c("arc_well_mda231", "tenx", "dop_merge","dlp_merge","dlp_plus","act_mda231","fresh_315a","formalin_315a",
              "fresh_mda231","formalin_mda231","bcis28_frozen_qc","bcis28_ffpe_qc")
tech_name <- c("Arc-well-MDAMB231", "10X-CNA", "DOP-PCR", "DLP", "DLP+", "ACT-MDA231","Arc-well-fresh-315A", 
               "Arc-well-formalin-315A","Arc-well-fresh-MDAMB231", "Arc-well-formalin-MDAMB231","BCIS28-frozen","BCIS28-FFPE")
raw_path0 <- paste0("./map_seg_output/", pro_name)

####----calculate over dispersion for cell line and QC data------
#-----create copykit and meta columns of filtered cells---
all_meta <- data.frame()
for (i in 1:length(pro_name)) {
print(paste0("Now is runnig: ", pro_name[i]))

filter_cells = read.table(paste0("./metrics/", pro_name[i], "_filtered_bincounts_newnormal.txt"), header = T, check.names = F)
meta_file <- read.table(paste0("./metrics/", pro_name[i], "_metadata.metrics_newnormal.txt"), header = T, check.names = F) 

#-----select filtered cells
varbin_mtx <- readVarbinCNA(paste0("./map_seg_output/", pro_name[i]), remove_Y = TRUE, clean_names = F)
filt_cells_names <- colnames(filter_cells)[4:ncol(filter_cells)] 
varbin_mtx_filter <- varbin_mtx[,filt_cells_names]
merged_meta <- as.data.frame(varbin_mtx_filter@colData) 

#---add meta data
rownames(meta_file) <- meta_file$sample
meta_file2 <- meta_file %>% mutate(dups_percentage = DupsRemoved/TotalReads) %>% mutate(sample_name = pro_name[i]) %>% 
  rownames_to_column() %>% 
  dplyr::select(c("rowname", "filter_corr_value","ReadsKept","MedianBinCount","dups_percentage", "sample_name"))
merged_meta2 <- merged_meta %>% left_join(meta_file2, by = c("sample" = "rowname"))
varbin_mtx_filter@colData <- cbind(varbin_mtx_filter@colData, 
                                   merged_meta2[,c("filter_corr_value","ReadsKept","MedianBinCount","dups_percentage", "sample_name")])

#----calculate over dispersion---
bin_count <- varbin_mtx_filter@assays@data$bin_counts
bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% dplyr::rename(over_disp = V1) %>% dplyr::select("over_disp")
varbin_mtx_filter@colData <- cbind(varbin_mtx_filter@colData, bin_count_overdisp2) 

write.table(varbin_mtx_filter@colData , paste0("./metrics/", pro_name[i],"_basicQC_meta.txt"), sep = "\t", quote = F, row.names = F)
all_meta <- rbind(all_meta, as.data.frame(varbin_mtx_filter@colData))
}
write.table(all_meta, file = "./metrics/QCsamples_merged_basicQC_meta.txt", sep = "\t", quote = F, row.names = F)
# all_meta <- read.table(file = "./metrics/QCsamples_merged_basicQC_meta.txt", sep = "\t", header = T)

all_meta$tech <- plyr::mapvalues(all_meta$sample_name, from = pro_name, to = tech_name)
write.table(all_meta, file = "./metrics/QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", quote = F, row.names = F)
# all_meta <- read.table(file = "./metrics/QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", header = T)

#----sampling cells ----
# sampling cells to make group sizes more similar
set.seed(32)
all_meta_s80 <- all_meta %>% group_by(sample_name) %>% sample_n(80) %>% ungroup()

write.table(all_meta_s80, "./metrics/sampled_80cells_QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", quote = F, row.names = F)
# all_meta_s80 <- read.table(file = "./metrics/sampled_80cells_QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", header = T)

#----tech comparing QC and Figure 2a-----
all_cov_s80_comp <- all_meta_s80 %>% dplyr::filter(!(sample_name %in% c("formalin_315a","formalin_mda231","bcis28_frozen_qc","bcis28_ffpe_qc"))) %>% 
  mutate(tech = as.factor(tech)) %>% 
  mutate(tech = fct_relevel(tech, c("Arc-well-fresh-315A","Arc-well-fresh-MDAMB231","Arc-well-MDAMB231","ACT-MDA231","10X-CNA", "DLP", "DLP+", "DOP-PCR")))

all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-315A") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-315A") %>% pull(over_disp) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-MDAMB231") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-MDAMB231") %>% pull(over_disp) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-MDAMB231") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-MDAMB231") %>% pull(over_disp) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "ACT-MDA231") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "ACT-MDA231") %>% pull(over_disp) %>% mad()


ggboxplot(all_cov_s80_comp, x = "tech", y = "over_disp",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "Arc-well-fresh-315A")     

ggboxplot(all_cov_s80_comp, x = "tech", y = "over_disp",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "Arc-well-MDAMB231")


p1 <- ggplot(all_cov_s80_comp) + 
  ggbeeswarm::geom_quasirandom(aes(x = tech, y = over_disp, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + scale_fill_manual(values = colorpal) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),legend.position = "none", 
        strip.background = element_rect(fill = "white"))  +
  ylab("Overdispersion") + xlab("") 
p1

cowplot::ggsave2("./figures/tech_comp_over_dispersion.pdf", p1, width = 4.5, height = 4.8)

####----calculate breadth of coverage for cell line and QC data------
pro_name <- c("arc_well_mda231", "tenx", "dop_merge","dlp_merge","dlp_plus","act_mda231","fresh_315a","formalin_315a",
              "fresh_mda231","formalin_mda231","bcis28_frozen_qc","bcis28_ffpe_qc")
tech_name <- c("Arc-well-MDAMB231", "10X-CNA", "DOP-PCR", "DLP", "DLP+", "ACT-MDA231","Arc-well-fresh-315A", 
               "Arc-well-formalin-315A","Arc-well-fresh-MDAMB231", "Arc-well-formalin-MDAMB231","BCIS28-frozen","BCIS28-FFPE")

#-----find bam files with over 500k reads---####
bam_path <- c("/volumes/seq/projects/CNA_projects/DT_CNA/cell_line/231/MDAMB231/MDAMB231/MDAMB231_P1_P2_P3/output/sort/",
              "/volumes/seq/external_data/laks_2019/dlp_plus/varbin_200kb/output/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/dlp_reprocess/dlp_merged_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/dop_pcr_reprocess/dop_pcr_merged_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/20210206_MDA231/data/MDA231-9x_200/res_200_k/sort/",
              "/volumes/seq/projects/CNA_projects/10X_CNA/Breast/TN7_TN17/TN17_CBS_output/CBS_output/output/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/20220515_315A_nofix/data/Arc_well_315A_nofix_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/20220517_315A_formalin/data/Arc_well_315A_formalin_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/20230327_mda231_nonfix/data/Arc_well_MDA231_nonfix_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/20230328_mda231_fixed/data/Arc_well_MDA231_fixed_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/20220516_BCIS28_frozen/data/BCIS28_frozen_200/res_200_k/sort/",
              "/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/20220518_BCIS28_FFPE_reseq/data/BCIS28_FFPE_200/res_200_k/sort/")

for (i in 1:12) {
  my_pro_name <- pro_name[i]
  my_bam_path <- bam_path[i]
  bin_f <- read.table(paste0("./metrics/", my_pro_name, "_filtered_bincounts_newnormal.txt"), header = T)
  output_path <- paste0("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/QC_downsample/", my_pro_name, "/")
  #----Find bam with over 500K reads---
  system(paste0("mkdir -p /volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/QC_downsample/", my_pro_name, "/data"))
  link_bam_files(bin_f, my_bam_path, output_path, my_pro_name, target_reads = 500000)
  #---link bam with over 500k reads---
  system(paste0("cd /volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/QC_downsample/", my_pro_name, 
                "/data; cat ../", my_pro_name, " | xargs -I % ln -s % ."))
}

#------downsample bam files to 500k reads and calculate the breadth of coverage---
#----run the snakemake pipeline in terminal--
# source ~/.bashrc
# source ~/.bash_profile
# conda activate snakemake
# snakemake --snakefile ./scripts/snakemake_files/ds_breadth.smk --cores 20

#----merge breadth of coverage data of all QC samples--
dir.create("./metrics/breadth_cov")
all_cov_breadth <- data.frame()
for (i in 1:length(pro_name)) {
  print(paste0("Now is runnig: ", pro_name[i]))
  my_cov_path <- paste0("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/data/QC_downsample/", pro_name[i], "/covfile/")
  my_cov <- calc_coverage(path = my_cov_path) %>% mutate(sample = pro_name[i])
  write.table(my_cov, paste0("./metrics/breadth_cov/",pro_name[i], "_coverage_breadth.txt"), sep = "\t", quote = F, row.names = F)
  all_cov_breadth <- rbind(all_cov_breadth, my_cov)
}

write.table(all_cov_breadth, "./metrics/QCsamples_merged_coverage_breadth.txt", sep = "\t", quote = F, row.names = F)
# all_cov_breadth <- read.table(file = "./metrics/QCsamples_merged_coverage_breadth.txt", sep = "\t", header = T)

all_cov_breadth$tech <- plyr::mapvalues(all_cov_breadth$sample, from = pro_name, to = tech_name)
write.table(all_cov_breadth, file = "./metrics/QCsamples_merged_coverage_breadth_with_tech.txt", sep = "\t", quote = F, row.names = F)
# all_cov_breadth <- read.table(file = "./metrics/QCsamples_merged_coverage_breadth_with_tech.txt", sep = "\t", header = T)

# sampling cells to make group sizes more similar
set.seed(14)
all_cov_s80 <- all_cov_breadth %>% group_by(sample) %>% sample_n(80) %>% ungroup()
write.table(all_cov_s80, "./metrics/sampled_80cells_QCsamples_merged_coverage_breadth_with_tech.txt", 
            sep = "\t", quote = F, row.names = F)
# all_cov_s80 <- read.table("./metrics/sampled_80cells_QCsamples_merged_coverage_breadth_with_tech.txt", sep = "\t", header = T)
#---plotting figure 2b----#####
all_cov_s80_comp <- all_cov_s80 %>% dplyr::filter(!(sample %in% c("formalin_315a","formalin_mda231","bcis28_frozen_qc","bcis28_ffpe_qc"))) %>% 
  mutate(tech = as.factor(tech)) %>% 
  mutate(tech = fct_relevel(tech, c("Arc-well-fresh-315A","Arc-well-fresh-MDAMB231","Arc-well-MDAMB231","ACT-MDA231","10X-CNA", "DLP", "DLP+", "DOP-PCR")))

p1 <- ggplot(all_cov_s80_comp) + ggbeeswarm::geom_quasirandom(aes(x = tech, y = breadth, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + scale_fill_manual(values = colorpal) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), 
        legend.position = "none", strip.background = element_rect(fill = "white"))  +
  # ggtitle("Index of Dispersion") +
  ylab("Breadth of coverage") +
  xlab("") 
p1
cowplot::ggsave2("./figures/tech_comp_breadth_cov.pdf", p1, width = 4.5, height = 4.8)

all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-315A") %>% pull(breadth) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-315A") %>% pull(breadth) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-MDAMB231") %>% pull(breadth) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-MDAMB231") %>% pull(breadth) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "ACT-MDA231") %>% pull(breadth) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "ACT-MDA231") %>% pull(breadth) %>% mad()


ggboxplot(all_cov_s80_comp, x = "tech", y = "breadth",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "Arc-well-fresh-315A")     

ggboxplot(all_cov_s80_comp, x = "tech", y = "breadth",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "Arc-well-MDAMB231")

#---plotting figure 2c--315A no-fix (fresh) vs formalin fix-------
arc315_meta <- all_meta %>% dplyr::filter(sample_name %in% c("formalin_315a","fresh_315a")) %>% 
  mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, c("Arc-well-fresh-315A","Arc-well-formalin-315A")))
arc315_cov <- all_cov_breadth %>% dplyr::filter(sample %in% c("formalin_315a","fresh_315a")) %>% 
  mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, c("Arc-well-fresh-315A","Arc-well-formalin-315A")))

sample_col <- c("#1D4E89","#fb5607")
p1 <- ggviolin(arc315_meta, x = "tech", y = "over_disp", fill = "tech", palette = sample_col, size = 0.1) + 
  geom_boxplot(outlier.shape = NA,width = 0.1) + stat_compare_means() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  ggtitle("315A non-fixed vs formalin fix")
p2 <- ggviolin(arc315_cov, x = "tech", y = "breadth", fill = "tech", palette = sample_col, size = 0.1) + 
  stat_compare_means()+ theme(legend.position = "none", axis.title.x = element_blank()) 
p3 <- plot_grid(plotlist=list(p1,p2), ncol=1, align='v')
p3
cowplot::ggsave2("./figures/315A_frozen_formalin_violin_overdisp_breadth.pdf", p3, width = 2.6, height = 4)

#------315A QC metrics---
median((arc315_meta %>% dplyr::filter(tech == "Arc-well-fresh-315A"))$over_disp)
mad((arc315_meta %>% dplyr::filter(tech == "Arc-well-fresh-315A"))$over_disp)
median((arc315_meta %>% dplyr::filter(tech == "Arc-well-formalin-315A"))$over_disp)
mad((arc315_meta %>% dplyr::filter(tech == "Arc-well-formalin-315A"))$over_disp)

median((arc315_cov %>% dplyr::filter(tech == "Arc-well-fresh-315A"))$breadth)
mad((arc315_cov %>% dplyr::filter(tech == "Arc-well-fresh-315A"))$breadth)
median((arc315_cov %>% dplyr::filter(tech == "Arc-well-formalin-315A"))$breadth)
mad((arc315_cov %>% dplyr::filter(tech == "Arc-well-formalin-315A"))$breadth)
ggboxplot(all_meta, x = "tech", y = "over_disp",color = "tech") + stat_compare_means(ref.group = "Arc-well-formalin-315A")     
ggboxplot(all_cov_breadth, x = "tech", y = "breadth",color = "tech") + stat_compare_means(ref.group = "Arc-well-formalin-315A")

#---plotting figure 2c--MDAMB231 no-fix (fresh) vs formalin fix-------
arcmda231_meta <- all_meta %>% dplyr::filter(sample_name %in% c("formalin_mda231","fresh_mda231")) %>% 
  mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, c("Arc-well-fresh-MDAMB231","Arc-well-formalin-MDAMB231")))
arcmda231_cov <- all_cov_breadth %>% dplyr::filter(sample %in% c("formalin_mda231","fresh_mda231")) %>% 
  mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, c("Arc-well-fresh-MDAMB231","Arc-well-formalin-MDAMB231")))

sample_col <- c("#1D4E89","#fb5607")
p11 <- ggviolin(arcmda231_meta, x = "tech", y = "over_disp", fill = "tech", palette = sample_col, size = 0.1) + 
  geom_boxplot(outlier.shape = NA,width = 0.1) + stat_compare_means() + 
  theme(legend.position = "none", axis.title.x = element_blank()) 

p12 <- ggviolin(arcmda231_cov, x = "tech", y = "breadth", fill = "tech", palette = sample_col, size = 0.1) + 
  stat_compare_means()+ theme(legend.position = "none", axis.title.x = element_blank()) 

p13 <- plot_grid(plotlist=list(p11,p12), ncol=2, align='v')
p13
cowplot::ggsave2("./figures/mda231_frozen_formalin_violin_overdisp_breadth.pdf", p13, width = 5, height = 2.5)

#------MDAMB231 QC metrics---
a <- median((arcmda231_meta %>% dplyr::filter(tech == "Arc-well-fresh-MDAMB231"))$over_disp)
mad((arcmda231_meta %>% dplyr::filter(tech == "Arc-well-fresh-MDAMB231"))$over_disp)
b <- median((arcmda231_meta %>% dplyr::filter(tech == "Arc-well-formalin-MDAMB231"))$over_disp)
mad((arcmda231_meta %>% dplyr::filter(tech == "Arc-well-formalin-MDAMB231"))$over_disp)

c <- median((arcmda231_cov %>% dplyr::filter(tech == "Arc-well-fresh-MDAMB231"))$breadth)
mad((arcmda231_cov %>% dplyr::filter(tech == "Arc-well-fresh-MDAMB231"))$breadth)
d <- median((arcmda231_cov %>% dplyr::filter(tech == "Arc-well-formalin-MDAMB231"))$breadth)
mad((arcmda231_cov %>% dplyr::filter(tech == "Arc-well-formalin-MDAMB231"))$breadth)

compar_group1 <- c("Arc-well-fresh-MDAMB231", "ACT-MDA231","10X-CNA", "DLP","DLP+","DOP-PCR")
compar_group2 <- c("Arc-well-formalin-MDAMB231","ACT-MDA231","10X-CNA", "DLP","DLP+","DOP-PCR")
all_meta_tmp <-  all_meta %>% dplyr::filter(tech %in% compar_group1) %>% mutate(tech = factor(tech, level = compar_group1))
all_meta_tmp2 <-  all_meta %>% dplyr::filter(tech %in% compar_group2) %>% mutate(tech = factor(tech, level = compar_group2))

all_cov_breadth_tmp <-  all_cov_breadth %>% dplyr::filter(tech %in% compar_group1) %>% mutate(tech = factor(tech, level = compar_group1))
all_cov_breadth_tmp2 <-  all_cov_breadth %>% dplyr::filter(tech %in% compar_group2) %>% mutate(tech = factor(tech, level = compar_group2))

colorpal_fresh <- c("#1D4E89","#F6CF71","#F89C74","#DCB0F2","#87C55F","#9EB9F3")
colorpal_formalin <- c("#fb5607","#F6CF71","#F89C74","#DCB0F2","#87C55F","#9EB9F3")

p1 <- ggplot(all_meta_tmp,  aes(x = tech, y = over_disp, color=tech)) + geom_boxplot() + geom_hline(yintercept = a, linetype="dashed", color="red") +
  theme(legend.position = "none") + scale_color_manual(values = colorpal_fresh)
p2 <- ggplot(all_meta_tmp2,  aes(x = tech, y = over_disp, color=tech)) + geom_boxplot() + geom_hline(yintercept = b, linetype="dashed", color="red") +
  theme(legend.position = "none") + scale_color_manual(values = colorpal_formalin)

p3 <- ggplot(all_cov_breadth_tmp,  aes(x = tech, y = breadth, color=tech)) + geom_boxplot() + geom_hline(yintercept = c, linetype="dashed", color="red") +
  theme(legend.position = "none") + scale_color_manual(values = colorpal_fresh)
p4 <- ggplot(all_cov_breadth_tmp2,  aes(x = tech, y = breadth, color=tech)) + geom_boxplot() + geom_hline(yintercept = d, linetype="dashed", color="red") +
  theme(legend.position = "none") + scale_color_manual(values = colorpal_formalin)

p5 <- plot_grid(plotlist=list(p1,p2,p3,p4), ncol=2, align='v')
p5
cowplot::ggsave2("./figures/mda231_frozen_formalin_vs_other_scDNA_overdisp_breadth.pdf", p5, width = 8, height = 5)

#---plotting figure 2c--bcis28 frozen vs ffpe-------
all_meta <- read.table(file = "./metrics/QCsamples_merged_basicQC_meta_with_tech.txt", sep = "\t", header = T)
all_cov_breadth <- read.table(file = "./metrics/QCsamples_merged_coverage_breadth_with_tech.txt", sep = "\t", header = T)

arcbcis28_meta <- all_meta %>% dplyr::filter(sample_name %in% c("bcis28_frozen_qc","bcis28_ffpe_qc")) %>% 
  mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, c("BCIS28-frozen","BCIS28-FFPE")))
arcbcis28_cov <- all_cov_breadth %>% dplyr::filter(sample %in% c("bcis28_frozen_qc","bcis28_ffpe_qc")) %>% 
  mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, c("BCIS28-frozen","BCIS28-FFPE")))

sample_col <- c("#1D4E89","#fb5607")
p4 <- ggviolin(arcbcis28_meta, x = "tech", y = "over_disp", fill = "tech", palette = sample_col, size = 0.1) + 
  geom_boxplot(outlier.shape = NA,width = 0.1) +
  stat_compare_means()+ theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ggtitle("BCIS28 Frozen vs FFPE")
p5 <- ggviolin(arcbcis28_cov, x = "tech", y = "breadth", fill = "tech", palette = sample_col, size = 0.1) + 
  stat_compare_means()+ theme(legend.position = "none", axis.title.x = element_blank()) 
p6 <- plot_grid(plotlist=list(p4,p5), ncol=1, align='v')
cowplot::ggsave2("./figures/BCIS28_frozen_vs_ffpe_violin_overdisp_breadth.pdf", p6, width = 2.6, height = 4)

p7 <- plot_grid(plotlist=list(p1,p11,p4,p2,p12,p5), ncol=3,  align='v')
cowplot::ggsave2("./figures/315A_BCIS28_frozen_vs_ffpe_violin_overdisp_breadth.pdf", p7, width = 7.5, height = 4)

#-----BCIS 28- QC- metrics
median((all_meta %>% dplyr::filter(tech == "BCIS28-frozen"))$over_disp)
mad((all_meta %>% dplyr::filter(tech == "BCIS28-frozen"))$over_disp)
median((all_meta %>% dplyr::filter(tech == "BCIS28-FFPE"))$over_disp)
mad((all_meta %>% dplyr::filter(tech == "BCIS28-FFPE"))$over_disp)

median((all_cov_breadth %>% dplyr::filter(tech == "BCIS28-frozen"))$breadth)
mad((all_cov_breadth %>% dplyr::filter(tech == "BCIS28-frozen"))$breadth)
median((all_cov_breadth %>% dplyr::filter(tech == "BCIS28-FFPE"))$breadth)
mad((all_cov_breadth %>% dplyr::filter(tech == "BCIS28-FFPE"))$breadth)
ggboxplot(all_meta, x = "tech", y = "over_disp",color = "tech") + stat_compare_means(ref.group = "BCIS28-FFPE")     
ggboxplot(all_cov_breadth, x = "tech", y = "breadth",color = "tech") + stat_compare_means(ref.group = "BCIS28-FFPE")

#------Overdispersion for all DCIS FFPE samples----#####
pro_name <- c("duke249_p1", "nki17_p2","duke248_p3","duke254_p4","nki23_p5","nki26_p6","nki28_p7",
              "nki19_p8","nki22_p9","nki15_p10","nki31_p11","nki12_p12")

#-----create copykit and meta columns of filtered cells-------
all_meta <- data.frame()
for (i in 1:length(pro_name)) {
  print(paste0("Now is runnig: ", pro_name[i]))
  varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name[i], c("_final_filtered_overdisp_copykit.rds")))
  
  my_meta <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% rownames_to_column() %>% 
    mutate(patient = paste0(pro_name[i],"_",timepoint)) %>% 
    dplyr::select(c("sample","reads_total","reads_assigned_bins","percentage_duplicates",
                    "dispense","peak","timepoint","subclones","tp_clst","over_disp","patient"))
  
  all_meta <- rbind(all_meta, my_meta)
}

write.table(all_meta, file = "./metrics/merged_QC_meta_all_DCIS_FFPE.txt", sep = "\t", quote = F, row.names = F)
all_meta <- read.table(file = "./metrics/merged_QC_meta_all_DCIS_FFPE.txt", sep = "\t", header = T)

#-----Overdispersion for all lung and prostate samples-----
pro_name2 <- c("lung30_p1","lung31_p2","prostate_p1_ax4bl","prostate_p2_pcf169")
#-----create copykit and meta columns of filtered cells---
all_meta_lp <- data.frame()
for (i in 1:length(pro_name2)) {
  print(paste0("Now is runnig: ", pro_name2[i]))
  varbin_mtx_tumor_log2 <- readRDS(paste0("./objects/", pro_name2[i], c("_final_filtered_overdisp_copykit.rds")))
  
  my_meta <- varbin_mtx_tumor_log2@colData %>% as.data.frame() %>% rownames_to_column() %>% 
    mutate(patient = paste0(pro_name2[i],"_",timepoint)) %>% 
    dplyr::select(c("sample","reads_total","reads_assigned_bins","percentage_duplicates",
                    "peak","timepoint","subclones","over_disp","patient"))
  all_meta_lp <- rbind(all_meta_lp, my_meta)
}

all_meta_lp$tp_clst <- paste(all_meta_lp$timepoint, all_meta_lp$subclones, sep="_")
write.table(all_meta_lp, file = "./metrics/merged_QC_meta_all_lung_prostate_FFPE.txt", sep = "\t", quote = F, row.names = F)
# all_meta_lp <- read.table(file = "./metrics/merged_QC_meta_all_lung_prostate_FFPE.txt", sep = "\t", header = T)
all_meta <- all_meta %>% dplyr::select(-("dispense")) %>% rbind(all_meta_lp)

write.table(all_meta, file = "./metrics/merged_QC_meta_all_DCIS_lung_prostate_FFPE.txt", sep = "\t", quote = F, row.names = F)
# all_meta <- read.table(file = "./metrics/merged_QC_meta_all_DCIS_lung_prostate_FFPE.txt", sep = "\t", header = T)

timepoint_col <- c("primary" = "#EA3291", "recurrence" = "#96C942", "single" = 'royalblue4')

pro_name2 <- c("duke248_p3","duke254_p4","nki23_p5","nki26_p6","nki28_p7","nki19_p8","nki22_p9","nki15_p10","nki31_p11","nki12_p12")
time_point <- c("primary","recurrence")
sample_levels <- c("duke249_p1_recurrence","nki17_p2_primary","lung30_p1_primary","lung31_p2_primary", "prostate_p1_ax4bl_primary",
                   "prostate_p2_pcf169_primary",paste(rep(pro_name2, each = length(time_point)), time_point, sep = "_"))

all_meta2 <- all_meta %>% dplyr::mutate(patient = fct_relevel(patient, sample_levels)) %>% 
  dplyr::mutate(class = ifelse(stringr::str_detect(patient, "duke249|nki17|lung|prostate"), "single", "paired")) %>% 
  dplyr::mutate(class = fct_relevel(class, c("single","paired"))) %>% 
  dplyr::mutate(timepoint2 = ifelse(stringr::str_detect(patient, "duke249|nki17|lung|prostate"), "single", timepoint))

p1 <- ggplot(all_meta2) + ggbeeswarm::geom_quasirandom(aes(x = patient, y = over_disp, fill = timepoint2), shape = 21, dodge.width = .8) +
  facet_grid(cols = vars(class), scales = "free_x", space = 'free_x') + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_fill_manual(values = timepoint_col) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),legend.position = "none", 
        strip.background = element_rect(fill = "white"))  +
  ylab("overdispersion") + xlab("") 

p1
# cowplot::ggsave2("./figures/QC_all_DCIS_FFPE_all_overdis.pdf", p1, width = 8, height = 4)
cowplot::ggsave2("./figures/QC_all_DCIS_lung_prostate_FFPE_all_overdis.pdf", p1, width = 8, height = 4.5)


#---Species mixture experiment A20 VS GM12878----#########
a20_gm <- read.table("./metrics/arc_well_a20_gm12878_hg19mm10_readcounts.txt")
a20_gm

library(tidyverse)
library(glue)
rc_cut <- 100000

count_tbl <- read_table("./metrics/arc_well_a20_gm12878_hg19mm10_readcounts.txt", col_names = F)
colnames(count_tbl) <- c("cellname","hg19", "mm10")

count_tbl %<>% mutate(p_hg19=hg19/(hg19+mm10),
                      p_mm10=mm10/(hg19+mm10),
                      tt=hg19+mm10) %>%
  mutate(class=if_else(p_mm10>0.95, "Mouse", if_else(p_hg19>0.95, "Human", "Collision"))) %>%
  filter(tt > rc_cut)

write.table(count_tbl,"./metrics/arc_well_a20_gm12878_hg19mm10_readcounts.txt")
# count_tbl <- read.table("./metrics/arc_well_a20_gm12878_hg19mm10_readcounts.txt")

cr <- round(sum(count_tbl$class=="Collision")/nrow(count_tbl)*100,2)

p1 <- ggplot(count_tbl, aes(x=hg19, y=mm10, color=class)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = c("gray","red","blue")) + 
  xlab("# of reads mapped to hg19") +
  ylab("# of reads mapped to mm10") +
  annotate("text", x=1e6, y=1e6, label= glue("{cr}% collision rate"), size = 5) +
  theme_bw()

cowplot::ggsave2("./figures/collison_plot_100K.pdf", p1, width = 6, height = 5)

#-----Poisson Distribution of cell dispensing---#########
pois_df <- data.frame()
for (i in 0:6) {
  tmp <- dpois(i, lambda=1)
  pois_df <- rbind(pois_df, c(i, tmp, 5184*tmp))
}
colnames(pois_df) <- c("cell_number_per_well","pos","expected_number_in_a_5184_wells_chip")

my_pale <- c("grey50", "#66bd63", "#d6604d", "#d6604d","#d6604d","#d6604d","#d6604d")
p1 <- ggplot(pois_df %>% mutate(cell_number_per_well = as.character(cell_number_per_well)), 
             aes(x=cell_number_per_well, y=expected_number_in_a_5184_wells_chip, fill=cell_number_per_well)) + 
  scale_fill_manual(values =  my_pale) +
  geom_bar(stat = "identity") + theme_classic() 
p1

cowplot::ggsave2("./figures/poisson_distribution_of_cell_dispensing.pdf", p1, width = 6, height = 4)

#------Doublets in A20-GM12878 species mixture experiment-----########
disp_file <- read.table("./metrics/waferDT_a20_gm12878_mix_fixed_5184.txt.txt", sep = "\t", header = T)
disp_file$State %>% table()

#-----Overdispersion comparison of unfiltered SCI-L3-WGS VS unfiltered Arc-well---########
pro_name <- c("fresh_315a","formalin_315a", "fresh_mda231", "formalin_mda231", "sci_l3_yi141")
tech_name <- c("Arc-well-fresh-315A", "Arc-well-formalin-315A", "Arc-well-fresh-MDAMB231", "Arc-well-formalin-MDAMB231",  "sci_L3_WGS")
raw_path0 <- paste0("./map_seg_output/", pro_name)

all_meta <- data.frame()
for (i in 1:length(pro_name)) {
  # i <- 5
  print(paste0("Now is runnig: ", pro_name[i]))
  meta_file <- read.table(paste0("./metrics/", pro_name[i], "_metadata.metrics_newnormal.txt"), header = T, check.names = F) 
  
  #-----select filtered cells
  varbin_mtx <- readVarbinCNA(paste0("./map_seg_output/", pro_name[i]), remove_Y = TRUE, clean_names = F)
  merged_meta <- as.data.frame(varbin_mtx@colData) 
  
  #---add meta data
  rownames(meta_file) <- meta_file$sample
  meta_file2 <- meta_file %>% mutate(dups_percentage = DupsRemoved/TotalReads) %>% mutate(sample_name = pro_name[i]) %>% 
    rownames_to_column() %>% 
    dplyr::select(c("rowname", "filter_corr_value","ReadsKept","MedianBinCount","dups_percentage", "sample_name"))
  merged_meta2 <- merged_meta %>% left_join(meta_file2, by = c("sample" = "rowname"))
  varbin_mtx@colData <- cbind(varbin_mtx@colData, 
                                     merged_meta2[,c("filter_corr_value","ReadsKept","MedianBinCount","dups_percentage", "sample_name")])
  
  #----calculate over dispersion---
  bin_count <- varbin_mtx@assays@data$bin_counts
  bin_count_overdisp <- map_dfr(bin_count, overdispersion) %>% t() %>% as.data.frame()
  bin_count_overdisp2 <- bin_count_overdisp %>% rownames_to_column() %>% dplyr::rename(over_disp = V1) %>% dplyr::select("over_disp")
  varbin_mtx@colData <- cbind(varbin_mtx@colData, bin_count_overdisp2) 
  
  write.table(varbin_mtx@colData , paste0("./metrics/", pro_name[i],"_unfiltered_cells_basicQC_meta.txt"), sep = "\t", quote = F, row.names = F)
  all_meta <- rbind(all_meta, as.data.frame(varbin_mtx@colData))
}
write.table(all_meta, file = "./metrics/arc_well_315a_mda231_sci_l3_unfiltered_cells_meta.txt", sep = "\t", quote = F, row.names = F)
# all_meta <- read.table(file = "./metrics/arc_well_315a_mda231_sci_l3_unfiltered_cells_meta.txt", sep = "\t", header = T)

all_meta$tech <- plyr::mapvalues(all_meta$sample_name, from = pro_name, to = tech_name)
write.table(all_meta, file = "./metrics/arc_well_315a_mda231_sci_l3_unfiltered_cells_meta_with_tech.txt", sep = "\t", quote = F, row.names = F)
# all_meta <- read.table(file = "./metrics/arc_well_315a_mda231_sci_l3_unfiltered_cells_meta_with_tech.txt", sep = "\t", header = T)

#----sampling cells ----
pro_name <- c("fresh_315a","formalin_315a", "fresh_mda231", "formalin_mda231", "sci_l3_yi141")
tech_name <- c("Arc-well-fresh-315A", "Arc-well-formalin-315A", "Arc-well-fresh-MDAMB231", "Arc-well-formalin-MDAMB231",  "sci_L3_WGS")
raw_path0 <- paste0("./map_seg_output/", pro_name)

# sampling cells to make group sizes more similar
set.seed(31)
all_meta_s80 <- all_meta %>% group_by(sample_name) %>% sample_n(80) %>% ungroup()

write.table(all_meta_s80, "./metrics/sampled_80cells_arc_well_315a_mda231_sci_l3_unfiltered_cells_meta_with_tech.txt", sep = "\t", quote = F, row.names = F)
# all_meta_s80 <- read.table(file = "./metrics/sampled_80cells_arc_well_315a_mda231_sci_l3_unfiltered_cells_meta_with_tech.txt", sep = "\t", header = T)

#----tech comparing QC with SCI_L3_WGS-----
all_cov_s80_comp <- all_meta_s80 %>%  mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, tech_name))

all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-315A") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-315A") %>% pull(over_disp) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-formalin-315A") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-formalin-315A") %>% pull(over_disp) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-MDAMB231") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-fresh-MDAMB231") %>% pull(over_disp) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-formalin-MDAMB231") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "Arc-well-formalin-MDAMB231") %>% pull(over_disp) %>% mad()
all_cov_s80_comp %>% dplyr::filter(tech == "sci_L3_WGS") %>% pull(over_disp) %>% median()
all_cov_s80_comp %>% dplyr::filter(tech == "sci_L3_WGS") %>% pull(over_disp) %>% mad()


ggboxplot(all_cov_s80_comp, x = "tech", y = "over_disp",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "sci_L3_WGS")     

ggboxplot(all_cov_s80_comp, x = "tech", y = "over_disp",color = "tech", palette = "colorpal")+
  stat_compare_means(ref.group = "sci_L3_WGS")


p1 <- ggplot(all_cov_s80_comp) + 
  ggbeeswarm::geom_quasirandom(aes(x = tech, y = over_disp, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + scale_fill_manual(values = colorpal) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),legend.position = "none", 
        strip.background = element_rect(fill = "white"))  +
  ylab("Overdispersion") + xlab("") 
p1

cowplot::ggsave2("./figures/sampled_80cells_arc_well_315a_mda231_sci_l3_unfiltered_cells_comp_over_dispersion.pdf", p1, width = 4, height = 4)


all_cov_comp <- all_meta %>%  mutate(tech = as.factor(tech)) %>% mutate(tech = fct_relevel(tech, tech_name))

p1 <- ggplot(all_cov_comp) + 
  ggbeeswarm::geom_quasirandom(aes(x = tech, y = over_disp, fill = tech), shape = 21, dodge.width = .8) +
  theme_cowplot() + scale_fill_manual(values = colorpal) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = 0.5),legend.position = "none", 
        strip.background = element_rect(fill = "white"))  +
  ylab("Overdispersion") + xlab("") 
p1


