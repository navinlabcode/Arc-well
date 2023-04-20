source("./scripts/0_arcwell_functions.R")
setwd("/volumes/USR2/wangkl/wafergen/DNA/ffpe_dcis/github_upload")
darlan_wd <- c("/volumes/lab/users/dminussi/projects/wafer_ffpe/")
n_reads_levels <- c('1M', '750k', '500k', '250k', '125k', '75k', '50k')

# Arc-well MDA-MB-231--(PCR duplicate rates: 16.43%)----
#-----get bam files that have over 1M reads---
# link_bam_files(bincounts = wafer_mda231p_bincounts,
#                path_bam_files = "/bam/file/path/",
#                output_path = "/downsample_corr/wafer231parental/",
#                output_name = "wafer_mda231p_dscorr.bamlist",
#                target_reads = 1000000)

#------downsample bam files ---
#----run the snakemake pipeline in terminal--
# source ~/.bashrc
# source ~/.bash_profile
# conda activate snakemake
# snakemake --snakefile ./scripts/snakemake_files/downsample_corr.smk --cores 20

#------run CNA_pipeline to get segmentation files on terminal--
#./run_CNA.sh --bam arc_well_231_bam_input_1M --sample arc_well_231_1M --output arc_well_231_1M --cpu 10
#./run_CNA.sh --bam arc_well_231_bam_input_750k --sample arc_well_231_750k --output arc_well_231_750k --cpu 10
#./run_CNA.sh --bam arc_well_231_bam_input_500k --sample arc_well_231_500k --output arc_well_231_500k --cpu 10
#./run_CNA.sh --bam arc_well_231_bam_input_250k --sample arc_well_231_250k --output arc_well_231_250k --cpu 10
#./run_CNA.sh --bam arc_well_231_bam_input_125k --sample arc_well_231_125k --output arc_well_231_125k --cpu 10
#./run_CNA.sh --bam arc_well_231_bam_input_75k --sample arc_well_231_75k --output arc_well_231_75k --filter_ReadCount=10000 --cpu 10
#./run_CNA.sh --bam arc_well_231_bam_input_50k --filter_ReadCount=10000 --sample arc_well_231_50k --output arc_well_231_50k --cpu 10

#----Plot figures---
# reading files segs
wafer_mda231p_original_seg <- read.table("/CNA_pipeline/outpath/final_result/uber.sample.seg.txt", header = T)
wafer_mda231p_1M_seg <- read.table("/CNA_pipeline/outpath/wafer231p_1M/final_result/uber.wafer231p_1M.seg.txt", header = T)
names(wafer_mda231p_1M_seg) <- str_remove(names(wafer_mda231p_1M_seg), '.sort')
wafer_mda231p_750k_seg <- read.table("/CNA_pipeline/outpath/wafer231p_750k/final_result/uber.wafer231p_750k.seg.txt", header = T)
names(wafer_mda231p_750k_seg) <- str_remove(names(wafer_mda231p_750k_seg), '.sort')
wafer_mda231p_500k_seg <- read.table("/CNA_pipeline/outpath/wafer231p_500k/final_result/uber.wafer231p_500k.seg.txt", header = T)
names(wafer_mda231p_500k_seg) <- str_remove(names(wafer_mda231p_500k_seg), '.sort')
wafer_mda231p_250k_seg <- read.table("/CNA_pipeline/outpath/wafer231p_250k/final_result/uber.wafer231p_250k.seg.txt", header = T)
names(wafer_mda231p_250k_seg) <- str_remove(names(wafer_mda231p_250k_seg), '.sort')
wafer_mda231p_125k_seg <- read.table("/CNA_pipeline/outpath/wafer231p_125k/final_result/uber.wafer231p_125k.seg.txt", header = T)
names(wafer_mda231p_125k_seg) <- str_remove(names(wafer_mda231p_125k_seg), '.sort')
wafer_mda231p_75k_seg <- read.table("/CNA_pipeline/outpath/wafer231p_75k/final_result/uber.wafer231p_75k.seg.txt", header = T)
names(wafer_mda231p_75k_seg) <- str_remove(names(wafer_mda231p_75k_seg), '.sort')
wafer_mda231p_50k_seg <- read.table("/CNA_pipeline/outpath/wafer231p_50k/final_result/uber.wafer231p_50k.seg.txt", header = T)
names(wafer_mda231p_50k_seg) <- str_remove(names(wafer_mda231p_50k_seg), '.sort')

wafer_mda231p_cor_1M <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_1M_seg)
wafer_mda231p_cor_750k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_750k_seg)
wafer_mda231p_cor_500k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_500k_seg)
wafer_mda231p_cor_250k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_250k_seg)
wafer_mda231p_cor_125k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_125k_seg)
wafer_mda231p_cor_75k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_75k_seg)
wafer_mda231p_cor_50k <- corr_cells(wafer_mda231p_original_seg, wafer_mda231p_50k_seg)

wafer_mda231p_corr_df <- create_corr_df(wafer_mda231p_cor_1M,
                                        wafer_mda231p_cor_750k,
                                        wafer_mda231p_cor_500k,
                                        wafer_mda231p_cor_250k,
                                        wafer_mda231p_cor_125k,
                                        wafer_mda231p_cor_75k,
                                        wafer_mda231p_cor_50k
)


p1 <- ggplot(wafer_mda231p_corr_df, aes(x = n_reads, y = correlation)) + 
  ylim(c(0, 1)) + 
  geom_jitter(aes(color = n_reads),position = position_jitter(seed = 1, width = 0.2), size = 0.01) +
  geom_violin(alpha=0, size = 0.2) + 
 stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
              fun.max = function(z) quantile(z, 0.75)) + 
cowplot::theme_cowplot() +
xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p1
cowplot::ggsave2("./figures/downsample_mda231p_spearman_cor_pvln.pdf", p1, width = 5, height = 4)  
write.csv(wafer_mda231p_corr_df, "./metrics/mda231_spearman_cor.csv")
# wafer_mda231p_corr_df <- read.csv("./metrics/mda231_spearman_cor.csv", row.names = 1)

# DUK249RE  (P1R: PCR duplicate rates: 22.76%)-----
# link_bam_files(bincounts = preduk249_bincounts,
#                path_bam_files = "/bam/file/path/",
#                output_path = "/downsample_corr/PREDUK249RE/",
#                output_name = "PREDUK249RE_cor_list",
#                target_reads = 1000000)

#./run_CNA.sh --bam preduk249_bam_input_1M --sample preduk249_1M --output preduk249_1M --cpu 10
#./run_CNA.sh --bam preduk249_bam_input_750k --sample preduk249_750k --output preduk249_750k --cpu 10
#./run_CNA.sh --bam preduk249_bam_input_500k --sample preduk249_500k --output preduk249_500k --cpu 10
#./run_CNA.sh --bam preduk249_bam_input_250k --sample preduk249_250k --output preduk249_250k --cpu 10
#./run_CNA.sh --bam preduk249_bam_input_125k --sample preduk249_125k --output preduk249_125k --cpu 10
#./run_CNA.sh --bam preduk249_bam_input_75k --sample preduk249_75k --output preduk249_75k --filter_ReadCount=10000 --cpu 10
#./run_CNA.sh --bam preduk249_bam_input_50k --filter_ReadCount=10000 --sample preduk249_50k --output preduk249_50k --cpu 10

#----Plot figures---
# reading files segs
preduk249_original_seg <- read.table("/CNA_pipeline/outpath/Duke27/final_result/uber.Duke27.seg.txt", header = T)
preduk249_1M_seg <- read.table("/CNA_pipeline/outpath/preduk249_1M/final_result/uber.preduk249_1M.seg.txt", header = T)
names(preduk249_1M_seg) <- str_remove(names(preduk249_1M_seg), '.sort.markdup')
preduk249_750k_seg <- read.table("/CNA_pipeline/outpath/preduk249_750k/final_result/uber.preduk249_750k.seg.txt", header = T)
names(preduk249_750k_seg) <- str_remove(names(preduk249_750k_seg), '.sort.markdup')
preduk249_500k_seg <- read.table("/CNA_pipeline/outpath/preduk249_500k/final_result/uber.preduk249_500k.seg.txt", header = T)
names(preduk249_500k_seg) <- str_remove(names(preduk249_500k_seg), '.sort.markdup')
preduk249_250k_seg <- read.table("/CNA_pipeline/outpath/preduk249_250k/final_result/uber.preduk249_250k.seg.txt", header = T)
names(preduk249_250k_seg) <- str_remove(names(preduk249_250k_seg), '.sort.markdup')
preduk249_125k_seg <- read.table("/CNA_pipeline/outpath/preduk249_125k/final_result/uber.preduk249_125k.seg.txt", header = T)
names(preduk249_125k_seg) <- str_remove(names(preduk249_125k_seg), '.sort.markdup')
preduk249_75k_seg <- read.table("/CNA_pipeline/outpath/preduk249_75k/final_result/uber.preduk249_75k.seg.txt", header = T)
names(preduk249_75k_seg) <- str_remove(names(preduk249_75k_seg), '.sort.markdup')
preduk249_50k_seg <- read.table("/CNA_pipeline/outpath/preduk249_50k/final_result/uber.preduk249_50k.seg.txt", header = T)
names(preduk249_50k_seg) <- str_remove(names(preduk249_50k_seg), '.sort.markdup')

preduk249_cor_1M <- corr_cells(preduk249_original_seg, preduk249_1M_seg)
preduk249_cor_750k <- corr_cells(preduk249_original_seg, preduk249_750k_seg)
preduk249_cor_500k <- corr_cells(preduk249_original_seg, preduk249_500k_seg)
preduk249_cor_250k <- corr_cells(preduk249_original_seg, preduk249_250k_seg)
preduk249_cor_125k <- corr_cells(preduk249_original_seg, preduk249_125k_seg)
preduk249_cor_75k <- corr_cells(preduk249_original_seg, preduk249_75k_seg)
preduk249_cor_50k <- corr_cells(preduk249_original_seg, preduk249_50k_seg)

preduk249_corr_df <- create_corr_df(preduk249_cor_1M,
                                    preduk249_cor_750k,
                                    preduk249_cor_500k,
                                    preduk249_cor_250k,
                                    preduk249_cor_125k,
                                    preduk249_cor_75k,
                                    preduk249_cor_50k)

p1 <- ggplot(preduk249_corr_df, aes(x = n_reads, y = correlation)) + 
  ylim(c(0, 1)) + 
  geom_jitter(aes(color = n_reads),position = position_jitter(seed = 1, width = 0.2), size = 0.01) +
  geom_violin(alpha=0, size = 0.2) + 
  stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75)) + 
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p1
cowplot::ggsave2("./figures/downsample_duke249_P1R_spearman_cor_pvln.pdf", p1, width = 5, height = 4)  

write.csv(preduk249_corr_df, "./metrics/duke249_P1R_spearman_cor.csv")
preduk249_corr_df <- read.csv("./metrics/duke249_P1R_spearman_cor.csv", row.names = 1)


# DUK254P  (P4P---PCR duplicate rates: 74.36%)-----
# link_bam_files(bincounts = preduk254_bincounts_p,
#                path_bam_files = "/bam/file/path/",
#                output_path = "/downsample_corr/PREDUK254P/",
#                output_name = "preduk254p_cor_list",
#                target_reads = 1000000)

#./run_CNA.sh --bam preduk254_bam_input_1M --sample preduk254_1M --output preduk254_1M --cpu 10
#./run_CNA.sh --bam preduk254_bam_input_750k --sample preduk254_750k --output preduk254_750k --cpu 10
#./run_CNA.sh --bam preduk254_bam_input_500k --sample preduk254_500k --output preduk254_500k --cpu 10
#./run_CNA.sh --bam preduk254_bam_input_250k --sample preduk254_250k --output preduk254_250k --cpu 10
#./run_CNA.sh --bam preduk254_bam_input_125k --sample preduk254_125k --output preduk254_125k --cpu 10
#./run_CNA.sh --bam preduk254_bam_input_75k --sample preduk254_75k --output preduk254_75k --cpu 10 --filter_ReadCount=10000 --filter_CellWithEmptyBin=1
#./run_CNA.sh --bam  preduk254_bam_input_50k --sample preduk254_50k --output preduk254_50k --cpu 10 --filter_ReadCount=10000 --filter_CellWithEmptyBin=1

#----Plot figures---
# reading files segs
preduk254_original_seg <- read.table("/CNA_pipeline/outpath/DUKE254_new_200/res_200_k/final_result/uber.sample.seg.txt", header = T)
preduk254_1M_seg <- read.table("/CNA_pipeline/outpath/preduk254_1M/final_result/uber.preduk254_1M.seg.txt", header = T)
names(preduk254_1M_seg) <- str_remove(names(preduk254_1M_seg), '.sort.markdup')
preduk254_750k_seg <- read.table("/CNA_pipeline/outpath/preduk254_750k/final_result/uber.preduk254_750k.seg.txt", header = T)
names(preduk254_750k_seg) <- str_remove(names(preduk254_750k_seg), '.sort.markdup')
preduk254_500k_seg <- read.table("/CNA_pipeline/outpath/preduk254_500k/final_result/uber.preduk254_500k.seg.txt", header = T)
names(preduk254_500k_seg) <- str_remove(names(preduk254_500k_seg), '.sort.markdup')
preduk254_250k_seg <- read.table("/CNA_pipeline/outpath/preduk254_250k/final_result/uber.preduk254_250k.seg.txt", header = T)
names(preduk254_250k_seg) <- str_remove(names(preduk254_250k_seg), '.sort.markdup')
preduk254_125k_seg <- read.table("/CNA_pipeline/outpath/preduk254_125k/final_result/uber.preduk254_125k.seg.txt", header = T)
names(preduk254_125k_seg) <- str_remove(names(preduk254_125k_seg), '.sort.markdup')
preduk254_75k_seg <- read.table("/CNA_pipeline/outpath/preduk254_75k/final_result/uber.preduk254_75k.seg.txt", header = T)
names(preduk254_75k_seg) <- str_remove(names(preduk254_75k_seg), '.sort.markdup')
preduk254_50k_seg <- read.table("/CNA_pipeline/outpath/preduk254_50k/final_result/uber.preduk254_50k.seg.txt", header = T)
names(preduk254_50k_seg) <- str_remove(names(preduk254_50k_seg), '.sort.markdup')

preduk254_cor_1M <- corr_cells(preduk254_original_seg, preduk254_1M_seg)
preduk254_cor_750k <- corr_cells(preduk254_original_seg, preduk254_750k_seg)
preduk254_cor_500k <- corr_cells(preduk254_original_seg, preduk254_500k_seg)
preduk254_cor_250k <- corr_cells(preduk254_original_seg, preduk254_250k_seg)
preduk254_cor_125k <- corr_cells(preduk254_original_seg, preduk254_125k_seg)
preduk254_cor_75k <- corr_cells(preduk254_original_seg, preduk254_75k_seg)
preduk254_cor_50k <- corr_cells(preduk254_original_seg, preduk254_50k_seg)

preduk254_corr_df <- create_corr_df(preduk254_cor_1M,
                                    preduk254_cor_750k,
                                    preduk254_cor_500k,
                                    preduk254_cor_250k,
                                    preduk254_cor_125k,
                                    preduk254_cor_75k,
                                    preduk254_cor_50k)

p1 <- ggplot(preduk254_corr_df, aes(x = n_reads, y = correlation)) + 
  ylim(c(0, 1)) + 
  geom_jitter(aes(color = n_reads),position = position_jitter(seed = 1, width = 0.2), size = 0.01) +
  geom_violin(alpha=0, size = 0.2) + 
  stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75)) + 
  cowplot::theme_cowplot() +
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p1
cowplot::ggsave2("./figures/downsample_duke254_P4P_spearman_cor_pvln.pdf", p1, width = 5, height = 4)  

write.csv(preduk254_corr_df, "./metrics/duke254_P4P_spearman_cor.csv")
# preduk254_corr_df <- read.csv("./metrics/duke254_P4P_spearman_cor.csv", row.names = 1)

#----Calculate median correlation score----
mda231 <- wafer_mda231p_corr_df %>% group_by(n_reads) %>% summarise(median_cor = median(correlation))
duke249 <- preduk249_corr_df %>% group_by(n_reads) %>% summarise(median_cor = median(correlation)) %>% pull(median_cor)
duke254 <- preduk254_corr_df %>% group_by(n_reads) %>% summarise(median_cor = median(correlation)) %>% pull(median_cor)
all_sum <- cbind(mda231,duke249, duke254) %>% column_to_rownames(var = "n_reads") 
all_sum <- all_sum[n_reads_levels,]

all_mean <- apply(all_sum, 1, mean)
all_sd <- apply(all_sum, 1, sd)

all_sum2 <- cbind(all_sum, all_mean, all_sd)

#----merged plots (Figure S4b)-----########
#--MDA231  (PCR duplicate rates: 16.43%)
wafer_mda231p_corr_df <- read.csv("./metrics/mda231_spearman_cor.csv", row.names = 1)
#--DUK249RE  (PCR duplicate rates: 22.76%)
preduk249_corr_df <- read.csv("./metrics/duke249_P1R_spearman_cor.csv", row.names = 1)
#--DUK254P  (PCR duplicate rates: 74.36%)
preduk254_corr_df <- read.csv("./metrics/duke254_P4P_spearman_cor.csv", row.names = 1)

wafer_mda231p_corr_df %>%  mutate(n_reads = fct_relevel(n_reads, n_reads_levels))
p1 <- wafer_mda231p_corr_df %>%  mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) %>% 
  ggplot(aes(x = n_reads, y = correlation)) + 
  ylim(c(0, 1)) + 
  geom_jitter(aes(color = n_reads),position = position_jitter(seed = 1, width = 0.2), size = 0.01) +
  geom_violin(alpha=0, size = 0.2) + 
  stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75)) + 
  cowplot::theme_cowplot() + ggtitle("MDA-MB-231 (PCR dup: 16.43%)") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = 'none', 
        plot.title = element_text(hjust = 0.5)) + 
  ylab("correlation with original profile") 

p2 <- preduk249_corr_df %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) %>% 
  ggplot(aes(x = n_reads, y = correlation)) + 
  ylim(c(0, 1)) + 
  geom_jitter(aes(color = n_reads),position = position_jitter(seed = 1, width = 0.2), size = 0.01) +
  geom_violin(alpha=0, size = 0.2) + 
  stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75)) + 
  cowplot::theme_cowplot() + ggtitle("P1R (PCR dup: 22.76%)") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = 'none', 
        plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()) 

p4 <- preduk254_corr_df %>% mutate(n_reads = fct_relevel(n_reads, n_reads_levels)) %>% 
  ggplot(aes(x = n_reads, y = correlation)) + 
  ylim(c(0, 1)) + 
  geom_jitter(aes(color = n_reads),position = position_jitter(seed = 1, width = 0.2), size = 0.01) +
  geom_violin(alpha=0, size = 0.2) + 
  stat_summary(fun = median, fun.min = function(z) quantile(z, 0.25),
               fun.max = function(z) quantile(z, 0.75)) + 
  cowplot::theme_cowplot() + ggtitle("P4P (PCR dup: 74.36%)") + 
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), axis.title.y = element_blank()) + 
  xlab("number of reads (downsampled)") + ylab("correlation with original profile")

p5 <- plot_grid(plotlist=list(p1,p2,p4), ncol=3, align='v')
p6 <- ggrastr::rasterize(p5, layers='Point', dpi=300)

cowplot::ggsave2("./figures/downsample_seg_correlation.pdf", p5, width = 10.5, height = 3)
cowplot::ggsave2("./figures/downsample_seg_correlation_raster.pdf", p6, width = 10.5, height = 3)


