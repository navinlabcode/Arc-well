# show colors -------
library(scales)
library(egg)
library(grid)
library(RColorBrewer)
library(tidyft)
library(tidyverse)
library(Rphenograph)
library(dplyr)
library(dbscan)
library(dendextend)
library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
require(RColorBrewer)
library(fastcluster)
library(ape)
library(ggtree)

#show_col(col)
ernst = c("#91323A", "#3A4960", "#6D7345", "#554540", "#D7C969")
escher = c("#C1395E", "#AEC17B", "#E07B42", "#89A7C2", "#F0CA50")
col_major_clones = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#FED439B2", "#709AE1B2", "#8A9197B2",
                     "#D2AF81B2", "#FD7446B2", "#D5E4A2B2", "#197EC0B2", "#F05C3BB2", "#46732EB2")
new_pal = c("#CC0C00B2", "#5C88DAB2", "#84BD00B2", "#FFCD00B2", "#7C878EB2", "#00B5E2B2", "#00AF66B2", "#D2AF81B2", "#FD7446B2", "#46732EB2", "#C1395E", "#E07B42",
            "#D4A2D9", "#8E72D5","#C0EDB9", "#364E4F", "#8EE5EE", "#FFA500", "#458B00", "#CD6090", "#FFAEB9", "#90EE90", "#5f9EA0",
            "#E6E6FA", "#8B7E66", ernst)
colr = new_pal
colourCount = 50
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
colors_set1 = getPalette(12)
getPalette3 = colorRampPalette(brewer.pal(9, "Set3"))
colors_set3 = getPalette3(9);show_col(c(colors_set1, colors_set3))

col_minor_clones = c(colors_set3,colors_set1)

col_grp = c("Primary" = "deeppink", "Recurrence" = "chartreuse1") ##6A6599  
col_peak = c("a_well" = "grey", "d_well" = "seagreen", "pop_well" = "black")
