
# Supplementary Figure S2

# Load Libraries
suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(readxl))


pan_genes_test_impute <- read.csv("core_files/pan_genes_test_impute.csv") 

# select only the mode positions of all accessions and pan-gene pseudo positions

pan_genes_subset_chr <- pan_genes_test_impute %>% select(pseudo_position, oryza_sativa_lima_mode_position, oryza_sativa_ir64_mode_position,
                                                         oryza_sativa_khaoyaiguang_mode_position, oryza_sativa_larhamugad_mode_position, oryza_sativa_gobolsailbalam_mode_position,
                                                         oryza_sativa_liuxu_mode_position, oryza_sativa_mh63_mode_position, oryza_sativa_ZS97_mode_position, oryza_sativa_pr106_mode_position,
                                                         oryza_sativa_natelboro_mode_position, oryza_sativa_n22_mode_position, oryza_sativa_arc_mode_position, oryza_sativa_azucena_mode_position,
                                                         oryza_sativa_ketannangka_mode_position, oryza_sativa_chaomeo_mode_position, oryza_sativa_nipponbaremerged_mode_position, mode_chr)

colnames(pan_genes_subset_chr) <- c("pseudo_position", "Lima", "IR64", "KYG", "Lamu", "Gosa", 
                                    "Lixu", "MH63", "ZS97", "PR106", "Nabo", "N22", "ARC", "AZU", "Kena", "Chao", "Nippon", "mode_chr")

#Bring into longer format
melted_chr <- reshape2::melt(pan_genes_subset_chr, id.vars = c("pseudo_position", "mode_chr"))

#set colors
melted_chr1 <-melted_chr[melted_chr$mode_chr == "1",]
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(melted_chr1$variable))

#plot

# uncomment below to export image
#png("./pan_gene_dotplots_070225.png", height  = 8000, width = 8000, res = 600)

options(scipen = 999)

chr_dotplot <- ggplot(melted_chr) +
  geom_point(aes(x = pseudo_position, y = value, color = variable), size = 0.02)  + 
  theme_cowplot(font_size = 10) +
  scale_colour_manual(values = getPalette(colourCount)) +
  labs( color="accession", y="Mode position", x="Pseudo Position") + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  facet_wrap(~ mode_chr , nrow = 4, ncol = 3, scales = "free", switch = "x" )+
  theme(strip.text = element_text(size = 12))
chr_dotplot
#dev.off() 
