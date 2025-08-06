# Supplementary figure S3
# Pangene positions 

#Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(readxl))
suppressMessages(library(data.table))


# Read the data , also available as supplementary table S3
df <- fread("./pangenes_locations_with_IRGSP.tsv", sep="\t", na.strings="None") # Set the correct path to the file


# use the IRGSP position only i.e. those who have an equivalent in IRGSP get placed not those who are missing on it



pan_genes_subset_chr <- df %>% select(IRGSP_position, oryza_sativa_lima_mode_position, oryza_sativa_ir64_mode_position,
                                      oryza_sativa_khaoyaiguang_mode_position, oryza_sativa_larhamugad_mode_position, oryza_sativa_gobolsailbalam_mode_position,
                                      oryza_sativa_liuxu_mode_position, oryza_sativa_mh63_mode_position, oryza_sativa_ZS97_mode_position, oryza_sativa_pr106_mode_position,
                                      oryza_sativa_natelboro_mode_position, oryza_sativa_n22_mode_position, oryza_sativa_arc_mode_position, oryza_sativa_azucena_mode_position,
                                      oryza_sativa_ketannangka_mode_position, oryza_sativa_chaomeo_mode_position, oryza_sativa_nipponbaremerged_mode_position, mode_chr)

colnames(pan_genes_subset_chr) <- c("IRGSP_position", "Lima", "IR64", "KYG", "Lamu", "Gosa", 
                                    "Lixu", "MH63", "ZS97", "PR106", "Nabo", "N22", "ARC", "AZU", "Kena", "Chao", "Nippon", "mode_chr")

#Bring into longer format

melted_chr <- pan_genes_subset_chr %>% 
  pivot_longer(cols = -c(IRGSP_position, mode_chr), names_to = "genome", values_to = "mode_position")



#set colors

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique(melted_chr$genome))

#plot

# uncomment below to export image
#png("./june_2025/pan_gene_dotplots.png", height  = 8000, width = 8000, res = 600)

options(scipen = 999)

chr_dotplot <- ggplot(melted_chr) +
  geom_point(aes(x = IRGSP_position, y = mode_position, color = genome), size = 0.02)  + 
  theme_cowplot(font_size = 10) +
  scale_colour_manual(values = getPalette(colourCount)) +
  labs( color="accession", y="Mode position", x="IRGSP Position") + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  facet_wrap(~ mode_chr , nrow = 4, ncol = 3, scales = "free", switch = "x" )+
  theme(strip.text = element_text(size = 12))

chr_dotplot

#dev.off()

