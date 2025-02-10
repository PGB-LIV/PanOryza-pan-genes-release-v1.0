chooseCRANmirror(ind = 1)
setRepositories(ind = 1:2)

# List of packages you want to check and install 
packages <- c('tidyverse', 'readxl', 'InteractiveComplexHeatmap' , 'shiny', 'shinyjs', 'ComplexHeatmap') 

# Loop through each package 
for (pkg in packages) { 
  if (!requireNamespace(pkg, quietly = TRUE)) { 
    install.packages(pkg) } 
  }

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(InteractiveComplexHeatmap))


pan_genes_test_impute <- suppressMessages(read_excel("data/pan_genes_positions.xlsx")) ## latest version used 



pan_genes_subset_number <- pan_genes_test_impute %>% select(`Pan-gene identifier`,
                                                            oryza_sativa_nipponbaremerged_chr,
                                                            oryza_sativa_chaomeo_chr,
                                                            oryza_sativa_ketannangka_chr,
                                                            oryza_sativa_azucena_chr,
                                                            oryza_sativa_arc_chr,
                                                            oryza_sativa_n22_chr,
                                                            oryza_sativa_natelboro_chr,
                                                            oryza_sativa_pr106_chr,
                                                            oryza_sativa_ZS97_chr,
                                                            oryza_sativa_mh63_chr,
                                                            oryza_sativa_liuxu_chr,
                                                            oryza_sativa_gobolsailbalam_chr, 
                                                            oryza_sativa_larhamugad_chr,
                                                            oryza_sativa_khaoyaiguang_chr,
                                                            oryza_sativa_ir64_chr,
                                                            oryza_sativa_lima_chr, mode_chr, occupancy, pseudo_position)


pan_genes_subset_number <- pan_genes_subset_number[pan_genes_subset_number$occupancy > 2,] ## remove cloud occupancy genes
pan_genes_subset_number <- as.data.frame(pan_genes_subset_number)
rownames(pan_genes_subset_number) <- pan_genes_subset_number$`Pan-gene identifier`

pan_gene_annot <- data.frame(row.names = pan_genes_subset_number$`Pan-gene identifier`, 
                             occupancy = factor(pan_genes_subset_number$occupancy, levels = 1:16))



color_scheme <- read.delim("data/color_scheme.txt")

color_scheme$sample1 <- c("oryza_sativa_arc_chr", 
                          "oryza_sativa_n22_chr",
                          "oryza_sativa_natelboro_chr",
                          "oryza_sativa_khaoyaiguang_chr",
                          "oryza_sativa_larhamugad_chr",
                          "oryza_sativa_mh63_chr",
                          "oryza_sativa_ir64_chr",
                          "oryza_sativa_ZS97_chr",
                          "oryza_sativa_pr106_chr",
                          "oryza_sativa_lima_chr",
                          "oryza_sativa_gobolsailbalam_chr",
                          "oryza_sativa_liuxu_chr",
                          "oryza_sativa_nipponbaremerged_chr",
                          "oryza_sativa_chaomeo_chr",
                          "oryza_sativa_azucena_chr",
                          "oryza_sativa_ketannangka_chr",
                          "oryza_sativa_nipponbaremerged_chr",
                          "oryza_sativa_nipponbaremerged_chr")

sample_annot <- color_scheme %>% select(sample, subspecies) %>% distinct()

row.names(sample_annot) <- sample_annot$sample 
sample_annot <- sample_annot %>% select(-sample)

pan_genes_subset_number <- pan_genes_subset_number %>% 
  select(-`Pan-gene identifier`) 

my_col = c("#33608c", "#406588", "#4d6984", "#5a6e80", "#66737d", "#737779", "#807c75", ### cols for occupancy
           "#8d8171", "#9a856d", "#a78a69", "#b48f65", "#c19361", "#cd985e",
           "#da9d5a", "#e7a156", "#f4a652")

names(my_col) = levels(pan_gene_annot$occupancy)

my_col1 <- list(occupancy = my_col)

# my_col2 <- unique(color_scheme$color2)
# my_col2 <- my_col2[-4]
# 
# names(my_col2) = c("aromatic", "aus", "indica", "japonica")

col_palette_rev <- c("#2E7D32" , "#FF822E", "#0003C7", "#374140" )
names(col_palette_rev) <- c("aromatic" ,"aus", "indica" , "japonica" )



annot_colors <- list(occupancy = my_col, subspecies = col_palette_rev)

hap_cols <- c("#EEDFCC" , "#8B8378", "#0000FF", "#FF4040", "#7AC5CD", "#66CD00", "#8B3E2F", "#EE7600", "#BF3EFF", "#FF1493", "#EEC900","#CDC673")


###
old_colnames <- names(pan_genes_subset_number)[-17:-19]
new_colnames <- color_scheme$sample[match(old_colnames, color_scheme$sample1)]

names(pan_genes_subset_number) <- c(new_colnames, "mode_chr", "occupancy", "pseudo_position")

###########################

chr_annot <- list()
chromo = list()
plot_by_chr = list()

##########
suppressWarnings(for (chr in unique(pan_genes_subset_number$mode_chr)){
  chromo[[chr]] <- pan_genes_subset_number[pan_genes_subset_number$mode_chr == chr,]
  
  # chromo[[chr]] <- chromo[[chr]][order(row.names(chromo[[chr]])), ]
  chromo[[chr]] <- chromo[[chr]][order(chromo[[chr]]$pseudo_position), ]
  
  chr_annot[[chr]] <- data.frame(row.names = row.names(chromo[[chr]]), 
                                 occupancy = chromo[[chr]]$occupancy)
  
  chromo[[chr]] <- chromo[[chr]] %>% select(-occupancy, -pseudo_position, -mode_chr)
  chromo[[chr]] <- as.matrix(chromo[[chr]])
  
  plot_by_chr[[chr]] <- ComplexHeatmap::pheatmap(chromo[[chr]],na_col = "white", cluster_rows = FALSE, cluster_cols = FALSE , show_rownames	 = FALSE,
                                                 color = hap_cols, annotation_col = sample_annot, show_colnames = TRUE,
                                                 main = paste("chromosome", chr), annotation_row = chr_annot[[chr]],
                                                 annotation_colors = annot_colors, legend_breaks = 1:12, 
                                                 legend = FALSE, annotation_legend = FALSE, annotation_names_row = FALSE,
                                                 annotation_names_col = FALSE, use_raster = FALSE)
  
})
