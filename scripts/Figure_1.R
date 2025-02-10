suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(gridExtra))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(patchwork))
suppressMessages(library(ggplotify))
suppressMessages(library(ggbreak))
suppressMessages(library(forcats))


setwd("<DIR_PATH>")
list.files()
###load files
load(file = "core_workspace.RData") #  load the core set of files
clusters_dir <- "./Os4530.POR.1" # directory with version 1 of pan-genes

###colors
color_scheme <- read.delim("core_files/color_scheme.txt")

col_palette <- color_scheme$color

names(col_palette) <- color_scheme$genome

col_palette2 <- color_scheme$color2

names(col_palette2) <- color_scheme$subspecies  

#Plot 1a
# exclude MSU singletons that do not cluster with rest of the pan-genes

plot1_data <- clusters_summary %>% filter(!(taxa_size == "1" & genome == "Nipponbare(MSU)"))

plot1_data <- plot1_data %>% group_by(taxa_size) %>% summarise(counts = n_distinct(Pan.id))


plot1_data$occupancy = with(plot1_data, ifelse(taxa_size == "16", "core",
                                               ifelse(taxa_size == "15", "softcore",
                                                      ifelse(taxa_size == "1", "cloud", 
                                                             ifelse(taxa_size == "2", "cloud" , "shell")))))



plot1a <- ggplot(plot1_data, aes(x = as.factor(taxa_size), y = counts, fill = occupancy)) + 
  geom_col()  +
  theme_cowplot(font_size = 14 , font_family  ="arial") +
  labs(fill="Occupancy", y="Number of pan-gene clusters", x="Pan-gene occupancy")

#optional, uncomment below to output image
# svg("/mnt/hc-storage/users/esharma/rice/ms/plots/plot1a.svg", width=5.5,height=3.5)
# plot1a
# dev.off()

# Plot 1b 

# created using GET_PANGENES pipeline (https://github.com/Ensembl/plant-scripts/tree/master/pangenes)

#Plot 1c


df_singletons <- clusters_summary[clusters_summary$taxa_size == "1",] %>% select(-File, -region, -Length)

df_singletons <- df_singletons %>% filter(!(genome %in% c("Nipponbare(MSU)", "Nipponbare(RAPDB)")))

singletons_sum <- df_singletons %>% group_by(genome) %>% summarise(counts = n_distinct(Pan.id)) %>% arrange(desc(counts))

singletons_sum <- singletons_sum %>% left_join(color_scheme[,c(1,4)], by = "genome")

col_palette_rev <- c("#2E7D32" , "#FF822E", "#0003C7", "#374140" )
names(col_palette_rev) <- c("aromatic" ,"aus", "indica" , "japonica" )


singletons_plot <- ggplot(singletons_sum) + geom_col(aes(y = reorder(genome, counts), x = counts, fill = subspecies) ) + 
  theme_cowplot(font_size = 14) +  scale_fill_manual(values = col_palette_rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1)) + 
  xlab("Number of pan-genes") + ylab(" Genome") + labs(title = "singletons", fill = "Subpopulation") 


#optional, uncomment to output image
# svg("/mnt/hc-storage/users/esharma/rice/ms/plots/singletons_plot.svg")
# singletons_plot
# dev.off()



# Plot 1d

# 0.1 call packages silently
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library("ape"))

# 0.2 set options
options(expressions = 100000) #https://stat.ethz.ch/pipermail/r-help/2004-January/044109.html
opar<-par(no.readonly = TRUE)

# 1. read and process the table
tab <- read.table(file=paste0(clusters_dir, "/POCS.matrix.tab"), header=TRUE, sep = "\t")
tab <- tab[complete.cases(tab), ]

## use accession names as row and column  names
ena_species <- read.delim("core_files/ena_species.txt")

ena_species <- bind_rows(ena_species, data.frame(species_ena = "Oryza_sativa_nipponbaremerged", species ="Nipponbare(merged)"))


ena_species$species_ena <- str_replace_all(ena_species$species_ena, "^O", "o")
ena_species$species_ena <- str_replace_all(ena_species$species_ena, "liuxu", "liuxu_chr")
ena_species$species_ena <- str_replace_all(ena_species$species_ena, "zs97", "ZS97")

ena_species <- ena_species[-6,]
  
for (i in 1:nrow(tab)) {
  search_keyword <- ena_species[i,1]
  replace_keyword <- ena_species[i,2]
  tab$genomes <- str_replace_all(tab$genomes, search_keyword, replace_keyword)
}

colnames(tab) <- c("genomes", tab$genomes)

#tab <- na.omit(tab)

if( dim(tab)[1] < 5) stop('There are less than four complete data rows. Pleae revise your input table!') 

dfr.num <- tab[,2:ncol(tab)]
dfr.num <- droplevels.data.frame(dfr.num)

#  convert dfr.num to matrix; this one will not be rounded, 
#    as it will be used to compute ANDg and the
#    mean silhouette width statistic
dfr.num.mat <- as.matrix(dfr.num)

# this matrix will be rounded to plot as heatmaps
mat_dat <- data.matrix(tab[,2:ncol(tab)])
mat_dat <- round(mat_dat,0)

rnames <- tab[,1]
rownames(mat_dat) <- rnames


# 1.2 subset mat_dat and dfr.num with user-provided regex, if requested
if(0 > 0 ){
  include_list <- grep("", rownames(mat_dat))
  mat_dat <- mat_dat[include_list, ]
  mat_dat <- mat_dat[,include_list]
  mat_dat <- round(mat_dat,0)
  
  #subset the dfr.num
  dfr.num <- dfr.num[include_list,]
  dfr.num <- dfr.num[,include_list]
  
  # convert to matrix
  dfr.num.mat <- as.matrix(dfr.num)
  rownames(dfr.num.mat) <- names(dfr.num)
  colnames(dfr.num.mat) <- names(dfr.num)
}

# 1.2 subset mat_dat to avoid excessive redundancy
if(100 < 100)
{
  tmp_mat = mat_dat
  diag(tmp_mat) = NA
  rows <- (!apply( tmp_mat , 1 , function(x) any( x > 100 , na.rm=T) ) )
  mat_dat <- mat_dat[rows, rows]
}


sidecolors <- c( "Nipponbare(merged)" = "#374140" ,
                 "Khao Yai Guang" = "#0003C7",
                 "Larha Mugad" = "#0003C7",
                 "Minghui 63" = "#0003C7" ,
                 "Chao Meo" = "#374140",
                 "Azucena" = "#374140" ,
                 "N22" = "#FF822E",
                 "Natel Boro" = "#FF822E",
                 "IR64" = "#0003C7" ,
                 "ARC 10497" = "#2E7D32",
                 "Zhenshan 97" = "#0003C7",
                 "PR106" = "#0003C7",
                 "Lima" = "#0003C7",
                 "Gobol Sail Balam" = "#0003C7",
                 "Ketan Nangka"= "#374140",
                 "Liu Xu" = "#0003C7")



if(1 > 0){
  if(0 > 0){
    svg("POCS.matrix_heatmap.svg", width=15, height=10, pointsize=)  
    heatmap.2(mat_dat, main="POCS.matrix.tab", notecol="black", density.info="none", key.xlab="Percent Conserved Sequences (POCS)",
              trace="none", margins=c(18,18), lhei = c(1,5), labCol=F,
              cexRow=1.0, cexCol=1.0, srtRow=45, srtCol=45)
    dev.off()
  }
  else {
    svg("POCS.matrix_heatmap.svg", width=15, height=10, pointsize=)  
    heatmap.2(mat_dat, cellnote=mat_dat, main="POCS.matrix.tab", notecol="black", density.info="none", key.xlab="Percent Conserved Sequences (POCS)",
              trace="none", margins=c(18,18), lhei = c(1,5), 
              cexRow=2, cexCol=2, srtRow=45, srtCol=45, colRow = sidecolors, colCol = sidecolors, RowSideColors = sidecolors, ColSideColors = sidecolors, key = TRUE)
    dev.off()
  }
} else {
  if(0 > 0){
    svg("POCS.matrix_heatmap.svg", width=15, height=10, pointsize=)  
    heatmap.2(mat_dat, main="POCS.matrix.tab", notecol="black", density.info="none", labCol=F, key.xlab="Percent Conserved Sequences (POCS)",
              trace="none", margins=c(18,18), lhei = c(1,5), dendrogram = "row", Colv = FALSE, 
              cexRow=1.0, cexCol=1.0, srtRow=45, srtCol=45)
    dev.off()
  }
  else {
    svg("POCS.matrix_heatmap.svg", width=15, height=10, pointsize=)  
    heatmap.2(mat_dat, cellnote=mat_dat, main="POCS.matrix.tab", notecol="black", density.info="none", key.xlab="Percent Conserved Sequences (POCS)",
              trace="none", margins=c(18,18), lhei = c(1,5), dendrogram = "row", Colv = FALSE,
              cexRow=1.0, cexCol=1.0, srtRow=45, srtCol=45)
    dev.off()
  }
}  



