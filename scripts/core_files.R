## core files

library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(progress)
library(reshape2)
library(data.table)


magic18_interpro_clusters <- fread("core_files/magic18_interpro_updated_clusters.tab")

#magic18_interpro <- fread("core_files/magic18_interpro_updated.tab")

magic18_models <- fread("core_files/magic18_models.tab")

clusters_summary_final <- read.delim("core_files/clusters_summary_final.tab")

cluster_merged <- read.delim2("core_files/cluster_merged.tab")

#clusters_dir <- "core_files/Os4530.POR.1" 

pangene_list <- read.delim2("core_files/pangene_matrix_genes.table")

df_merged <- read.delim("core_files/df_merged.tab", sep = "\t")
df_merged$taxa_size <- as.numeric(df_merged$taxa_size)
df_merged$Cluster_size <- as.numeric(df_merged$Cluster_size)

gff_data_subset <- read.delim("core_files/nipponbaremerged_gff_subset.tab") 


clusters_summary <- clusters_summary_final %>% left_join(distinct(cluster_merged[,c(1,3,4)]), join_by(Pan.id == oryzasativanipponbaremerged))

clusters_summary$occupancy = with(clusters_summary, ifelse(taxa_size == "16", "core",
                                                           ifelse(taxa_size == "15", "softcore",
                                                                  ifelse(taxa_size == "1", "cloud", 
                                                                         ifelse(taxa_size == "2", "cloud" , "shell")))))

clusters_summary <- clusters_summary %>% arrange(., Pan.id) 
clusters_summary <- clusters_summary %>% select(Pan.id, Cluster_size, taxa_size, occupancy, transcript, merged.or.gene_id ,region, Length, genome , source, File)

clusters_summary <- clusters_summary %>% left_join(magic18_models , join_by("transcript" == "transcipt")) %>% select(-species) %>% distinct()

clusters_summary$accession <- clusters_summary$genome

clusters_summary <- clusters_summary %>% 
  mutate(., accession = recode(accession, "Nipponbare(RAPDB)" = "Nipponbare", "Nipponbare(MSU)" = "Nipponbare" , "Nipponbare(Gramene)" = "Nipponbare"))


magic18_interpro_clusters$accession <- magic18_interpro_clusters$genome

magic18_interpro_clusters <- magic18_interpro_clusters %>% 
  mutate(., accession = recode(accession, "Nipponbare(RAPDB)" = "Nipponbare", "Nipponbare(MSU)" = "Nipponbare" , "Nipponbare(Gramene)" = "Nipponbare"))

