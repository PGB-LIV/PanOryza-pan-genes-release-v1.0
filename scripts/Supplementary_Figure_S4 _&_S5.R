# Figure S4 and S5

suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(ggh4x))


load(file = "core_workspace.RData") #  load the core set of files

clusters_dir <- "./Os4530.POR.1" # directory with version 1 of pan-genes

color_scheme <- read.delim("core_files/color_scheme.txt")

suppressMessages(source("scripts/add_consensus.R"))


# Fig S4
# 1. NAC family

NAC_total_m18_clusters <- read.delim("core_files/NAC_total_m18_clusters.tsv")

NAC_total_m18_clusters <- NAC_total_m18_clusters %>% dplyr::rename(genome = "species")

NAC_subset <- NAC_total_m18_clusters %>% select(Pan.id, merged.or.gene_id, genome, taxa_size, accession) %>% distinct() # using merged gene id for nipponbare, same for rest

clusters_summary_NAC <- clusters_summary %>% 
  filter(Pan.id %in% NAC_subset$Pan.id) %>% 
  select(Pan.id , Cluster_size, taxa_size , merged.or.gene_id, genome, gene_id, accession)

clusters_summary_NAC <- clusters_summary_NAC %>% 
  mutate(., accession = recode(accession, "Nipponbare" = "Nipponbare(merged)"))
         
clusters_summary_NAC <- clusters_summary_NAC %>% mutate(isNAC = ifelse(merged.or.gene_id %in% NAC_subset$merged.or.gene_id , "yes", "no"))


NAC_final_table <- add_consensus(clusters_summary_NAC, "Pan.id", "isNAC") # add consensus column


NAC_final_table <- NAC_final_table %>% left_join(color_scheme[, c(1,4)], join_by(genome == genome)) # add varietal group information


# plot the heatmap

color5 <- c("#2E7D32","#FF822E","#0003C7","#D90000","#374140")
names(color5) <- c("aromatic" , "aus", "indica", "japonica", "Nipponbare(merged)")


color10 <- c("#0000FF","#0beffe","#ff6961","#ff69b4","#800080","#ffa500","#ff4500","#ffff00","#004225","#89C35C")
strip <- strip_themed(background_x = elem_list_rect(fill = color5))

cols = c("darkgreen", "lightgreen", "white")
names(cols) = c("yes", "no" , NA)


p_1_15 <- ggplot(data = NAC_final_table %>% filter(taxa_size <16), aes(x = accession, y = reorder(Pan.id,dplyr::desc(taxa_size)), fill = consensus)) +
  geom_tile(color = "black", height = 1) +
  scale_fill_manual(values = cols) +
  ylab("Pangene IDs")+
  xlab("Rice") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(color = "white", size = 0)) +
  facet_grid2(rows = vars(taxa_size), cols = vars(subspecies), scales = "free", space = "free", strip = strip)

p_16 <- ggplot(data = NAC_final_table %>% filter(taxa_size >15), aes(x = accession, y = reorder(Pan.id,dplyr::desc(taxa_size)), fill = consensus)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = cols)+
  ylab("Pangene IDs")+
  xlab("Rice") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.x = element_text(color = "white", size = 0),
        axis.text.y = element_text(size = 4)) +
  facet_grid2(rows = vars(taxa_size), cols = vars(subspecies), scales = "free", space = "free", strip = strip)


ggarrange(p_1_15, p_16, ncol = 2, common.legend = T)
# uncomment below to export image
# png(filename = "NAC_summary.png", height = 9, width = 9, units = "in", res = 300)
# ggarrange(p_1_15, p_16, ncol = 2, common.legend = T)
# dev.off()



# Supplementary Figure S5

# other family : NB-ARC

#read in files
nb_arc_dir <- ("core_files/oyrza_magic16_NB_ARC")
nb_arc_files <- list.files("core_files/oyrza_magic16_NB_ARC", pattern = "ARC.list")


nb_arc_list <- list()

for (i in 1:length(nb_arc_files)) {
  nb_arc_list[[i]] <- read.delim(paste0(nb_arc_dir, "/", nb_arc_files[i]), header = FALSE)
}

# merge all lists of individual files
nb_arc_m16 <- bind_rows(nb_arc_list)

nb_arc_m16 <- distinct(nb_arc_m16)

length(unique(nb_arc_m16$V1)) # 8401 NB-ARC genes 

nb_arc_m16$V1 <- str_replace_all(nb_arc_m16$V1, "gene:", "") # correction for Osnip models

# get the clusters for nb arc genes

nb_arc_clusters <- magic18_interpro_clusters %>% filter(gene_id %in% nb_arc_m16$V1)

nb_arc_cluster_subset <- nb_arc_clusters %>% select(Pan.id, original_transcript, gene_id, genome , merged.or.gene_id , taxa_size, accession) %>% distinct()


clusters_summary_NBARC <- clusters_summary %>% 
  filter(Pan.id %in% nb_arc_cluster_subset$Pan.id) %>% 
  select(Pan.id , Cluster_size, taxa_size , merged.or.gene_id, genome, gene_id, accession)

clusters_summary_NBARC <- clusters_summary_NBARC %>% 
  mutate(., accession = recode(accession, "Nipponbare" = "Nipponbare(merged)"))

clusters_summary_NBARC <- clusters_summary_NBARC %>% mutate(isNBARC = ifelse(merged.or.gene_id %in% nb_arc_cluster_subset$merged.or.gene_id , "yes", "no"))

NBARC_final_table <- add_consensus(clusters_summary_NBARC, "Pan.id", "isNBARC") # add consensus column


NBARC_final_table <- NBARC_final_table %>% left_join(color_scheme[, c(1,4)], join_by(genome == genome)) # add varietal group information


# plot heatmap

color5 <- c("#2E7D32","#FF822E","#0003C7","#D90000","#374140")
names(color5) <- c("aromatic" , "aus", "indica", "japonica", "Nipponbare(merged)")


color10 <- c("#0000FF","#0beffe","#ff6961","#ff69b4","#800080","#ffa500","#ff4500","#ffff00","#004225","#89C35C")
strip <- strip_themed(background_x = elem_list_rect(fill = color5))

cols = c("darkgreen", "lightgreen", "white")
names(cols) = c("yes", "no" , NA)

nb_p_cloud <- ggplot(data = NBARC_final_table %>% filter(taxa_size < 3 & taxa_size > 0 ), aes(x = accession, y = reorder(Pan.id,dplyr::desc(taxa_size)), fill = consensus)) +
  geom_tile(color = "black", height = 1) +
  scale_fill_manual(values = cols) +
  ylab("Pangene IDs")+
  xlab("Rice") +
  theme_bw() +
  labs(title = "Cloud") +
  theme(plot.title = element_text(size = 14, hjust = 0.5), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 2),
        strip.text.x = element_text(color = "white", size = 0)) +
  facet_grid2(rows = vars(taxa_size), cols = vars(subspecies), scales = "free", space = "free", strip = strip)

nb_p_shell <- ggplot(data = NBARC_final_table %>% filter(taxa_size > 2 & taxa_size < 15), aes(x = accession, y = reorder(Pan.id,dplyr::desc(taxa_size)), fill = consensus)) +
  geom_tile(color = "black", height = 1) +
  scale_fill_manual(values = cols) +
  ylab("Pangene IDs")+
  xlab("Rice") +
  theme_bw() +
  labs(title = "Shell") +
  theme(plot.title = element_text(size = 14, hjust = 0.5), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 2),
        strip.text.x = element_text(color = "white", size = 0)) +
  facet_grid2(rows = vars(taxa_size), cols = vars(subspecies), scales = "free", space = "free", strip = strip)

nb_p_core <- ggplot(data = NBARC_final_table %>% filter(taxa_size >=15), aes(x = accession, y = reorder(Pan.id,dplyr::desc(taxa_size)), fill = consensus)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = cols)+
  ylab("Pangene IDs")+
  xlab("Rice") +
  theme_bw() +
  labs(title = "Softcore and Core") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text.x = element_text(color = "white", size = 0), 
    axis.text.y = element_text(size = 2)) +
  facet_grid2(rows = vars(taxa_size), cols = vars(subspecies), scales = "free", space = "free", strip = strip)

ggarrange(nb_p_cloud, nb_p_shell, nb_p_core, ncol = 3, common.legend = T)

# uncomment below to export image
# png(filename = "NB_ARC_summary.png", height = 9, width = 12, units = "in", res = 300)
# ggarrange(nb_p_cloud, nb_p_shell, nb_p_core, ncol = 3, common.legend = T)
# dev.off()


