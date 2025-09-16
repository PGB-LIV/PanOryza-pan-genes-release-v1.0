# Supplementary figure S12

library(tidyverse)
library(readxl)
library(cowplot)
library(data.table)
library(ggpubr)

load(file = "core_workspace.RData") #  load the core set of files with PanOryza results

#PanOryza results simplified
POR_long <- clusters_summary %>% select(Pan.id, taxa_size, transcript , genome) %>% distinct()

POR_long <- POR_long %>% group_by(Pan.id) %>% mutate(POR_size = n_distinct(transcript)) %>% ungroup()


# GENESPACE results 
# Load pangenes identified by GENESPACE
# '_pangenes.txt' file has pgID: the unique position-by-orthogroup combination identifier. 
# Each row in the pangenome-annotation is a unique value of this field. But it varies depending on the reference genome. Using Nipponbaremerged here.

os_nipmerged_pangenes_txt <- fread("core_files/oryza_sativa_nipponbaremerged_pangenes.txt")



## take only pgID and gene model ; remove not syntenic orthologs 

gp_long <- os_nipmerged_pangenes_txt %>%
  filter(!(flag == "NSOrtho")) %>%
  group_by(pgID) %>%
  mutate(GP_size = n_distinct(id) ,
         n_genomes = n_distinct(genome)) %>% 
  ungroup() %>%
  select(pgID , id , GP_size , n_genomes)


gp_long$pgID <- as.character(gp_long$pgID)




# Compare the two : Genespace vs GET_Pangenes 

# Find matching IDs

matching_ids <- inner_join(distinct(POR_long[,-4]), gp_long, join_by("transcript" == "id"))

all_ids <- full_join(distinct(POR_long[,c(1,3)]), gp_long, join_by("transcript" == "id"))

# Count matching IDs per group in both dataframes
matching_counts <- matching_ids %>%
  group_by(Pan.id, pgID) %>%
  summarise(matching_count = n())

# add cluster size info for both sets of pangenes

matching_counts <- matching_counts %>% left_join(distinct(POR_long[, c(1, 5)]) , by = "Pan.id")
matching_counts <- matching_counts %>% left_join(distinct(gp_long[, c(1, 3)]) , by = "pgID")

matching_counts <- matching_counts %>% left_join(distinct(POR_long[, c(1, 3)]) , by = "Pan.id" , relationship = "many-to-many")
matching_counts <- matching_counts %>% left_join(distinct(gp_long[, c(1, 2)]) , by = "pgID" , relationship = "many-to-many")


# smplified summary
matching_counts_summary <- matching_counts %>%
  group_by(Pan.id, pgID , matching_count) %>%
  pivot_longer(cols = c("transcript", "id"), names_to = "set_type", values_to = "all_transcripts") %>%
  summarise(POR_size = first(POR_size),
            GP_size = first(GP_size),
            total_count = n_distinct(all_transcripts)) %>%
  ungroup()

# keep only one row with max matches
matching_counts_summary <- matching_counts_summary %>% 
  group_by(Pan.id ) %>%
  filter(matching_count == max(matching_count)) %>%
  ungroup()

# Only unique values
matching_counts_summary <- matching_counts_summary %>% distinct()


matching_counts_summary <- matching_counts_summary %>% mutate(similarity = (matching_count/total_count)*100) %>% 
  arrange(similarity)

# Plot A
gp_por_density_plot <- ggplot(matching_counts_summary, aes(x = similarity )) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.5, fill = "darkgrey", color = "black") + theme_cowplot() + 
  labs(x = "Similarity in %", y = "Pangene counts", fill = "", title = "GENESPACE vs POR")  +
  scale_fill_manual(values = c("red", "blue" ))




# Compare the size of clusters in PanOryza vs GENESPACE for two families: NAC and NB-ARC

# read in NAC genes


NAC_total_m18_clusters <- read.delim("core_files/NAC_total_m18_clusters.tsv")

NAC_genes <- NAC_total_m18_clusters %>% select(gene_id) %>% distinct() # select only genes



# other family : NB-ARC

# read in files
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



# get pangenes for NAC and NB-ARC families

NAC_pangenes <- clusters_summary %>% 
  filter(gene_id %in% NAC_genes$gene_id) %>%
  select(Pan.id, gene_id, genome, taxa_size) %>%
  distinct()

nb_arc_pangenes <- clusters_summary %>% 
  filter(gene_id %in% nb_arc_m16$V1) %>%
  select(Pan.id, gene_id, genome, taxa_size) %>%
  distinct()



# filter matching_counts for NAC and NB-ARC families

matching_counts_NAC <- matching_counts %>% 
  filter(Pan.id %in% NAC_pangenes$Pan.id) %>%
  select(-transcript , -id) %>% 
  distinct()

matching_counts_nb_arc <- matching_counts %>% 
  filter(Pan.id %in% nb_arc_pangenes$Pan.id) %>%
  select(-transcript , -id) %>% 
  distinct()  


# get original sizes of clusters in GENESPACE and PanOryza for NAC and NB-ARC families

matching_counts_NAC_sum <- matching_counts_NAC %>% 
  mutate(group = "NAC") %>% select(group, GP_size, POR_size) %>%
  pivot_longer(cols = c(GP_size, POR_size), names_to = "type", values_to = "size") %>%
  distinct()


matching_counts_nb_arc_sum <- matching_counts_nb_arc %>%
  mutate(group = "NB-ARC") %>% select(group, GP_size, POR_size) %>%
  pivot_longer(cols = c(GP_size, POR_size), names_to = "type", values_to = "size") %>%
  distinct()

# recode GP_size as GP and POR_size as POR
matching_counts_NAC_sum$type <- recode(matching_counts_NAC_sum$type, 
                                       GP_size = "GS", 
                                       POR_size = "POR")
matching_counts_nb_arc_sum$type <- recode(matching_counts_nb_arc_sum$type, 
                                          GP_size = "GS", 
                                          POR_size = "POR")



#plot individually

cbPalette <- c( "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")

cluster_size_plot_NAC <- ggplot(matching_counts_NAC_sum, aes(x = type, y = size, fill = type)) +
  geom_boxplot(width = 0.4, alpha = 1, na.rm = TRUE , outliers = FALSE) +
  geom_violin(width=1, alpha = 0.3, color = NA) +  
  theme_cowplot(font_size = 14) +
  geom_jitter( alpha = 0.1 , width = 0.1) +
  labs(x = "NAC family",
       y = "Cluster Size" ,
       fill = "Resource type") +
  scale_fill_manual(values=cbPalette) +
  theme(legend.position = "top")

# plot NB-ARC


cluster_size_plot_nbarc <- ggplot(matching_counts_nb_arc_sum, aes(x = type, y = size, fill = type)) +
  geom_boxplot(width = 0.4, alpha = 1, na.rm = TRUE , outliers = FALSE) +
  geom_violin(width=1, alpha = 0.3, color = NA) +  
  geom_jitter( alpha = 0.1 , width = 0.1) +
  theme_cowplot(font_size = 14) +
  labs(x = "NB-ARC family", 
       y = "Cluster Size" ,
       fill = "Resource type") +
  scale_fill_manual(values=cbPalette) +
  theme(legend.position = "top")


plots_grid <- ggarrange(
  cluster_size_plot_NAC , cluster_size_plot_nbarc , 
  labels = c("B", "C") ,
  ncol=2,
  nrow=1)


plots_grid1 <- ggarrange(
  gp_por_density_plot , plots_grid , 
  labels = c("A") ,
  widths = c(2, 1),
  ncol=1,
  nrow=2)



# Uncomment below to export the figure

# svg(filename = "figure_S.svg", width=7,height=8.5, pointsize = 9)
# plots_grid1
# dev.off()



# stats for genespace vs POR


matching_counts_NAC_sum %>% group_by(type) %>% summarise(median_val = median(size),
                                                         IQR = IQR(size))

matching_counts_NAC_sum %>% filter(size == 1) %>% group_by(type) %>% summarise(counts = n_distinct(Pan.id))


