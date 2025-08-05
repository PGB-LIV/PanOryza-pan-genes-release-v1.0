# Supplementary figure 
# RGI vs PanOryza

# method: to compare both the set of clusters, first prepared a list of clusters in both sets that contain RAPDB ids. 
# Then filtered the full set of clusters only for pan-genes that contain a RAPDB id
# then removed the incomparable ids from both set of pangene clusters. So effectively comparing without 3 sets of genome annotations. 


suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))

load(file = "core_workspace.RData") #  load the core set of files


# Determine the maximum number of columns for RGI upload
max_cols <- max(count.fields("core_files/RGI_clusters_sorted.txt", sep = "\t"))
RGI <- read.table("core_files/RGI_clusters_sorted.txt", sep = "\t", fill = TRUE, header = FALSE, stringsAsFactors = FALSE , col.names = paste0("V", 1:max_cols))

# PanOryza clusters
POR_long <- clusters_summary %>% select(Pan.id, taxa_size, gene_id , genome) %>% distinct()

POR_long <- POR_long %>% group_by(Pan.id) %>% mutate(original_POR_size = n_distinct(gene_id)) %>% ungroup()

POR_long_rapdb <- POR_long %>% filter(genome == "Nipponbare(RAPDB)") %>% select(Pan.id) %>% distinct() # pangenes with only rapdb ids


# Compare
RGI_long <- as.data.frame(pivot_longer(RGI, cols = -V1, names_to = "id_col", values_to = "id")) %>% select(-id_col)

RGI_long <- RGI_long %>% filter_all(all_vars(. != ""))

RGI_long <- RGI_long %>% distinct()

RGI_clusters_occupancy <- RGI_long %>% group_by(V1) %>% summarise(original_RGI_size = n_distinct(id))



# take only rapdb based clusters


RGI_long_rapdb <- magic18_models %>% filter(species == "Nipponbare(RAPDB)" ) %>% inner_join(RGI_long, join_by(gene_id == id)) %>% select(V1) %>% distinct()


### filter the clusters for only those clusters that contain a RAPDB gene


POR_long_filtered <- POR_long %>% filter(Pan.id %in% POR_long_rapdb$Pan.id)

##remove ids not compatible match in both i.e mh63 , zs97 and osnip since the ids are following a different rule in RGI clusters , Number of characters in gene ids are also different. 


POR_long_filtered1 <- POR_long_filtered %>% filter(!(grepl("^OsMH63|^OsZS97|^OsNip" , gene_id)))

RGI_long_filtered <- RGI_long %>% filter(V1 %in% RGI_long_rapdb$V1)

RGI_long_filtered1 <- RGI_long_filtered %>% filter(!(grepl("^OsMH63|^OsZS97|^OsNip" , id)))



# Find matching IDs


matching_ids <- inner_join(distinct(POR_long_filtered1[,c(1,3)]), RGI_long_filtered1, join_by("gene_id" == "id"))

all_ids <- full_join(distinct(POR_long_filtered1[,c(1,3)]), RGI_long_filtered1, join_by("gene_id" == "id"))

# Count matching IDs per group in both dataframes
matching_counts <- matching_ids %>%
  group_by(Pan.id, V1) %>%
  summarise(matching_count = n())


matching_counts <- matching_counts %>% left_join(RGI_clusters_occupancy, by = "V1")
matching_counts <- matching_counts %>% left_join(distinct(POR_long[,c(1,5)]), by = "Pan.id")



# corrected way 

matching_counts <- matching_counts %>% left_join(distinct(POR_long_filtered1[, c(1, 3)]) , by = "Pan.id" , relationship = "many-to-many")
matching_counts <- matching_counts %>% left_join(distinct(RGI_long_filtered1) , by = "V1" , relationship = "many-to-many")

matching_counts_summary <- matching_counts %>%
  group_by(Pan.id, V1 , matching_count) %>%
  pivot_longer(cols = c("gene_id", "id"), names_to = "set_type", values_to = "all_genes") %>%
  summarise(POR_size = first(original_POR_size),
            RGI_size = first(original_RGI_size),
            total_count = n_distinct(all_genes)) %>%
  ungroup()

matching_counts_summary <- matching_counts_summary %>% # keep only one row with max matches
  group_by(Pan.id ) %>%
  filter(matching_count == max(matching_count)) %>%
  ungroup()

matching_counts_summary <- matching_counts_summary %>% distinct()


matching_counts_summary <- matching_counts_summary %>% mutate(similarity = (matching_count/total_count)*100) %>% 
  arrange(similarity)

# Plot A
rgi_por_rapdb_density_plot <- ggplot(matching_counts_summary, aes(x = similarity )) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.5, fill = "darkgrey", color = "black") + theme_cowplot() + 
  labs(x = "Similarity in %", y = "Pangene counts", fill = "", title = "RGI vs POR")  +
  scale_fill_manual(values = c("red", "blue" ))



# Compare the size of clusters in PanOryza vs RGI for two families: NAC and NB-ARC

# read in NAC genes


NAC_total_m18_clusters <- read.delim("core_files/NAC_total_m18_clusters.tsv")

NAC_genes <- NAC_total_m18_clusters %>% select(gene_id) %>% distinct() ## select only genes



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


### use matching_counts to filter for NAC and NB-ARC genes as they are filtered out for three genomes with non-matching ids

# get pangenes for NAC and NB-ARC families

NAC_pangenes <- clusters_summary %>% 
  filter(gene_id %in% NAC_genes$gene_id) %>%
  select(Pan.id, gene_id, genome, taxa_size) %>%
  distinct()

nb_arc_pangenes <- clusters_summary %>% 
  filter(gene_id %in% nb_arc_m16$V1) %>%
  select(Pan.id, gene_id, genome, taxa_size) %>%
  distinct()


head(matching_counts)

# filter for NAC and NB-ARC families

matching_counts_NAC <- matching_counts %>% 
  filter(Pan.id %in% NAC_pangenes$Pan.id) %>%
  distinct()

matching_counts_nb_arc <- matching_counts %>% 
  filter(Pan.id %in% nb_arc_pangenes$Pan.id) %>%
  distinct()  

# get original sizes of clusters in RGI and PanOryza for NAC and NB-ARC families

matching_counts_NAC_sum <- matching_counts_NAC %>% 
  mutate(group = "NAC") %>% select(group, original_RGI_size, original_POR_size) %>%
  pivot_longer(cols = c(original_RGI_size, original_POR_size), names_to = "type", values_to = "size") %>%
  distinct()


matching_counts_nb_arc_sum <- matching_counts_nb_arc %>%
  mutate(group = "NB-ARC") %>% select(group, original_RGI_size, original_POR_size) %>%
  pivot_longer(cols = c(original_RGI_size, original_POR_size), names_to = "type", values_to = "size") %>%
  distinct()

# recode original_RGI_size as RGI and original_POR_size as POR
matching_counts_NAC_sum$type <- recode(matching_counts_NAC_sum$type, 
                                        original_RGI_size = "RGI", 
                                        original_POR_size = "POR")
matching_counts_nb_arc_sum$type <- recode(matching_counts_nb_arc_sum$type, 
                                           original_RGI_size = "RGI", 
                                           original_POR_size = "POR")



#plot individually B and C

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
  theme(legend.position = "top") + ylim(0,250)


plots_grid <- ggarrange(
  cluster_size_plot_NAC , cluster_size_plot_nbarc , 
  labels = c("B", "C") ,
  ncol=2,
  nrow=1)


plots_grid1 <- ggarrange(
  rgi_por_rapdb_density_plot , plots_grid , 
  labels = c("A") ,
  widths = c(2, 1),
  ncol=1,
  nrow=2)

#  Uncomment to export the figure
# svg(filename = "RGI_VS_POR_figure_S.svg", width=7,height=8.5, pointsize = 9)
# plots_grid1
# dev.off()

