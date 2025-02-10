# Figure S6

# method: to compare both the set of clusters, first prepared a list of clusters in both sets that contain RAPDB ids. 
# Then filtered the full set of clusters only for pan-genes that contain a RAPDB id
# then removed the incomparable ids from both set of pangene clusters. So effectively comparing without 3 sets of genome annotations. 


suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
suppressMessages(library(cowplot))

load(file = "core_workspace.RData") #  load the core set of files


# Determine the maximum number of columns for RGI upload
max_cols <- max(count.fields("core_files/RGI_clusters_sorted.txt", sep = "\t"))
RGI <- read.table("core_files/RGI_clusters_sorted.txt", sep = "\t", fill = TRUE, header = FALSE, stringsAsFactors = FALSE , col.names = paste0("V", 1:max_cols))

# PanOryza clusters
POR_long <- clusters_summary %>% select(Pan.id, taxa_size, gene_id , genome) %>% distinct()

POR_long <- POR_long %>% group_by(Pan.id) %>% mutate(POR_size = n_distinct(gene_id)) %>% ungroup()

POR_long_rapdb <- POR_long %>% filter(genome == "Nipponbare(RAPDB)") %>% select(Pan.id) %>% distinct() # pangenes with only rapdb ids


# Compare
RGI_long <- as.data.frame(pivot_longer(RGI, cols = -V1, names_to = "id_col", values_to = "id")) %>% select(-id_col)

RGI_long <- RGI_long %>% filter_all(all_vars(. != ""))

RGI_long <- RGI_long %>% distinct()

RGI_clusters_occupancy <- RGI_long %>% group_by(V1) %>% summarise(RGI_size = n_distinct(id))



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


matching_counts <- matching_counts %>% # keep only one row with max matches
  group_by(Pan.id ) %>%
  filter(matching_count == max(matching_count)) %>%
  ungroup()

# for all ids against which similarity would be calculated 

df1_list <- list()
df2_list <- list()
combined_list <- list()

for (i in 1:nrow(matching_counts)){
  
  x <- matching_counts$Pan.id[i]
  y <- matching_counts$V1[i]
  
  df1_list[matching_counts$Pan.id[i]] = POR_long_filtered1 %>% filter(Pan.id == matching_counts$Pan.id[i]) %>% select(gene_id)
  
  df2_list[matching_counts$V1[i]] = RGI_long_filtered1 %>% filter(V1 == matching_counts$V1[i]) %>% select(id)
  
  combined_list[[x]] <- unique(c(df1_list[[x]], df2_list[[y]]))
}



combined_df <- do.call(rbind, lapply(names(combined_list), function(name) {
  data.frame(
    pan_gene = name,
    gene_id = combined_list[[name]],
    stringsAsFactors = FALSE
  )
}))

combined_df1 <- combined_df %>%
  group_by(pan_gene) %>%
  summarise(Total_count = n_distinct(gene_id))


matching_counts1 <- matching_counts %>% left_join(combined_df1, join_by(Pan.id == pan_gene))

matching_counts2 <- matching_counts1 %>% mutate(similarity = (matching_count/Total_count)*100) # calculate similarity based on counts

matching_counts2 <- matching_counts2 %>% arrange(similarity)

# plot

rgi_por_rapdb_density_plot <- ggplot(matching_counts2, aes(x = similarity )) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.5, fill = "darkgrey", color = "black") + theme_cowplot() + 
  labs(x = "Similarity in %", y = "Pangene counts", fill = "", title = "RGI vs POR")  +
  scale_fill_manual(values = c("red", "blue" ))

