## Figure 4

suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(progress))
suppressMessages(library(cowplot))
suppressMessages(library(purrr))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(ggrepel))


load(file = "core_workspace.RData") #  load the core set of files
clusters_dir <- "./Os4530.POR.1" # directory with version 1 of pan-genes

# get all pfam and interpro domains in the clusters
pfam_annot <- magic18_interpro_clusters %>% filter(V4 == "Pfam") %>% select(V5, V6) %>% distinct() 
ipr_annot <- magic18_interpro_clusters %>% select(V12, V13) %>% distinct()


# To perform various analyses on interpro and pfam annotated pan-genes

#either load analysed data as a list of data frames object.
load("interpro_pfam_check.RData")

#or to calculate fresh using # external script (takes about 10-15 minutes on our systems) # uncomment to run
# source("./scripts/summarise_interpro_pfam_tables.R")
# updated_list <- summarise_interpro_pfam(magic18_interpro_clusters) 


## separate the results and view the tables

gene_count_table_interpro <- updated_list$gene_count_table 

easy_interpro_summary <- updated_list$easy_interpro_summary

merged_summary_interpro <- updated_list$merged_summary_interpro_pfam

merged_summary_interpro2 <- as_tibble(updated_list$merged_summary_interpro_extended)

gene_count_table_pfam <- updated_list$gene_count_table 

merged_summary_pfam2 <- updated_list$merged_summary_pfam_extended

# no. of pan-genes
merged_summary_pfam2 %>% filter(!(pangene_occupancy == 1 )) %>% filter(accessions_lacking_pfam_domain == 0) %>% summarise(pan_gene_counts = length(unique(Pan.id)))



## Analysis of those lacking domain
ALL <- merged_summary_pfam2 %>% filter(pangene_occupancy == 16) %>% filter(accessions_lacking_pfam_domain == 0) %>% select(accesion_names_with_common_pfam_domain) %>% distinct()
ALL <- ALL$accesion_names_with_common_pfam_domain[1]
ALL <- str_replace_all(unlist(str_split(ALL, ",")) , "^ " , "")

not_consistent_pfam <- merged_summary_pfam2  %>% filter(!(pangene_occupancy == 1 )) %>% filter(!(accessions_lacking_pfam_domain == 0) )

not_consistent_interpro <- merged_summary_interpro2 %>% filter(!(pangene_occupancy == 1 )) %>% filter(!(accessions_lacking_interpro_domain == 0) )

find_missing_names <- function(names_in_A, all_names) {
  names_in_A <- unlist(str_split(names_in_A, ","))
  names_in_A <- str_replace_all(names_in_A, "^ " , "")
  missing_names <- setdiff(all_names, names_in_A)
  return(paste(missing_names, collapse = ","))
}

not_consistent_pfam <- not_consistent_pfam %>%
  rowwise() %>%
  mutate(accession_names_lacking_pfam_domain = find_missing_names(accesion_names_with_common_pfam_domain, ALL))

not_consistent_pfam_genomes <- not_consistent_pfam %>% select(Pan.id, accession_names_lacking_pfam_domain) %>% separate_wider_delim(., cols = accession_names_lacking_pfam_domain, delim = "," , names_sep = "" , too_few = "align_start") %>% 
  pivot_longer(cols = starts_with("accession_names"))

not_consistent_pfam_genomes %>% na.omit() %>% group_by(value) %>% summarise(counts = length(unique(Pan.id)))

not_consistent_interpro <- not_consistent_interpro %>%
  rowwise() %>%
  mutate(accession_names_lacking_interpro_domain = find_missing_names(accesion_names_with_common_interpro_domain, ALL))

not_consistent_interpro_genomes <- not_consistent_interpro %>% select(Pan.id, accession_names_lacking_interpro_domain) %>% separate_wider_delim(., cols = accession_names_lacking_interpro_domain, delim = "," , names_sep = "" , too_few = "align_start") %>% 
  pivot_longer(cols = starts_with("accession_names"))

not_consistent_interpro_genomes %>% na.omit() %>% group_by(value) %>% summarise(counts = length(unique(Pan.id)))


magic18_interpro_clusters %>% filter(!(V5 == "-" | V12 == "-")) %>% group_by(genome) %>% summarise(gene_counts = n_distinct(gene_id)) ## matches total gene counts with gene counts table


#Fig 4a
# plot 1 core vs cloud pangenes interpro significant terms


# remove TE and plot core vs cloud terms

msu_locus_data <- read.delim("core_files/all.locus_brief_info.7.0")

MSU_TE <- msu_locus_data %>% filter(., is_TE == "Y" )  ## TE from MSU annotations

iproscan_TE_magic18 <- magic18_interpro_clusters %>% filter(grepl('transposon',V6) | grepl('transposon',V13))  ## interproscan based results , 4271 clusters with TE in them

# 73259 clusters after filtering TE based clusters

clusters_summary_TE_excl <- clusters_summary %>% filter(!(Pan.id %in% unique(iproscan_TE_magic18$Pan.id)))  ## TE excluded based on interpro annotations

# if taxa size=1, remove those clusters where transcript labelled as TE in msu list or interproscan

clusters_summary_size1 <- clusters_summary[clusters_summary$taxa_size == "1",] 

filtered_pan_ids_msu <- clusters_summary_size1 %>%
  group_by(Pan.id) %>%
  summarize(any_TE = any(transcript %in% MSU_TE$model)) %>%
  filter(!any_TE) %>%
  select(Pan.id)

filtered_pan_ids_ipr <- clusters_summary_size1 %>%
  group_by(Pan.id) %>%
  summarize(any_TE_ipr = any(transcript %in% iproscan_TE_magic18$original_transcript)) %>%
  filter(!any_TE_ipr) %>%
  select(Pan.id)

filtered_pan_ids <- bind_rows(filtered_pan_ids_msu, filtered_pan_ids_ipr ) %>% distinct()

TE_fil_clusters_summary_size1 <- clusters_summary %>% #### taxa size 1 clusters with TE excluded
  semi_join(filtered_pan_ids, by = "Pan.id")


# rest taxa sizes
# if taxa size >1, remove only those clusters where all transcripts labelled as TE in interproscan

clusters_summary_rest <- clusters_summary[!(clusters_summary$taxa_size == "1"),] 

filtered_pan_ids_ipr_all <- clusters_summary_rest %>%  ### is MSU transcript TE
  group_by(Pan.id, transcript, genome) %>%
  filter(genome == "Nipponbare(MSU)") %>% 
  summarize(is_MSU_TE  = ifelse(transcript %in% iproscan_TE_magic18$original_transcript | transcript %in% MSU_TE$model, "yes", "no")) %>% ungroup()


filtered_pan_ids_ipr_all1 <- clusters_summary_rest %>%  # if any other transcript other than MSU is TE ? 
  group_by(Pan.id, transcript, genome) %>%
  filter(!(genome == "Nipponbare(MSU)")) %>%
  summarize(is_TE = ifelse(transcript %in% iproscan_TE_magic18$original_transcript | transcript %in% MSU_TE$model, "yes", "no")) %>% ungroup()


filtered_pan_ids_ipr_all_merged <- bind_rows(filtered_pan_ids_ipr_all , filtered_pan_ids_ipr_all1)

df_list <- split(filtered_pan_ids_ipr_all_merged, filtered_pan_ids_ipr_all_merged$Pan.id) # split into list of pan-gene groups

df_list1 <- lapply(df_list, function(x) { # if both MSU transcript and  transcripts in other genomes other than MSU are TE , then label the cluster to remove. 
  x %>% group_by(Pan.id) %>% summarise(consensus = ifelse(any(is_MSU_TE == "yes") & any(is_TE == "yes") , "remove", "keep"))
})

pan_ids_TE_consensus <- bind_rows(df_list1) %>% filter(consensus == "remove")

TE_fil_clusters_summary_rest <- clusters_summary_rest %>% filter(!(Pan.id %in% pan_ids_TE_consensus$Pan.id))  ## remove consensus "remove" pan-gene ids

TE_fil_clusters_summary <- bind_rows(TE_fil_clusters_summary_rest, TE_fil_clusters_summary_size1) %>% distinct() ### join taxa size 1 and the rest data

magic18_interpro_clusters_TE_fil1 <- magic18_interpro_clusters %>% filter(Pan.id %in% TE_fil_clusters_summary$Pan.id)

ipr_clusters_TE_fil <- magic18_interpro_clusters_TE_fil1 %>% 
  select(Pan.id, V12) %>% 
  filter(!(V12 == "-")) %>% 
  na.omit() %>%
  distinct()


ipr_clusters_TE_fil <- ipr_clusters_TE_fil %>% left_join(distinct(magic18_interpro_clusters_TE_fil1[,c(1,11)]), by = "Pan.id")

ipr_clusters_TE_fil$occupancy <- with(ipr_clusters_TE_fil, ifelse(taxa_size == "16", "core",
                                                                          ifelse(taxa_size == "15", "softcore",
                                                                                 ifelse(taxa_size == "1", "cloud",
                                                                                        ifelse(taxa_size == "2", "cloud" , "shell")))))

# Separate out clusters based on occupancy

core_ipr_clusters <- ipr_clusters_TE_fil[ipr_clusters_TE_fil$occupancy == "core",]


cloud_ipr_clusters <- ipr_clusters_TE_fil[ipr_clusters_TE_fil$occupancy == "cloud",]


softcore_ipr_clusters <- ipr_clusters_TE_fil[ipr_clusters_TE_fil$occupancy == "softcore",]


shell_ipr_clusters <- ipr_clusters_TE_fil[ipr_clusters_TE_fil$occupancy == "shell",]


ipr_clusters <- magic18_interpro_clusters %>% 
  select(Pan.id, V12) %>% 
  filter(!(V12 == "-")) %>% 
  na.omit() %>%
  distinct()


ipr_annot <- magic18_interpro_clusters %>% select(V12, V13) %>% distinct()

## prepare for clusterprofiler 

term2name <- ipr_annot %>% filter(grepl('^IPR', V12))
term2gene <- ipr_clusters %>% select(V12 , Pan.id)

genelist_core = unique(core_ipr_clusters$Pan.id)
genelist_cloud = unique(cloud_ipr_clusters$Pan.id)
genelist_shell = unique(shell_ipr_clusters$Pan.id)
genelist_softcore = unique(softcore_ipr_clusters$Pan.id)

library(clusterProfiler)

core_test <- enricher(genelist_core, TERM2GENE = term2gene , TERM2NAME = term2name , pvalueCutoff = 0.05, pAdjustMethod = "BH")

cloud_test <- enricher(genelist_cloud, TERM2GENE = term2gene , TERM2NAME = term2name , pvalueCutoff = 0.05, pAdjustMethod = "BH")


options(scipen = 0)

core_res <- core_test@result %>% select(-geneID)

cloud_res <- cloud_test@result %>% select(-geneID)


core_res <- core_res %>% mutate(gene_ratio_num = Count / 20913)

cloud_res <- cloud_res %>% mutate(gene_ratio_num = Count / 13649)


core_res_30 <- core_res %>% arrange(p.adjust) %>% slice(1:30) %>% mutate(group  = "core")

cloud_res_30 <- cloud_res %>% arrange(p.adjust) %>% slice(1:30) %>% mutate(group  = "cloud")


# collapse redundant terms using text 

library(dplyr)
library(stringdist)


# Function to find near duplicates and keep only one entry
find_and_filter_duplicates <- function(df, text_column, value_column, threshold = 0.15) {
  duplicates <- data.frame(index1 = integer(), index2 = integer(), distance = numeric())
  for (i in 1:(nrow(df) - 1)) {
    for (j in (i + 1):nrow(df)) {
      distance <- stringdist::stringdist(df[i, text_column], df[j, text_column], method = "jw")
      if (distance <= threshold) {
        duplicates <- rbind(duplicates, data.frame(index1 = i, index2 = j, distance = distance))
      }
    }
  }
  
  # Keep only one entry based on the value column
  to_remove <- c()
  for (row in 1:nrow(duplicates)) {
    index1 <- duplicates[row, "index1"]
    index2 <- duplicates[row, "index2"]
    if (df[index1, value_column] >= df[index2, value_column]) {
      to_remove <- c(to_remove, index2)
    } else {
      to_remove <- c(to_remove, index1)
    }
  }
  
  df_filtered <- df[-unique(to_remove), ]
  return(df_filtered)
}

# Identify and filter near duplicates

core_res_30_f <- find_and_filter_duplicates(core_res_30, "Description", "gene_ratio_num")

cloud_res_30_f <- find_and_filter_duplicates(cloud_res_30, "Description", "gene_ratio_num")


combined_30_f <- bind_rows(core_res_30_f,cloud_res_30_f ) # combine cloud and core

#Plot 4a
top30_dotplot_v <- ggplot(combined_30_f, aes(x= group, y= reorder(Description,p.adjust) , size= gene_ratio_num, color=p.adjust)) +
  geom_point()   +
  scale_color_viridis_c(name = 'p-adjusted', direction =1, option = "turbo") + 
  theme_bw() + theme(axis.line  = element_blank()) +
  labs(x = "Pan-gene occupancy", y = " Interpro terms", size = "Gene ratio") +
  theme(axis.text.y = element_text(size = 10,family  = "arial"), axis.text.x = element_text( vjust = 0.5, hjust=1)) +
  theme(axis.ticks = element_blank()) 

#Uncomment to export image
# svg("./plots/interpro_dotplot_corevscloud_filtered30_repeat.svg", height = 13 , width = 9)
# top30_dotplot_v
# dev.off()
# 


# Fig 4b
# calculate domain occupancy ; excl MSU &  RAPDB

magic16_interpro_clusters <- magic18_interpro_clusters %>% 
  filter(!(genome %in% c("Nipponbare(MSU)" , "Nipponbare(RAPDB)"))) 


domain_occupancy_table_pfam <- magic16_interpro_clusters %>% 
  filter(V4 == "Pfam") %>% select(Pan.id, V5, genome) %>% distinct() %>%
  group_by(V5, Pan.id ) %>%
  summarize(domain_occupancy = length(genome))  

mean_occupancy_pfam <- domain_occupancy_table_pfam %>%
  group_by(V5) %>%
  summarise(count_pangenes = length(Pan.id) , 
    mean_occupancy  =  mean(domain_occupancy)) %>% left_join(pfam_annot, by = "V5")


domain_occupancy_table_interpro <- magic16_interpro_clusters %>% 
  filter(!(V12 == "-")) %>%
  select(Pan.id, V12, genome) %>% distinct() %>%
  group_by(V12, Pan.id ) %>%
  summarize(domain_occupancy = length(genome))   


mean_occupancy_interpro <- domain_occupancy_table_interpro %>%
  group_by(V12) %>%
  summarise(count_pangenes = length(Pan.id),
            mean_occupancy  =  mean(domain_occupancy)) %>% left_join(ipr_annot, by = "V12")


# optional
# write.table(domain_occupancy_table_pfam, file = "domain_occupancy_table_pfam.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(domain_occupancy_table_interpro, file = "domain_occupancy_table_interpro.tsv", sep = "\t", quote = FALSE, row.names = FALSE)  


mean_occupancy_pfam <- mean_occupancy_pfam %>% rename(id = "V5")
mean_occupancy_interpro <- mean_occupancy_interpro %>% rename(id = "V12")

mean_occupancy_pfam$type = "Pfam"
mean_occupancy_interpro$type = "Interpro"

combined_domain_occupancy <-  bind_rows(mean_occupancy_interpro, mean_occupancy_pfam)

# filter for presence of domain in at least 5 pan-genes
combined_domain_occupancy_fil <- combined_domain_occupancy %>% filter(count_pangenes > 5)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")


# plot Fig 4b
filtered_domain_occupancy_histogram <- ggplot(combined_domain_occupancy_fil, aes(x = mean_occupancy , fill = type)) +
  geom_histogram(position = "identity", alpha = 0.5) + theme_cowplot() + labs(x = "Mean occupancy", y = "Domain counts", fill = "", subtitle = "pangene counts > 5") +
  scale_fill_manual(values = c("red", "blue" ))


# add_dotted_line
filtered_domain_occupancy_histogram <- filtered_domain_occupancy_histogram + geom_vline(xintercept = c(10, 15), linetype = "dotted", color = "black")

filtered_domain_occupancy_histogram <- filtered_domain_occupancy_histogram + 
  annotate("text", x = c(5,12.5,16), y = 330, label = c("hv" , "pv","inv" ) , colour = "black")


#uncomment to export image
# svg("./plots/filtered_domain_occupancy_histogram.svg")
# filtered_domain_occupancy_histogram
# dev.off()


# Fig 4c

# 1. NAC family

NAC_total_m18_clusters <- read.delim("core_files/NAC_total_m18_clusters.tsv")

NAC_genes <- NAC_total_m18_clusters %>% select(gene_id) %>% distinct() ## select only genes

# subset by NAC genes
NAC_clusters <- magic18_interpro_clusters %>% filter(gene_id %in% NAC_genes$gene_id)

# shorten df
NAC_cluster_subset <- NAC_clusters %>% select(Pan.id, original_transcript, gene_id, genome , merged.or.gene_id , taxa_size, accession) %>% distinct()

NAC_cluster_subset$occupancy = with(NAC_cluster_subset, ifelse(taxa_size == "16", "core", 
                                                               ifelse(taxa_size == "15", "softcore",
                                                                      ifelse(taxa_size == "1", "cloud",
                                                                             ifelse(taxa_size == "2", "cloud" , "shell")))))

# counts of genes and pan-genes
TF_number_summary <- NAC_cluster_subset %>%
  group_by(accession) %>%
  summarize(gene_count = length(unique(merged.or.gene_id)),
            pan_gene_count = length(unique(Pan.id))) %>%
  melt(id.vars = "accession") %>% na.omit()

TF_number_summary1 <- NAC_cluster_subset %>%
  group_by(occupancy, accession) %>%
  summarize(pan_gene_count = length(unique(Pan.id)))

nac_occupancy_status <- NAC_cluster_subset %>%
  group_by(taxa_size, occupancy) %>%
  summarize(pan_gene_count = length(unique(Pan.id)))


cbPalette <- c( "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7")

plot_NAC_occupancy <- ggplot(nac_occupancy_status, aes(x=as.factor(taxa_size), y = pan_gene_count, fill=occupancy)) +
  geom_col( width = 0.4)+
  labs(title = "NAC family",fill="occupancy", y="Number of clusters", x="Pan-gene occupancy") +
  theme_cowplot(font_size = 14) 


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


nb_arc_cluster_subset$occupancy = with(nb_arc_cluster_subset, ifelse(taxa_size == "16", "core", 
                                                                     ifelse(taxa_size == "15", "softcore",
                                                                            ifelse(taxa_size == "1", "cloud",
                                                                                   ifelse(taxa_size == "2", "cloud" , "shell")))))

# counts of genes and pan-genes
nb_arc_number_summary <- nb_arc_cluster_subset %>%
  group_by(accession) %>%
  summarize(gene_count = length(unique(merged.or.gene_id)),
            pan_gene_count = length(unique(Pan.id))) %>%
  melt(id.vars = "accession") %>% na.omit()

nbarc_occupancy_status <- nb_arc_cluster_subset %>%
  group_by(taxa_size, occupancy) %>%
  summarize(pan_gene_count = length(unique(Pan.id)))


plot_NBARC_occupancy <- ggplot(nbarc_occupancy_status, aes(x=as.factor(taxa_size), y = pan_gene_count, fill=occupancy)) +
  geom_col( width = 0.4)+
  labs(title = "NB-ARC family",fill="occupancy", y="Number of clusters", x="Pan-gene occupancy") +
  theme_cowplot(font_size = 14) 



# compare both families together


combined_counts <- bind_rows(TF_number_summary, nb_arc_number_summary, .id = "Family")

combined_counts$Family<- str_replace_all(combined_counts$Family, c("1" = "NAC" , "2" = "NB-ARC")) 

combined_gene_counts <- combined_counts[combined_counts$variable == "gene_count",]

plot_combined_counts <- ggplot(combined_gene_counts, aes(x=accession, y = value, fill=Family)) +
  geom_bar(width = 0.8, stat = "identity", position = position_dodge(0.9)) +
  labs( fill="Gene family", y="Number of genes", x="Source") +
  theme_cowplot(font_size = 12) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  scale_fill_manual(labels = c("NAC", "NB-ARC"), values=cbPalette)


# since numbers are very different range, calculate percentage out of total numbers of clusters

fam_proportion_nac <- nac_occupancy_status %>% 
  mutate(proportion = pan_gene_count / sum(nac_occupancy_status$pan_gene_count))

fam_proportion_nb_arc <- nbarc_occupancy_status %>%
  mutate(proportion = pan_gene_count / sum(nbarc_occupancy_status$pan_gene_count))

fam_proportion_nac %>% group_by(occupancy) %>% summarise(proportions = sum(proportion))
fam_proportion_nb_arc %>% group_by(occupancy) %>% summarise(proportions = sum(proportion))

fam_prop_cb <- bind_rows(fam_proportion_nac,fam_proportion_nb_arc , .id = "Family")
fam_prop_cb$Family<- str_replace_all(fam_prop_cb$Family, c("1" = "NAC" , "2" = "NB-ARC")) 

fam_prop_cb$occupancy = with(fam_prop_cb, ifelse(taxa_size == "16", "core",
                                                 ifelse(taxa_size == "15", "softcore",
                                                        ifelse(taxa_size == "1", "cloud", 
                                                               ifelse(taxa_size == "2", "cloud" , "shell")))))

# plot Fig 4c

proportion_fam_plot <- ggplot(fam_prop_cb, aes(x=as.factor(taxa_size), y = proportion *100, fill = Family)) +
  geom_bar(stat = "identity",  position = position_dodge(1))+ scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  labs(fill="Gene family", y="Proportion of clusters", x="Pan-gene occupancy") +
  theme_cowplot(font_size = 12) 


# Fig 4d
# violin plot using interpro domains

magic18_interpro <- fread("core_files/magic18_interpro_updated.tab")

# TF first

TF_table <- read.delim("core_files/TF (1).list", header = FALSE)

TFdb_interpro <- TF_table %>% left_join(magic18_interpro, join_by(V1 == V1)) 

counts_TF <- TFdb_interpro %>% group_by(V2.x) %>% summarise(transcript_counts = length(unique(V1)))

TFdb_interpro_list <- split(TFdb_interpro, TFdb_interpro$V2.x)

find_common_ipr_value <- function(x){
  df <- x %>% group_by(V2.x) %>% filter(!(V12 == "-")) %>% summarise(most_common_value = names(which.max(table(V12)))) %>% na.omit()
  return(df)
}

# essential interpro domain for each TF family

essential_TF_ipr <- lapply(TFdb_interpro_list , FUN = find_common_ipr_value) %>% 
  bind_rows() %>% 
  rename(TF = "V2.x" , interpro_id = 'most_common_value')


# all interpro domains that are present in about 80% of the genes of each TF family

table_all_TF_interpro <- lapply(TFdb_interpro_list , FUN = function(x){
  x <- x %>% select(V1, V12) %>% distinct() %>% filter(!(V12 == "-")) 
  x <- table(x$V12)
  x <- as.data.frame(x)
  x <- x %>% filter(Freq >= 0.8*max(Freq))
})


df_all_TF_interpro <- bind_rows(table_all_TF_interpro, .id = "id")

ipr_annot <- magic18_interpro %>% select(V12, V13) %>% distinct()


essential_TF_ipr <- essential_TF_ipr %>% left_join(ipr_annot, join_by(interpro_id == V12))
df_all_TF_interpro <- df_all_TF_interpro %>% left_join(ipr_annot, join_by(Var1 == V12))

# filter for more than 5 pan-genes 
mean_occupancy_interpro <- mean_occupancy_interpro %>% filter(count_pangenes > 5)

essential_TF_domain_occupancy <- essential_TF_ipr %>% left_join(mean_occupancy_interpro[,-4], join_by(interpro_id == id))

all_TF_domain_occupancy <- df_all_TF_interpro %>% left_join(mean_occupancy_interpro[,-4], join_by(Var1 == id))

# or

TF_interpro_occupancy <- mean_occupancy_interpro %>% filter(id %in% df_all_TF_interpro$Var1) %>% mutate(class= "Transcription factor") 



# F-box
fb_domains <- data.frame(ipr = c( "IPR041567" , "IPR044997" , "IPR036047", "IPR001810" , "IPR017451", "IPR013187" ,  "IPR006527", "IPR044595"), 
                         pfam = c("PF18511", "PF07734", "PF08268", "PF00646", "PF12937", "PF19270", rep("-",2)))


fbox_mean_interpro_occupancy <- mean_occupancy_interpro %>% filter(id %in% fb_domains$ipr) %>% mutate(class= "F-box domain")



LRR_interpro <- c("IPR001611", "IPR003591", "IPR006553" , "IPR013101", "IPR011713" , "IPR025875", "IPR026906", "IPR000372" , "IPR000483" , "IPR004830" , 
                  "IPR013210" , "IPR031782" , "IPR041283" , "IPR041302" , "IPR041403" , "IPR032675" , "IPR025875")

LRR_interpro_occupancy <- mean_occupancy_interpro %>% filter(id %in% LRR_interpro) %>% mutate(class= "Leucine-rich repeats (LRR)") 


DUF_interpro <- mean_occupancy_interpro %>% filter(grepl("DUF" , V13) | grepl("Domain of unknown function", V13) ) %>% mutate(class= "Domain of Unknown Functions (DUF)")


Kinase_mean_interpro_occupancy <- mean_occupancy_interpro %>% filter(grepl("kinase", V13 , ignore.case = TRUE)) %>% mutate(class= "Protein kinase")

disease_res <- c( "IPR002182" , "IPR044974" ,"IPR039203")     ##  2182 is NB-ARC , IPR044974 and IPR039203 disease resistance

resistance_genes_interpro <- mean_occupancy_interpro %>% filter(id %in% disease_res) %>% mutate(class= "Disease resistance domains") 

selected_interpro_occupancy <- as.data.frame(bind_rows(Kinase_mean_interpro_occupancy, DUF_interpro , 
                                                       LRR_interpro_occupancy, Kinase_mean_interpro_occupancy, 
                                                       fbox_mean_interpro_occupancy,  TF_interpro_occupancy , resistance_genes_interpro))


# numbers

selected_interpro_occupancy %>% group_by(class) %>% summarise(median_val  =median(mean_occupancy),
                                                              IQR = IQR(mean_occupancy))


#Fig 4D
plot_ipr_violin <- ggplot(selected_interpro_occupancy, aes(x = reorder(class, mean_occupancy), y = mean_occupancy)) +
  geom_violin(aes(fill = class), width=1, alpha = 0.3 ) + ylim(0,16) + 
  theme_cowplot(font_size = 14) +
  labs( y="Mean domain occupancy", x="Terms",fill = "Interpro description" ) + 
  scale_colour_manual(values = cbPalette) +
  stat_summary(fun="median", geom = "point",  size = 3, colour = "black")  + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(legend.position = "none")


# Uncomment to export image
# svg("./plots/family_pfam_violin_new_interpro.svg")
# plot_ipr_violin
# dev.off()



# Fig 4e
# DUF analysis

DUF_interpro <- mean_occupancy_interpro %>% filter(grepl("DUF" , V13) | grepl("Domain of unknown function", V13)) %>% mutate(class= "DUF")

DUF_label <- DUF_interpro %>% 
  filter(count_pangenes > 40) %>%
  mutate(plotname = as.character(V13))

DUF_interpro <- DUF_interpro %>% mutate(plotname = ifelse(V13 %in% DUF_label$V13, V13 , ""))

DUF_interpro$plotname <- str_replace_all(DUF_interpro$plotname, 
                                         c("Protein of unknown function" = "", 
                                           "Domain of unknown function" = "",
                                           "Domain unknown function" = "", 
                                           ", Oryza sativa" = "",
                                           "^ " = ""  , ", plant" = "") )


color_range = range(DUF_interpro$count_pangenes) / max(DUF_interpro$count_pangenes) *  6

duf_interpro_scatter_plot <- ggplot(DUF_interpro, aes(x = count_pangenes, y = mean_occupancy , color = mean_occupancy)) +
  geom_point(size = 2.5 ) + theme_cowplot() + 
  # scale_fill_continuous(range = color_range) +
  labs(x = "Number of pan-genes", y = "Mean domain occupancy" ) +
  geom_text_repel(aes(x = count_pangenes, label = plotname), colour = "black", size = 4, hjust = 0) + 
  scale_colour_gradientn(colours = c("blue", "orange")) +
  theme(legend.position = c(0.8, 0.9), legend.direction	 = "vertical")

# uncomment to export image
# svg("./plots/duf_interpro_scatter_plot.svg")
# duf_interpro_scatter_plot
# dev.off()


