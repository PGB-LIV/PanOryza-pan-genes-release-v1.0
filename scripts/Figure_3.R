## Figure 3

#Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(gridExtra))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(reshape2))


load(file = "core_workspace.RData") #  load the core set of files
clusters_dir <- "./Os4530.POR.1" # directory with version 1 of pan-genes
source("scripts/figure2_accessory_function.R")

###colors
color_scheme <- read.delim("core_files/color_scheme.txt")

col_palette <- color_scheme$color

names(col_palette) <- color_scheme$genome

col_palette2 <- color_scheme$color2

names(col_palette2) <- color_scheme$subspecies  


# get the data
peptide_mapping <- read.delim("core_files/peptide_mapping.tsv", header=FALSE)

# remove DECOY
decoy_filtered <- peptide_mapping[!(grepl("DECOY", peptide_mapping$V3)),]

## change old rice magic 16 ids to new
gramene_old_new <- read.csv("core_files/gramene_old_new.csv")

decoy_filtered$new_id <- decoy_filtered$V3

for (i in 1:nrow(gramene_old_new)) {
  search_keyword <- gramene_old_new[i,1]
  replace_keyword <- gramene_old_new[i,2]
  decoy_filtered$new_id <- str_replace_all(decoy_filtered$new_id, search_keyword, replace_keyword)
}

decoy_filtered_1 <- decoy_filtered %>% left_join(clusters_summary, join_by("new_id" == "transcript")) %>% select(-V3, -File)

decoy_filtered_1 <- decoy_filtered_1 %>% filter(!(grepl("_Ung", new_id) | grepl("_UnG", new_id))) ### remove rows with unmapped transcripts # _UnG containing ids

decoy_filtered_1 <- decoy_filtered_1[!(is.na(decoy_filtered_1$Pan.id)),] ## removed NA pan-gene ids



### or if continuing from Figure 2


head(decoy_filtered_1)

decoy_filtered_1$accession <- decoy_filtered_1$genome

decoy_filtered_1 <- decoy_filtered_1 %>% 
  mutate(., accession = recode(accession, "Nipponbare(RAPDB)" = "Nipponbare(merged)", "Nipponbare(MSU)" = "Nipponbare(merged)" , "Nipponbare(Gramene)" = "Nipponbare(merged)"))

## shortened df
decoy_filtered_2 <- decoy_filtered_1 %>% select(V1, new_id, Pan.id, taxa_size, occupancy, merged.or.gene_id , genome, gene_id, accession) %>%
  rename(peptide = "V1", transcript = "new_id")


### taking only the first transcript
new_df <- decoy_filtered_1[,c(1,2,7,8, 15)] %>% distinct() 


new_df1 <- new_df %>% arrange(V2, new_id) %>% group_by(V1, genome) %>% slice(1)

new_df1 <- new_df1 %>% na.omit()

new_df2 <- new_df1 %>% group_by(new_id, genome) %>% summarise(peptide_count_per_protein  = n_distinct(V1))


peptide_count_summary_1st_entry <- new_df2 %>% filter(peptide_count_per_protein >= 2) %>% group_by(genome) %>% summarise(counts_by_genome = n_distinct(new_id))

peptide_count_summary_1st_entry <- peptide_count_summary_1st_entry %>% left_join(color_scheme[,c(1,4)], by = "genome")


#plot Fig 3a
plot_evidence <- ggplot(peptide_count_summary_1st_entry, aes(x=reorder(genome, desc(counts_by_genome)),y=counts_by_genome, fill=subspecies)) +
  geom_bar(stat="identity", width = 0.8) + scale_fill_manual(values = col_palette2) +
  labs( y="Number of proteins", x="genome") +
  theme_cowplot(font_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs( fill = "Subpopulation") 

# uncomment to export figure
# svg("figure3a.svg", width=7,height=5.5, pointsize = 10)
# plot_evidence_new
# dev.off()



# plot Fig 3b

### numbers of distinct peptides per gene

peptide_count_summary4 <- decoy_filtered_2 %>% group_by(merged.or.gene_id, accession) %>% summarise(peptide_count_per_gene  = n_distinct(peptide))

peptide_count_summary4 <- peptide_count_summary4 %>% left_join(clusters_summary[, c(1,3:6)], by = "merged.or.gene_id" )

pep_plot_occupancy <- boxplot_pangenes_with_occupancy(peptide_count_summary4, 
                                                      y_axis = "peptide_count_per_gene", 
                                                      x_axis = "taxa_size",
                                                      fill = "occupancy", y_label = "log10(Number of peptides per gene)",
                                                      x_label = "Pan-gene occupancy") 

pep_plot_occupancy <- pep_plot_occupancy +  scale_y_continuous(trans = 'log10')

pep_plot_occupancy <- add_significance_bar(pep_plot_occupancy, peptide_count_summary4, x_arg = "occupancy", 
                                           y_arg = "peptide_count_per_gene", pos_a = 300, pos_b = 400 , vjust = 50 , vjust_text = 50 )


# uncomment to export figure
# svg("figure3b.svg", width=7,height=6, pointsize = 10)
# pep_plot_occupancy
# dev.off()


#Plot Fig 3c

# number of peptides mapping to number of genomes and accessions

peptide_map_plot2_data <- decoy_filtered_2 %>% group_by(peptide) %>% summarise(counts = n_distinct(accession)) ### number of accessions(total 16) a peptide maps to

peptide_map_plot2_data <- peptide_map_plot2_data %>% group_by(counts) %>% summarise(peptide_number = length(peptide))


peptide_count_plot2 <- bar_plot_pangenes_with_occupancy(df = peptide_map_plot2_data , y_axis = "peptide_number", x_axis = "counts", 
                                                        fill = "counts", y_label = "count of peptides", x_label = "Number of accessions") + labs(fill = "accessions")

#or

peptide_count_plot2 <- ggplot(peptide_map_plot2_data) +
  geom_bar(aes(x = as.factor(counts), y = peptide_number ) , stat="summary", width = 0.8, alpha = 1, na.rm = TRUE) +
  theme_cowplot(font_size = 14) + scale_fill_manual(values = "darkgrey") +
  labs( y="count of peptides", x="Number of accessions")

# uncomment to export figure
# svg("figure3c.svg", width=7,height=6, pointsize = 10)
# peptide_count_plot2
# dev.off()


# Plot Fig 3d

#aed plot 

aed_scores <- read.delim("core_files/rice_trace_reorder_with_AED_rerankTranscript_unique_peptide.out", sep = "\t")

clusters_summary <- clusters_summary %>% mutate(accession  = genome ) %>% 
  mutate(., accession = recode(accession, "Nipponbare(RAPDB)" = "Nipponbare(merged)", "Nipponbare(MSU)" = "Nipponbare(merged)" , "Nipponbare(Gramene)" = "Nipponbare(merged)"))


Nip_cluster_pep_aed <- clusters_summary %>% select(Pan.id, taxa_size, occupancy, transcript, merged.or.gene_id, genome, gene_id, accession) %>% 
  filter(accession == "Nipponbare(merged)") 

Nip_cluster_pep_aed <- Nip_cluster_pep_aed %>%
  left_join(aed_scores[,c(2,3,10,12,13)], join_by("transcript" == "Transcript")) %>%
  left_join(decoy_filtered_2[,c(1,2)], by = "transcript")



Nip_cluster_pep_aed$pep_evidence <- ifelse(is.na(Nip_cluster_pep_aed$peptide) , paste0("no"), paste0("yes"))

# 60,400 pangenes containing  62225 Nipponbare transcripts total : 23288 pangenes with evidence, 23742 genes with evidence.

# 16701 PANGENES WITH EVIDENCE IN rapdb ; 22155 in MSU


Nip_cluster_pep_aed_summary <- Nip_cluster_pep_aed %>% group_by(taxa_size, pep_evidence) %>% summarise(mean_min_aed = mean(min_aed))

plot_aed <- ggplot(Nip_cluster_pep_aed_summary, aes(x = as.factor(taxa_size), y = mean_min_aed , color = pep_evidence)) + geom_point(size = 4) + 
  labs(color="peptide evidence", y="minimum AED score", x="Pan-gene occupancy") +
  scale_color_manual(values = c("red", "blue")) +theme_cowplot() + ylim(0,1)
  
  
# ADD LM 
  
proteins_wo_pep_support <- Nip_cluster_pep_aed[Nip_cluster_pep_aed$pep_evidence == "no",] 

proteins_with_pep_support <- Nip_cluster_pep_aed[Nip_cluster_pep_aed$pep_evidence == "yes",] 


model_pep <- lm(min_aed ~ taxa_size, data = proteins_with_pep_support)
model_wo_pep <- lm(min_aed ~ taxa_size, data = proteins_wo_pep_support)

summary(lm(min_aed ~ 0+pep_evidence, data = Nip_cluster_pep_aed))

summary(model_pep)
summary(model_wo_pep)

xmin = min(model_pep$model$taxa_size)
xmax = max(model_pep$model$taxa_size)

predicted_pep <- data.frame(taxa_size = seq(xmin, xmax, length.out = 100))
predicted_pep$count <- predict(model_pep, predicted_pep)
predicted_pep$pep_evidence <- "yes"

predicted_pep_wo <- data.frame(taxa_size = seq(min(model_wo_pep$model$taxa_size), max(model_wo_pep$model$taxa_size), length.out = 100))
predicted_pep_wo$count <- predict(model_wo_pep, predicted_pep_wo)
predicted_pep_wo$pep_evidence <- "no"

predicted <- bind_rows(predicted_pep,predicted_pep_wo )

plot_aed <- plot_aed + 
  geom_line(aes(as.numeric(taxa_size), y = count , colour = pep_evidence), data = predicted, size = 0.75, linetype = 2) 


# uncomment to export figure
# svg("/mnt/hc-storage/users/esharma/rice/ms/plots/figure3/figure3_aed_plot.svg", width=7,height=6, pointsize = 10)
# plot_aed
# dev.off()


plots_grid <- ggarrange(
  plot_evidence_new,pep_plot_occupancy, peptide_count_plot2, plot_aed,
  labels = c("A", "B", "C", "D") ,
  ncol=2,
  nrow=2,
  common.legend = FALSE, legend = "right"
)

# uncomment to export figure
# svg("figure 3.svg", width=9.5,height=7, pointsize = 8)
# plots_grid
# dev.off()

