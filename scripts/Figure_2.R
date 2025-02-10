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

### take out RAPDB transcripts contained in the pan-gene clusters and calculate various counts

#first add original gene ids , merged ids already present () already added in core workspace)

# clusters_summary <- clusters_summary %>% left_join(magic18_models , join_by("transcript" == "transcipt")) %>% select(-species) %>% distinct()

head(clusters_summary) #check

rapdb_clusters_counts <- clusters_summary %>% filter(genome == "Nipponbare(RAPDB)") %>% 
  select(Pan.id, taxa_size , occupancy, gene_id, transcript) %>% 
  group_by(taxa_size) %>% summarise(pan_counts = n_distinct(Pan.id),
                                    transcript_counts = n_distinct(transcript),
                                    gene_counts = n_distinct(gene_id))

### here 42474 transcripts of RAPDB grouped into 32688 clusters
sum(rapdb_clusters_counts$pan_counts) # 32688  , 5244 pan-genes are cloud and 21767 are core. 
sum(rapdb_clusters_counts$gene_counts) # 35694 genes 
sum(rapdb_clusters_counts$transcript_counts) # 42472 transcripts ; 5391 in cloud & 29977 in core


head(magic18_interpro_clusters) ## interpro annotated clusters file

# check
magic18_interpro_clusters %>% filter(genome == "Nipponbare(RAPDB)") %>% 
  select(Pan.id, taxa_size , gene_id, original_transcript) %>% 
  group_by(taxa_size) %>% summarise(pan_counts = n_distinct(Pan.id),
                                    transcript_counts = n_distinct(original_transcript),
                                    gene_counts = n_distinct(gene_id))

# take out RAPDB ones
rapdb_clusters <- magic18_interpro_clusters %>% filter(genome == "Nipponbare(RAPDB)") 


rapdb_clusters$occupancy = with(rapdb_clusters, ifelse(taxa_size == "16", "core",
                                                       ifelse(taxa_size == "15", "softcore",
                                                              ifelse(taxa_size == "1", "cloud", 
                                                                     ifelse(taxa_size == "2", "cloud" , "shell")))))

rapdb_clusters_interpro <- rapdb_clusters %>%
  select(Pan.id, V12, V13, taxa_size, occupancy, gene_id, original_transcript) %>% distinct()

total_pan_counts <- rapdb_clusters_interpro %>%
  group_by(taxa_size, occupancy) %>%
  summarise(pan_counts = n_distinct(Pan.id),
            transcript_count = n_distinct(original_transcript))

IPR_filter_pan_counts <- rapdb_clusters_interpro %>%
  group_by(taxa_size, occupancy) %>% filter(grepl("^IPR", V12)) %>%
  summarise(pan_counts_f = n_distinct(Pan.id),
            transcript_count_f = n_distinct(original_transcript))

rap_ipr_sum <- merge(total_pan_counts, IPR_filter_pan_counts, by = "taxa_size")

rap_ipr_sum <- rap_ipr_sum %>% mutate(proportion_pan_with_ipr = pan_counts_f/pan_counts ,
                                      proportion_transcripts_with_ipr = transcript_count_f/transcript_count)


# plot of pan-genes with an interpro id
plot_interpro <- ggplot(rap_ipr_sum, aes(x=as.factor(taxa_size),y=(100*proportion_pan_with_ipr), fill=occupancy.x)) +
  geom_bar(stat="identity", width = 0.8)+
  labs(fill="occupancy", y="Percentage of Pan-genes\nwith an Interpro id", x="Pan-gene occupancy") +
  theme_cowplot(font_size = 14)

#

## get aed scores  

aed_scores <- read.delim("core_files/rice_trace_reorder_with_AED_rerankTranscript_unique_peptide.out", sep = "\t")

rapdb_clusters_aed <- merge(rapdb_clusters, aed_scores, by.x = "original_transcript", by.y = "Transcript")

rapdb_clusters_aed_red <- rapdb_clusters_aed %>% select(Pan.id, genome, gene_id,Protein_length, original_transcript, taxa_size, occupancy, Avg_AED_Gene, min_aed)


### table of results based on taxa size
aed_result_df1 <- rapdb_clusters_aed_red %>%
  group_by(taxa_size) %>%
  summarise(gene_count = length(unique(gene_id)),
            transcript_count = length(unique(original_transcript)),
            sd = sd(min_aed, na.rm = TRUE),
            min_aed = mean(min_aed)) %>%
  mutate(occupancy = ifelse(taxa_size == "16", "core",
                            ifelse(taxa_size == "15", "softcore",
                                   ifelse(taxa_size == "1", "cloud", 
                                          ifelse(taxa_size == "2", "cloud" , "shell")))))
### table of results based on occupancy class
aed_result_df2 <- rapdb_clusters_aed_red %>%
  group_by(occupancy) %>%
  summarise(gene_count = length(unique(gene_id)),
            transcript_count = length(unique(original_transcript)),
            sd = sd(min_aed, na.rm = TRUE),
            min_aed = mean(min_aed))


##plot aed scores vs pan-genes

plot_aed <- bar_plot_pangenes_with_occupancy(rapdb_clusters_aed_red, y_axis = "min_aed", x_axis = "taxa_size", fill = "occupancy", y_label = "Minimum AED score", x_label = "Pan-gene occupancy") 

plot_aed <- plot_aed +  stat_summary(aes(x=as.factor(taxa_size),y=min_aed), data = rapdb_clusters_aed_red, fun.data = mean_cl_boot, colour="black", geom="errorbar", width=0.1) 

plot_aed <- add_significance_bar(plot = plot_aed, df = rapdb_clusters_aed_red, x_arg = "occupancy", y_arg = "min_aed", pos_a= 0.95, pos_b = 1.0 , vjust = NA , vjust_text = NA )

plot_aed <- plot_aed + ylim(0,1.05)



### plot 2c protein length

#get data
aaplot_data <- rapdb_clusters_aed_red[rapdb_clusters_aed_red$Protein_length <= 2000, ]

plot3_aa_boxplot <- boxplot_pangenes_with_occupancy(aaplot_data, y_axis = "Protein_length", x_axis = "taxa_size", fill = "occupancy", y_label = "Number of Amino acids", x_label = "Pan-gene occupancy" )

plot3_aa_boxplot <- add_significance_bar(plot = plot3_aa_boxplot, df = rapdb_clusters_aed_red, x_arg = "occupancy", y_arg = "Protein_length", pos_a= 2100, pos_b = 2300 , vjust = NA , vjust_text = NA )


# plot 2d and 2e : orthologs and paralog counts

rice_homologies <- read.delim("core_files/oryza_sativa_from_homologies.tsv", sep = "\t")

taxa <- fread("core_files/Species.csv", select = c("Name","Classification"))
taxa <- as.data.frame(taxa)
taxa$Name <- tolower(taxa$Name)
taxa$Name <- str_replace_all(taxa$Name, ' ', '_')

rice_homologies$homo_species_name <- rice_homologies$homology_species

split_compara <-  split(rice_homologies, rice_homologies$homology_type)



# repeating the analyses using tidy.

split_compara_summary <- lapply(split_compara, function(x) {
  
  x <- x %>% left_join(distinct(rapdb_clusters[,c(1,2,11,13,14)]), join_by("protein_stable_id" == "original_transcript"))
  x_summary <- x %>% group_by(gene_stable_id, taxa_size, occupancy) %>% summarise(counts = n_distinct(homology_gene_stable_id),
                                                                                  counts_prot = n_distinct(homology_protein_stable_id),
                                                                                  counts_homology_id = n_distinct(homology_id))
  return(x_summary)
  
})


names(split_compara_summary ) <- paste0(names(split_compara), "_summary")

one2one <- split_compara_summary$ortholog_one2one_summary %>% na.omit()

species_paralog <- split_compara_summary$within_species_paralog_summary %>% na.omit()


### tabular results

one2one %>% group_by(taxa_size) %>% summarise(median_val = median(counts),
                                              IQR = IQR(counts))

species_paralog %>% group_by(taxa_size) %>% summarise(median_val = median(counts),
                                                      IQR = IQR(counts))


### plots 2d and 2e

plot4_ortho <- boxplot_pangenes_with_occupancy(one2one, y_axis = "counts", x_axis = "taxa_size", fill = "occupancy", y_label = "Orthologs count", x_label = "Pan-gene occupancy" )

plot4_ortho <- add_significance_bar(plot = plot4_ortho, df = one2one, x_arg = "occupancy", y_arg = "counts", pos_a= 105, pos_b = 120 , vjust = NA , vjust_text = NA )


## paralog plots

paralog_plot_data <- species_paralog[species_paralog$counts <= 30, ]  

plot5_para <- boxplot_pangenes_with_occupancy(paralog_plot_data, y_axis = "counts", x_axis = "taxa_size", fill = "occupancy", y_label = "Paralogs count", x_label = "Pan-gene occupancy" )

plot5_para <- add_significance_bar(plot = plot5_para, df = species_paralog, x_arg = "occupancy", y_arg = "counts", pos_a= 32, pos_b = 36 , vjust = NA , vjust_text = NA )




#### Plot 6 AF

rice_AF <- read.csv("core_files/rice_Alphafold_mean_plDDT.csv", row.names=1) 

uniprotkb_proteome_ids <- read.delim("core_files/uniprotkb_proteome_ids.txt", na.strings = "")

uniprotkb_rice_reshaped <- melt(uniprotkb_proteome_ids, id.vars = c("Entry", "Entry.Name" ))

uniprotkb_rice_reshaped <- uniprotkb_rice_reshaped %>% dplyr::rename(gene_id = "value") %>% na.omit()

uniprotkb_rice_reshaped <- uniprotkb_rice_reshaped %>% left_join(distinct(magic18_interpro_clusters[,c(1,8,10:12)]), join_by(gene_id == gene_id)) %>% na.omit()

uniprotkb_rice_reshaped$occupancy = with(uniprotkb_rice_reshaped, ifelse(taxa_size == "16", "core",
                                                                         ifelse(taxa_size == "15", "softcore",
                                                                                ifelse(taxa_size == "1", "cloud",
                                                                                       ifelse(taxa_size == "2", "cloud" , "shell")))))

uniprot_clusters_AF <- merge(uniprotkb_rice_reshaped, rice_AF, by.x = "Entry", by.y =  "AF_id")

# tabular result AF
uniprot_clusters_AF %>% group_by(occupancy) %>% summarise(median_val = median(mean_prediction_score),
                                                          IQR = IQR(mean_prediction_score))


##plot 2g

plot6_af <- boxplot_pangenes_with_occupancy(uniprot_clusters_AF, y_axis = "mean_prediction_score", x_axis = "taxa_size", fill = "occupancy", y_label = "Mean prediction score", x_label = "Pan-gene occupancy" )

plot6_af <- add_significance_bar(plot = plot6_af, df = uniprot_clusters_AF, x_arg = "occupancy", y_arg = "mean_prediction_score", pos_a= 102, pos_b = 109 ,  vjust = NA , vjust_text = NA)




##### expression data plot : Fig. 2f

rice.FPKM <- fread("core_files/rice.FPKM.csv")

RAP.MSU <- read.delim("core_files/RAP.MSU.tab")

RAP.MSU <- RAP.MSU %>% select(-msu_transcript) %>% distinct()

print(rice.FPKM[1:4,1:10])

rice.FPKM$mean_FPKM <- rowMeans(rice.FPKM[, 2:8284])

rice.FPKM.mean <- rice.FPKM %>% select(Sample, mean_FPKM)

### all FPKM rae MSU ids, need to get rapdb ids as well
rice.FPKM.mean <- rice.FPKM.mean %>% left_join(RAP.MSU, join_by("Sample" == "msu_gene"))

msu_rap_panclusters <- magic18_interpro_clusters %>% filter(genome %in% c("Nipponbare(MSU)", "Nipponbare(RAPDB)")) %>% 
  select(Pan.id, merged.or.gene_id, genome, taxa_size, gene_id) %>% distinct()

msu_values <- msu_rap_panclusters %>% right_join(rice.FPKM.mean, join_by("gene_id" == "Sample")) %>% dplyr::rename(other_id = "RAP")

rap_values <- msu_rap_panclusters %>% right_join(rice.FPKM.mean, join_by("gene_id" == "RAP")) %>% dplyr::rename(other_id = "Sample")

expn_clusters <- bind_rows(msu_values,rap_values) %>% select(Pan.id, mean_FPKM, taxa_size) %>% distinct() %>% na.omit()

expn_clusters2 <- expn_clusters %>% mutate(log2_fpkm = log2(mean_FPKM + 1))

expn_clusters2$occupancy = with(expn_clusters2, ifelse(taxa_size == "16", "core",
                                                       ifelse(taxa_size == "15", "softcore",
                                                              ifelse(taxa_size == "1", "cloud",
                                                                     ifelse(taxa_size == "2", "cloud" , "shell")))))

expn_clusters3 <- expn_clusters2 %>% group_by(Pan.id, taxa_size, occupancy) %>% summarise(mean_log2_fpkm = mean(log2_fpkm))

expn_clusters3_fil <- expn_clusters3 %>% filter(mean_log2_fpkm < 9)


#plot

plot7_expn <- boxplot_pangenes_with_occupancy(expn_clusters3_fil, y_axis = "mean_log2_fpkm", x_axis = "taxa_size", 
                                              fill = "occupancy", y_label = "Gene expression\n(in log2(FPKM+1)", x_label = "Pan-gene occupancy" )

plot7_expn <- add_significance_bar(plot = plot7_expn, df = expn_clusters3, x_arg = "occupancy", y_arg = "mean_log2_fpkm", pos_a= 9, pos_b = 10 ,  vjust = NA , vjust_text = NA)


# next plot peptide data

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

## now 28658 pan ids and 435196 transcripts , 329739 individual gene ids(for 18 genomes) and 293844 merged gene ids(for 16 genomes)

# calculate number of unique counts of peptides per gene id for 16 genomes i.e. merged ids need to be used. 

peptide_count_summary <- decoy_filtered_1 %>% group_by(merged.or.gene_id, Pan.id) %>% summarise(peptide_count_per_gene  = n_distinct(V1))

# add taxa and occupancy cols
peptide_count_summary <- clusters_summary %>% select(Pan.id, taxa_size, occupancy) %>% 
  distinct() %>% 
  right_join(peptide_count_summary, by = "Pan.id")

### taking total pan.ids as number of total clusters to calculate proportion against those with peptide evidence 
# i.e evidence of more than 2 peptides per gene of any cluster 

total_cluster_counts <- clusters_summary %>% group_by(taxa_size, occupancy) %>% summarise(total_cluster_counts = n_distinct(Pan.id))

cluster_counts_minimum_2_pep_per_gene <- peptide_count_summary %>% filter(peptide_count_per_gene >= 2) %>% group_by(taxa_size, occupancy) %>% summarise(counts_per_gene_min_2_pep = n_distinct(Pan.id))


proportion_clusters_2pep_evidence <- merge(total_cluster_counts, cluster_counts_minimum_2_pep_per_gene[,c(1,3)], by = "taxa_size")

proportion_clusters_2pep_evidence <- proportion_clusters_2pep_evidence %>% mutate(proportion_percent = 100*(counts_per_gene_min_2_pep/total_cluster_counts))




plot8_peptide <- bar_plot_pangenes_with_occupancy(proportion_clusters_2pep_evidence, "proportion_percent", "taxa_size", 
                                                  fill = "occupancy", y_label = "Percentage of clusters\nwith peptide count > 2", x_label = "Pan-gene occupancy")



plots_grid <- ggarrange(
  plot_aed, plot_interpro,
  plot3_aa_boxplot, 
  plot4_ortho, 
  plot5_para, 
  plot7_expn, plot6_af, 
  plot8_peptide,
  labels = c("A", "B", "C", "D", "E", "F", "G", "H") ,
  ncol=2,
  nrow=4,
  common.legend = TRUE, legend = "bottom" ) +
  theme(legend.position = "bottom", 
        legend.justification = "center")

# (Optional) Uncomment to export the figure
# svg("figure2.svg", width=9,height=12, pointsize = 9)
# plots_grid
# dev.off()
