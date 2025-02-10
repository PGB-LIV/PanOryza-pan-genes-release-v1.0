### functions file
 
### functions to plot

boxplot_pangenes_with_occupancy <- function(df, y_axis, x_axis, fill, y_label, x_label){ ## defaults to fill by occupancy and calculate stats based on occupancy
  plot <- ggplot(df) +
    geom_boxplot(aes(x = as.factor(.data[[x_axis]]), y = .data[[y_axis]] , fill = .data[[fill]] ), width = 0.4, alpha = 1, na.rm = TRUE) +
    geom_violin(aes(x = as.factor(.data[[x_axis]]), y = .data[[y_axis]] , fill = .data[[fill]] ), width=1, alpha = 0.3, color = NA) +  
    theme_cowplot(font_size = 14) +
    labs(fill="occupancy", y=y_label, x=x_label)
  return(plot)
}

bar_plot_pangenes_with_occupancy <- function(df, y_axis, x_axis, fill, y_label, x_label){ ## defaults to fill by occupancy and calculate stats based on occupancy
  plot <- ggplot(df) +
    geom_bar(aes(x = as.factor(.data[[x_axis]]), y = .data[[y_axis]] , fill = .data[[fill]] ), stat="summary", width = 0.8, alpha = 1, na.rm = TRUE) +
    theme_cowplot(font_size = 14) +
    labs(fill="occupancy", y=y_label, x=x_label)
  return(plot)
}


## function to add segment lines showing significance

add_significance_bar <- function(plot, df, x_arg, y_arg, pos_a, pos_b , vjust , vjust_text){ # a and b are y-axis values for the two segments , p-value based on tukey test, x_arg and y_arg are x and y for stat test
  
  model <- aov(df[[y_arg]] ~ df[[x_arg]], df)
  tukey_res <- TukeyHSD(model)
  
  pv_core <- ifelse(tukey_res$`df[[x_arg]]`["core-cloud",][["p adj"]] == 0 , "0", format(round(tukey_res$`df[[x_arg]]`["core-cloud",][["p adj"]],3)))
  #pv_soft <- ifelse(tukey_res$`df[[x_arg]]`["softcore-cloud",][["p adj"]] == 0, "0", format(round(tukey_res$`df[[x_arg]]`["softcore-cloud",][["p adj"]],3)))
  
  vjust = ifelse(is.na(vjust) , (2.0 * (pos_a + pos_b /2))/100 , vjust)
  vjust_text = ifelse(is.na(vjust_text), (3.0 * (pos_a + pos_b /2))/100 , vjust_text )
  
  plot + 
    #annotate("segment", x = 1.5, xend = 15, y = pos_a, colour = "black", size=0.5, alpha=0.7) +
    #annotate("text", x = 8, y = pos_a+vjust_text, label = paste("p :", pv_soft ) , colour = "black")  +
    #annotate("segment", x = 1.5, xend = 1.5, y = pos_a-vjust, yend = pos_a, colour = "black", size=0.5, alpha=0.7) +
    #annotate("segment", x = 1, xend = 2, y = pos_a-vjust, yend = pos_a-vjust, colour = "black", size=0.5, alpha=0.7) +
    #annotate("segment", x = 15, xend = 15, y = pos_a-vjust, yend = pos_a, colour = "black", size=0.5, alpha=0.7) +
    annotate("segment", x = 1.5, xend = 16, y = pos_b, colour = "black", size=0.5, alpha=0.7) +
    annotate("text", x = 8, y = pos_b+vjust_text, label = paste("p :", pv_core ) , colour = "black")  +
    annotate("segment", x = 1.5, xend = 1.5, y = pos_b-vjust, yend = pos_b, colour = "black", size=0.5, alpha=0.7) +
    annotate("segment", x = 1, xend = 2, y = pos_b-vjust, yend = pos_b-vjust, colour = "black", size=0.5, alpha=0.7) +
    annotate("segment", x = 16, xend = 16, y = pos_b-vjust, yend = pos_b, colour = "black", size=0.5, alpha=0.7) 
  
}



## function to calculate proportion of interpro annotated pan-genes 

get_percent_interpro <- function(x , genome){ # x is interpro annotated file with pan.ids and genome column 
  ifelse(genome %in% unique(magic18_interpro_updated_clusters1$genome), 
         x_clusters <- x[x$genome == genome,] %>%
    select(Pan.id, V12, V13, taxa_size,  original_transcript) %>% distinct() , 
    print("genome not present"))
  
  total_pan_counts <- x_clusters %>%
    group_by(taxa_size) %>%
    summarise(pan_counts = n_distinct(Pan.id),
              transcript_count = n_distinct(original_transcript))
  
  IPR_filter_pan_counts <- x_clusters %>%
    group_by(taxa_size) %>% filter(grepl("^IPR", V12)) %>%
    summarise(pan_counts_f = n_distinct(Pan.id),
              transcript_count_f = n_distinct(original_transcript))
  
  total_ipr_sum <- merge(total_pan_counts, IPR_filter_pan_counts, by = "taxa_size")
  
  total_ipr_sum <- total_ipr_sum %>% mutate(proportion_pan_with_ipr = pan_counts_f/pan_counts ,
                                        proportion_transcripts_with_ipr = transcript_count_f/transcript_count)
  
  total_ipr_sum$occupancy = with(total_ipr_sum, ifelse(taxa_size == "16", "core",
                                                        ifelse(taxa_size == "15", "softcore",
                                                               ifelse(taxa_size == "1", "cloud", 
                                                                      ifelse(taxa_size == "2", "cloud" , "shell")))))
  
  return(total_ipr_sum)
}

# 
# msu_ipr_sum <- get_percent_interpro(magic18_interpro_updated_clusters1, "Nipponbare(MSU)")
# 
# msu_ipr_sum$genome  = "Nipponbare(MSU)"
# 
# ggplot(msu_ipr_sum, aes(x=as.factor(taxa_size),y=(100*proportion_pan_with_ipr), fill=occupancy)) +
#   geom_bar(stat="identity", width = 0.8)+
#   labs(fill="occupancy", y="Interpro ID %", x="Pan-gene occupancy", title = "MSU") +
#   theme_cowplot(font_size = 14)
# 
# 
# gramene_ipr_sum <- get_percent_interpro(magic18_interpro_updated_clusters1, "Nipponbare(Gramene)")
# gramene_ipr_sum$genome = "Nipponbare(Gramene)"
# 
# ggplot(gramene_ipr_sum, aes(x=as.factor(taxa_size),y=(100*proportion_pan_with_ipr), fill=occupancy)) +
#   geom_bar(stat="identity", width = 0.8)+
#   labs(fill="occupancy", y="Interpro ID %", x="Pan-gene occupancy") +
#   theme_cowplot(font_size = 14)
# 
# N22_ipr_sum <- get_percent_interpro(magic18_interpro_updated_clusters1, "N22")
# 
# ggplot(N22_ipr_sum, aes(x=as.factor(taxa_size),y=(100*proportion_pan_with_ipr), fill=occupancy)) +
#   geom_bar(stat="identity", width = 0.8)+
#   labs(fill="occupancy", y="Interpro ID %", x="Pan-gene occupancy") +
#   theme_cowplot(font_size = 14)
# 
# rap_ipr_sum <- get_percent_interpro(magic18_interpro_updated_clusters1, "Nipponbare(RAPDB)")
# 
# rap_ipr_sum$genome = "Nipponbare(RAPDB)"


## do a plot shwoing all 3, what is the proportion of cloud genes with and without interpor/ pfam domains. 

# nip_merged <- bind_rows(rap_ipr_sum, gramene_ipr_sum , msu_ipr_sum ) %>% mutate(pan_gene_without_interpro = pan_counts - pan_counts_f)
# 
# pivot_longer(nip_merged, cols = c(pan_counts, pan_gene_without_interpro) ) %>%  arrange(name, desc(value)) %>%
# ggplot(., aes(x=as.factor(taxa_size),y=value, fill=name)) +
#   geom_col(position = "stack", width = 0.6)+
#   labs(fill="occupancy", y="pan_gene_counts", x="Pan-gene occupancy") +
#   theme_cowplot(font_size = 14) + facet_wrap(~ genome)

# write.table(rapdb_clusters_interpro, "../data/fig2_rapdb_clusters_interpro.tsv", sep = "\t",  quote = FALSE, row.names = FALSE) # only for testing

# rapdb_clusters_with_pfam
# rapdb_clusters_with_pfam <- rapdb_clusters_with_pfam %>% rename(original_transcript = Transcript_ID ,  Pan.id =  cluster_id)
# rapdb_clusters_with_pfam$genome <- "Nipponbare(RAPDB)"
# 
# get_percent_interpro(rapdb_clusters_with_pfam, "Nipponbare(RAPDB)")
# 
# 
# 
# genome_ipr_list <- split(magic18_interpro_updated_clusters1, magic18_interpro_updated_clusters1$genome)

# 
# 
# ### check the code
# #genome_ipr_list1 <- for (name in names(genome_ipr_list)){
#   lapply(genome_ipr_list, FUN = function(x , genome = name){ # x is interpro annotated file with pan.ids and genome column 
#     ifelse(genome %in% unique(magic18_interpro_updated_clusters1$genome), 
#            x_clusters <- x[x$genome == genome,] %>%
#              select(Pan.id, V12, V13, taxa_size,  original_transcript) %>% distinct() , 
#            print("genome not present"))
#     
#     total_pan_counts <- x_clusters %>%
#       group_by(taxa_size) %>%
#       summarise(pan_counts = n_distinct(Pan.id),
#                 transcript_count = n_distinct(original_transcript))
#     
#     IPR_filter_pan_counts <- x_clusters %>%
#       group_by(taxa_size) %>% filter(grepl("^IPR", V12)) %>%
#       summarise(pan_counts_f = n_distinct(Pan.id),
#                 transcript_count_f = n_distinct(original_transcript))
#     
#     total_ipr_sum <- merge(total_pan_counts, IPR_filter_pan_counts, by = "taxa_size")
#     
#     total_ipr_sum <- total_ipr_sum %>% mutate(proportion_pan_with_ipr = pan_counts_f/pan_counts ,
#                                               proportion_transcripts_with_ipr = transcript_count_f/transcript_count)
#     
#     total_ipr_sum$occupancy = with(total_ipr_sum, ifelse(taxa_size == "16", "core",
#                                                          ifelse(taxa_size == "15", "softcore",
#                                                                 ifelse(taxa_size == "1", "cloud", 
#                                                                        ifelse(taxa_size == "2", "cloud" , "shell")))))
#     
#     return(total_ipr_sum)
#   })
# }
# 
