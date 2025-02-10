summarise_interpro_pfam <- function(x){ ## x is magic18_interpro_clusters style df with Pan.id, original_transcript,V12, V13, gene_id,genome, Merged.gene.ID, taxa_size etc.
  
  if (!"accession" %in% colnames(x)) {
    x$accession <- x$genome
  x <- x %>% 
    mutate(., accession = recode(accession, "Nipponbare(RAPDB)" = "Nipponbare", "Nipponbare(MSU)" = "Nipponbare" , "Nipponbare(Gramene)" = "Nipponbare"))
  }
  
  #### Generate summary of pfam domains collpases one line per pan-gene. each pfam domain and corresponding transcripts, accessions for that pfam  separated by ";" . Those associated with same
  # pfam id are separated by comma ",". same way for interpo ids. 
  
  pfam_summary <- x %>% filter(V4 == "Pfam") %>%
    group_by(Pan.id, V5) %>%
    summarise(pfam_description = unique(V6) %>% toString(),
              proteins_with_pfam = unique(original_transcript) %>% toString(),
              accessions_with_pfam = unique(accession) %>% toString())
  
  #print(head(pfam_summary))
  
  #### calculate pfam domain mean occupancy and then add to the pfam_summary
  x_msu_rap_excl <- x %>% 
    filter(!(genome %in% c("Nipponbare(MSU)" , "Nipponbare(RAPDB)"))) 
  
  
  
  domain_occupancy_table_pfam <- x_msu_rap_excl %>%
    filter(V4 == "Pfam") %>% select(Pan.id, V5, genome) %>% distinct() %>%
    group_by(V5, Pan.id ) %>%
    summarize(domain_occupancy = length(genome))
  
  mean_occupancy_pfam <- domain_occupancy_table_pfam %>%
    group_by(V5) %>%
    summarise(mean_occupancy_pfam  =  mean(domain_occupancy)) %>% left_join(pfam_annot, by = "V5")
  
  pfam_summary <- pfam_summary %>% left_join(mean_occupancy_pfam[,1:2], by = "V5")
  
  
  pfam_summary_list <- split(pfam_summary, pfam_summary$Pan.id)  
  
  pb <- progress_bar$new(
    format = "  working [:bar] :percent eta: :eta",
    total = length(pfam_summary_list), clear = FALSE, width= 60)
  
  pfam_summary1 <- lapply(pfam_summary_list, function(dataframe) {
    pb$tick()
    dataframe <- apply(dataframe, 2, function(x) paste(unique(x), collapse = ";"))
  }) %>% bind_rows()
  
  pfam_summary_counts <- x %>% filter(V4 == "Pfam") %>% group_by(Pan.id,) %>% summarise(count_pfam_domains = n_distinct(V5))
  
  pfam_summary2 <- pfam_summary1 %>% full_join(pfam_summary_counts, by = "Pan.id")
  
  ## same for interpro
  
  IPR_summary <- x %>% 
    group_by(Pan.id, V12) %>% filter(!(V12 == "-")) %>%
    summarise(interpro_description = unique(V13) %>% toString(),
              proteins_with_interpro = unique(original_transcript) %>% toString(),
              accessions_with_interpro = unique(accession) %>% toString())
  
  #print(head(IPR_summary))
  
  ####### calculate interpro domain mean occupancy and add to the table above
  
  domain_occupancy_table_interpro <- x_msu_rap_excl %>%
    filter(!(V12 == "-")) %>%
    select(Pan.id, V12, genome) %>% distinct() %>%
    group_by(V12, Pan.id ) %>%
    summarize(domain_occupancy = length(genome))
  
  
  mean_occupancy_interpro <- domain_occupancy_table_interpro %>%
    group_by(V12) %>%
    summarise(mean_occupancy_interpro  =  mean(domain_occupancy)) %>% left_join(ipr_annot, by = "V12")
  
  
  IPR_summary <- IPR_summary %>% left_join(mean_occupancy_interpro[,1:2], by = "V12")
  
  IPR_summary_list <- IPR_summary %>% split(., .$Pan.id) 
  
  pb <- progress_bar$new(
    format = "  working [:bar] :percent eta: :eta",
    total = length(IPR_summary_list), clear = FALSE, width= 60)
  
  IPR_summary1 <- lapply(IPR_summary_list, function(dataframe) {
    pb$tick()
    dataframe <- apply(dataframe, 2, function(x) paste(unique(x), collapse = ";"))
  }) %>% bind_rows()
  
  IPR_summary_counts <- x %>% filter(grepl("^IPR", V12)) %>% group_by(Pan.id) %>% summarise(count_interpro_domains = n_distinct(V12))
  
  IPR_summary2 <- IPR_summary1 %>% full_join(IPR_summary_counts, by = "Pan.id")
  
  # merge both pfam and interpro
  merged_summary_interpro <- list(IPR_summary2, pfam_summary2, distinct(x %>% select(Pan.id, taxa_size ) )) %>% reduce(full_join, by='Pan.id')
  merged_summary_interpro <- merged_summary_interpro %>% rename(pfam_ids = "V5", interpro_ids = "V12" , pangene_occupancy = "taxa_size")
  
  merged_summary_interpro <- merged_summary_interpro[!(is.na(merged_summary_interpro$Pan.id)),] # remove those pan-genes that are not present in any pan-gene cluster
  
  
  #### analyse interpro using most common domain
  
  df <- x %>% select(Pan.id, V12, gene_id, genome, taxa_size) %>% distinct()
  
  # identify the most common domain in all genes of a clusters
  res <- as.data.frame(df %>%
                         group_by(Pan.id) %>% filter(!(V12 == "-")) %>% 
                         summarise(most_common_value = names(which.max(table(V12)))) %>% na.omit())
  
  print(head(res))
  
  split_list <- split(x, x$Pan.id)
  split_list <- split_list[names(split_list) %in% res$Pan.id]
  
  pb <- progress_bar$new(
    format = "  working [:bar] :percent eta: :eta",
    total = length(res$Pan.id), clear = FALSE, width= 60)
  
  
  filtered_list <- lapply(names(split_list), function(x) {
    pb$tick()
    filter_value <- res %>% filter(Pan.id == x) %>% pull(most_common_value)
    split_list[[x]] %>% filter(V12 == filter_value) 
  })
  
  pb <- progress_bar$new(
    format = "  working [:bar] :percent eta: :eta",
    total = length(res$Pan.id), clear = FALSE, width= 60)
  
  filtered_list1 <- lapply(filtered_list, function(x) {
    pb$tick()
    x %>%
      group_by(gene_id) %>%
      summarise(count = n_distinct(gene_id))
  }) 
  
  pb <- progress_bar$new(
    format = "  working [:bar] :percent eta: :eta",
    total = length(res$Pan.id), clear = FALSE, width= 60)
  
  options(dplyr.summarise.inform = FALSE)
  
  filtered_list2 <- lapply(filtered_list, function(x) {
    pb$tick()
    x %>%
      group_by(Pan.id, V12) %>%
      summarise(accession_counts_with_most_common_interpro = n_distinct(accession) ,
                genomes_with_most_common_interpro = n_distinct(genome) ,
                accesion_names_with_common_interpro_domain = unique(accession) %>% toString())
  })
  
  
  # Assign names back to the list
  names(filtered_list1) <- names(split_list)
  names(filtered_list2) <- names(split_list)
  
  joined_df1 <- bind_rows(filtered_list1)  #### 1st calc
  summary_most_common <- bind_rows(filtered_list2) ## 2nd calc
  
  if ("species" %in% colnames(magic18_models)){
    magic18_models <- magic18_models %>% rename(genome = "species")
  }
  
  # 1st calc
  joined_df2 <- merge(joined_df1, distinct(magic18_models[,c(1,3)]) , by.x = "gene_id", by.y = "gene_id", all.x = TRUE) %>% 
    group_by(genome) %>%
    summarise(gene_ids_in_cluster_with_common_domain = sum(count))
  
  interpro_annotated <- x[!(is.na(x$Pan.id)),] 
  
  interpro_annotated1 <- interpro_annotated[!(interpro_annotated$V12 == "-"),] 
  
  interpro_annotated_genes <- interpro_annotated1 %>% group_by(genome) %>% summarise(Total_genes_with_interpro_domain = n_distinct(gene_id)) %>% na.omit()
  
  gene_count_table <- merge(joined_df2, interpro_annotated_genes, by = "genome") %>% mutate(gene_ids_in_cluster_without_common_domain = Total_genes_with_interpro_domain - gene_ids_in_cluster_with_common_domain) %>% 
    mutate(percent_missing = (gene_ids_in_cluster_without_common_domain / Total_genes_with_interpro_domain)*100 )
  
  
  ##2nd calc
  
  merged_summary_interpro2 <- merged_summary_interpro %>% full_join(summary_most_common, by = "Pan.id")
  
  ipr_annot <- x %>% select(V12, V13) %>% distinct()
  
  merged_summary_interpro2 <- merged_summary_interpro2 %>% left_join(ipr_annot, by = "V12") %>% 
    rename(most_common_interpro_id  = "V12", most_common_Interpro_domain_name  = "V13" )
  
  
  ###rename, arrange columns etc. 
  
  merged_summary_interpro2 <- merged_summary_interpro2 %>% mutate(accessions_lacking_interpro_domain  = pangene_occupancy - accession_counts_with_most_common_interpro)
  
  merged_summary_interpro2 <- merged_summary_interpro2 %>% relocate(., most_common_interpro_id , .after = genomes_with_most_common_interpro)
  
  
  #### do the same for Pfam
  df_pfam <- x %>% filter(V4 == "Pfam") %>% select(Pan.id, V5, gene_id, genome, taxa_size) %>% distinct()
  
  # identify the most common domain in all genes of a clusters
  res_pfam <- as.data.frame(df_pfam %>%
                              group_by(Pan.id) %>% filter(!(V5 == "-")) %>% 
                              summarise(most_common_value = names(which.max(table(V5)))) %>% na.omit())
  
  print(head(res_pfam))
  
  split_list <- split(x, x$Pan.id)
  split_list_pfam <- split_list[names(split_list) %in% res_pfam$Pan.id]
  
  pb <- progress_bar$new(
    format = "  working [:bar] :percent eta: :eta",
    total = length(res_pfam$Pan.id), clear = FALSE, width= 60)
  
  filtered_list_pfam <- lapply(names(split_list_pfam), function(x) {
    pb$tick()
    filter_value <- res_pfam %>% filter(Pan.id == x) %>% pull(most_common_value)
    split_list_pfam[[x]] %>% filter(V5 == filter_value) 
  })
  
  pb <- progress_bar$new(
    format = "  working [:bar] :percent eta: :eta",
    total = length(res_pfam$Pan.id), clear = FALSE, width= 60)
  
  filtered_list1_pfam <- lapply(filtered_list_pfam, function(x) {
    pb$tick()
    x %>%
      group_by(gene_id) %>%
      summarise(count = n_distinct(gene_id))
  }) 
  
  
  pb <- progress_bar$new(
    format = "  working [:bar] :percent eta: :eta",
    total = length(res_pfam$Pan.id), clear = FALSE, width= 60)
  
  options(dplyr.summarise.inform = FALSE)
  
  filtered_list2_pfam <- lapply(filtered_list_pfam, function(x) {
    pb$tick()
    x %>%
      group_by(Pan.id, V5) %>%
      summarise(accession_counts_with_most_common_pfam = n_distinct(accession) ,
                genomes_with_most_common_pfam = n_distinct(genome) ,
                accesion_names_with_common_pfam_domain = unique(accession) %>% toString())
  })
  
  
  # Assign names back to the list
  names(filtered_list1_pfam) <- names(split_list_pfam)
  joined_df1_pfam <- bind_rows(filtered_list1_pfam) 
  
  # combine summary of all clusters and merge with previous dataframes
  
  summary_most_common_pfam <- bind_rows(filtered_list2_pfam)
  
  # 1st calc
  joined_df2_pfam <- merge(joined_df1_pfam, distinct(magic18_models[,c(1,3)]) , by.x = "gene_id", by.y = "gene_id", all.x = TRUE) %>% 
    group_by(genome) %>%
    summarise(gene_ids_in_cluster_with_common_domain = sum(count))
  
  pfam_annotated <- x[!(is.na(x$Pan.id)),] %>% filter(V4 == "Pfam")
  
  pfam_annotated1 <- pfam_annotated[!(pfam_annotated$V5 == "-"),] 
  
  
  pfam_annotated_genes <- pfam_annotated1 %>% group_by(genome) %>% summarise(Total_genes_with_pfam_domain = n_distinct(gene_id)) %>% na.omit()
  
  gene_count_table_pfam <- merge(joined_df2_pfam, pfam_annotated_genes, by = "genome") %>% 
    mutate(gene_ids_in_cluster_without_common_domain = Total_genes_with_pfam_domain - gene_ids_in_cluster_with_common_domain)
  
  gene_count_table_pfam <- gene_count_table_pfam %>% 
    mutate(percent_missing = (gene_ids_in_cluster_without_common_domain / Total_genes_with_pfam_domain)*100 )
  
  
  # 2nd calc
  
  merged_summary_pfam <- Reduce(function(x, y) merge(x, y, by = "Pan.id", all = TRUE), 
                                list(merged_summary_interpro,  summary_most_common_pfam))
  
  
  
  pfam_annot <- x %>% filter(V4 == "Pfam") %>% select(V5, V6) %>% distinct()
  
  merged_summary_pfam2 <- merge(merged_summary_pfam, pfam_annot, by = "V5", all.x = TRUE) %>% rename(most_common_pfam_id  = "V5", most_common_pfam_domain_name  = "V6")
  
  
  ###rename, arrange columns etc. 
  
  merged_summary_pfam2 <- merged_summary_pfam2 %>% mutate(accessions_lacking_pfam_domain  = pangene_occupancy - accession_counts_with_most_common_pfam)
  
  merged_summary_pfam2 <- merged_summary_pfam2 %>% relocate(., most_common_pfam_id , .after = genomes_with_most_common_pfam)
  
  
  #### list all the dataframes to be output
  
  
  list_df <- list(gene_count_table, pfam_summary, IPR_summary , merged_summary_interpro, merged_summary_interpro2, gene_count_table_pfam, merged_summary_pfam2)
  
  names(list_df) <- c("gene_count_table", "easy_pfam_summary", "easy_interpro_summary", "merged_summary_interpro_pfam", "merged_summary_interpro_extended", "gene_count_table_pfam", "merged_summary_pfam_extended")
  
  return(list_df)
  
}