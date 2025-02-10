
# Extended data
# Supplementary figure S3 

# recommended to run analyses in Fig 4 before plotting this one

interpro2go <- read.table("core_files/interpro2go.txt", sep=";", quote="", comment.char="")
interpro2go_split <- interpro2go %>% separate_wider_delim(., V1, delim = " > ", names = c("interpro_name", "GO_term"))

interpro2go_split$interpro_name <-  str_replace_all(interpro2go_split$interpro_name, "InterPro:" , "")

interpro2go_split$V2 <-  str_replace_all(interpro2go_split$V2, " " , "")

interpro2go_split <- interpro2go_split %>% separate_wider_delim(., interpro_name , " ", names = c("id", "name"), too_many  = "merge")

ipr_clusters2go <- merge(ipr_clusters , interpro2go_split , by.x = "V12", by.y = "id")

go_term2gene <- ipr_clusters2go %>% select(V2, Pan.id)
go_term2name <- interpro2go_split[,c(4,3)]

mean_occupancy_interpro_go <- mean_occupancy_interpro %>% left_join(interpro2go_split[,-2] , by = "id")

# filter for at least 5 pan-genes ? 

mean_occupancy_interpro_go_unfil <- mean_occupancy_interpro_go %>% filter(!(count_pangenes > 5)) 

mean_occupancy_interpro_go <- mean_occupancy_interpro_go %>% filter(count_pangenes > 5) 


# domain groups by occupancy


highly_variable = mean_occupancy_interpro_go %>% filter(mean_occupancy < 10)

partially_variable <-  mean_occupancy_interpro_go %>% filter(mean_occupancy > 10 & mean_occupancy < 15)

completely_invariable <- mean_occupancy_interpro_go %>% filter(mean_occupancy > 15)


# get pan-genes with these domains

highly_variable_pan_id <- domain_occupancy_table_interpro %>% filter(V12 %in% highly_variable$id) %>% ungroup() %>% dplyr::select(Pan.id)

partially_variable_pan_id <- domain_occupancy_table_interpro %>% filter(V12 %in% partially_variable$id) %>% ungroup() %>% dplyr::select(Pan.id)

completely_invariable_pan_id <- domain_occupancy_table_interpro %>% filter(V12 %in% completely_invariable$id) %>% ungroup() %>% dplyr::select(Pan.id)


# run enricher function for over-representation analyses

go_hv_domain <- enricher(highly_variable_pan_id$Pan.id, TERM2GENE = go_term2gene , TERM2NAME = go_term2name , pvalueCutoff = 0.05, pAdjustMethod = "BH") # highly variable

go_pv_domain <- enricher(partially_variable_pan_id$Pan.id, TERM2GENE = go_term2gene , TERM2NAME = go_term2name , pvalueCutoff = 0.05, pAdjustMethod = "BH") # partially variable

go_ci_domain <- enricher(completely_invariable_pan_id$Pan.id, TERM2GENE = go_term2gene , TERM2NAME = go_term2name , pvalueCutoff = 0.05, pAdjustMethod = "BH") # completely invariable


# Plot indiviodual groups

highly_variable_plot <- dotplot(go_hv_domain, 
                                color = "p.adjust", 
                                showCategory = 20, 
                                size = "GeneRatio",
                                font.size = 8,
                                title = "highly variable")

partially_variable_plot <- dotplot(go_pv_domain, 
                                   color = "p.adjust", 
                                   showCategory = 20, 
                                   size = "GeneRatio",
                                   font.size = 8,
                                   title = "partially variable")

invariable_plot <- dotplot(go_ci_domain, 
                           color = "p.adjust", 
                           showCategory = 20, 
                           size = "GeneRatio",
                           font.size = 8,
                           title = "Invariable(conserved)")


# optional uncomment to export image
# svg("./plots/GO_variable_plots_new.svg", width = 15, height = 7)
# ggarrange(highly_variable_plot, partially_variable_plot, invariable_plot , nrow = 1, ncol = 3)
# dev.off()


#### exploring top GO terms

# for GO:cell surface receptor signaling pathway , the IPR036537 and IPR045274 are the two interpro terms in our set.

receptor_pan_genes <-  magic18_interpro_clusters %>% filter(grepl("IPR045274", V12) | grepl("IPR036537", V12))  %>% select(Pan.id) %>% distinct()

#about 100 other domains also present in these pan-genes
receptor_iprs <- magic18_interpro_clusters %>% filter(Pan.id %in% receptor_pan_genes$Pan.id) %>% select(V12, V13) %>% distinct() %>% na.omit()

receptor_iprs_mean_domain_occupancy <- receptor_iprs %>% filter(!(V12 == "-")) %>% left_join(mean_occupancy_interpro, join_by(V12 == id))

receptor_iprs_mean_domain_occupancy %>% filter(mean_occupancy < 10 & count_pangenes > 10)






