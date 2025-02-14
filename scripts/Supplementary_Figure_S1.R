# Supplement Fig 1
## check QC of pan-genes based on the output of get_pangenes pipeline and all 3 Nipponbare annotations therein. 

suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))

load(file = "core_workspace.RData") #  load the core set of files

MSA_OsNip <- read.delim("core_files/magic clusters - magic15MSA1stisof.OsNip.tsv")
log.quality.cDNA.fix <- read.delim("core_files/log.quality.cDNA.fix.tsv", header=FALSE)
pep_global_dist.max <- read.delim("core_files/pep_global_dist.max.tsv", header=FALSE)


## add pan-genes notations

MSA_OsNip$ids <- MSA_OsNip$file

MSA_OsNip$ids <- str_replace_all(MSA_OsNip$ids, ".cds.faa", "")

MSA_OsNip$ids <- str_replace_all(MSA_OsNip$ids, "gene:", "")

MSA_OsNip <- MSA_OsNip %>% left_join(clusters_summary[,c(1,6)], join_by(ids == merged.or.gene_id))

MSA_OsNip_final <- MSA_OsNip %>% select(-ids, -file) %>% relocate(Pan.id, .before = X1stisof) 

# remove NA 
MSA_OsNip_final <- MSA_OsNip_final %>% na.omit() %>% distinct()
pep_global_dist.max <- pep_global_dist.max %>% na.omit() %>% distinct()
log.quality.cDNA.fix <- log.quality.cDNA.fix %>% na.omit() %>% distinct()


plot_A <- MSA_OsNip_final %>% 
  ggplot(., aes(x = SE_exons )) +
  geom_histogram(position = "identity", alpha = 0.5, fill = "#56B4E9", color = "black") + theme_cowplot() + 
  labs(x = "SE exons number", y = "Counts of pan-genes", fill = "", 
       subtitle = paste("n =", length(unique(MSA_OsNip_final$Pan.id)))) +
  scale_fill_manual(values = "blue" )

plot_B <- MSA_OsNip_final %>% 
  ggplot(., aes(x = SE_len )) +
  geom_histogram(position = "identity", alpha = 0.5, fill = "#56B4E9", color = "black") + theme_cowplot() + 
  labs(x = "SE length", y = "Counts of pan-genes", fill = "", 
       subtitle = paste("n =", length(unique(MSA_OsNip_final$Pan.id)))) +
  scale_fill_manual(values = "blue" )

plot_C <- MSA_OsNip_final %>% 
  ggplot(., aes(x = SE_dist )) +
  geom_histogram(position = "identity", alpha = 0.5, fill = "#56B4E9", color = "black") + theme_cowplot() + 
  labs(x = "SE distance", y = "Counts of pan-genes" , fill = "", 
       subtitle = paste("n =", length(unique(MSA_OsNip_final$Pan.id)))) +
  scale_fill_manual(values = "blue" )


plot_D <- MSA_OsNip_final %>% 
  ggplot(., aes(x = Ca )) +
  geom_histogram(position = "identity", alpha = 0.5, fill = "#56B4E9", color = "black") + theme_cowplot() + 
  labs(x = "MSA completeness" , y = "Counts of pan-genes", fill = "", 
       subtitle = paste("n =", length(unique(MSA_OsNip_final$Pan.id)))) +
  scale_fill_manual(values = "blue" )

plot_E <- pep_global_dist.max %>% 
  ggplot(., aes(x = V2 )) +
  geom_histogram(position = "identity", alpha = 0.5, fill = "#56B4E9", color = "black") + theme_cowplot() + 
  labs(x = "maximum distance in\n protein alignment" , y = "Count of pan-genes", fill = "", subtitle = paste("n =", length(unique(pep_global_dist.max$V1)))) +
  scale_fill_manual(values = "blue" )

plot_F <- log.quality.cDNA.fix %>% 
  ggplot(., aes(x = V2 )) +
  geom_histogram(position = "identity", alpha = 0.5, fill = "#56B4E9", color = "black") + theme_cowplot() + 
  labs(x = "Minimum BLASTN cDNA \n alignment identity(in %)" , y = "Count of pan-genes", fill = "", subtitle = paste("n =", length(unique(log.quality.cDNA.fix$V1)))) +
  scale_fill_manual(values = "blue" )


# uncomment below to export
# svg("supplement_QC.svg", width=9,height=6, pointsize = 8)
# ggarrange(plot_A, plot_B, plot_C, plot_D , plot_E,plot_F, labels = "AUTO")
# dev.off()


# Supplementary table
# write.table(MSA_OsNip_final, file = "Supplementary table S2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

