# Supplementary figure S4
# Example of enriched domain in cloud pan-genes

suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(progress))
suppressMessages(library(cowplot))
suppressMessages(library(purrr))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(ggrepel))


load(file = "core_workspace.RData") #  load the core set of files


# COLOR SCHEME
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

color_scheme <- read.delim("./core_files/color_scheme.txt")

color_scheme <- color_scheme %>% mutate(sample = recode(sample, 
                                                        "Nipponbare(merged)" = "Nipponbare"))


# number of genes in each accession with the domain PF05132/IPR007811 # removed MSU and RAPDB for fair comparison

# Uncomment to export the plot
# png("../plots/Rpc4_genes_per_genome.png" , height  = 6000, width = 8000, res = 600 , pointsize = 14)
magic18_interpro_clusters %>% 
  filter(V5 == "PF05132" | V12 == "IPR007811") %>% 
  filter(!(genome %in% c("Nipponbare(MSU)" , "Nipponbare(RAPDB)"))) %>%
  group_by(genome) %>% summarise(gene_counts = n_distinct(gene_id)) %>% arrange(gene_counts) %>%
  left_join(color_scheme[,c(1,4)], by = c("genome" = "genome")) %>%
  ggplot() + 
  geom_col(aes(x = reorder(genome,gene_counts), y = gene_counts , fill = subspecies), position = "stack") + theme_cowplot() +
  theme(axis.text.x = element_text(size = 12,family  = "arial", angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = cbPalette) + labs(y = "Number of genes containing the domain", x = "Genome", subtitle = "PF05132/IPR007811" , fill = "Subpopulation"   )
# dev.off()

# number of pan-genes for each accession containing a gene with the domain PF05132/IPR007811

# Uncomment to export the plot
# png("../plots/Rpc4_pan_genes_per_genome.png" , height  = 6000, width = 8000, res = 600 , pointsize = 14 )

magic18_interpro_clusters %>% 
  filter(V5 == "PF05132" | V12 == "IPR007811") %>% 
  filter(!(genome %in% c("Nipponbare(MSU)" , "Nipponbare(RAPDB)"))) %>%
  group_by(genome) %>% summarise(pan_gene_counts = n_distinct(Pan.id)) %>% arrange(pan_gene_counts) %>%
  left_join(color_scheme[,c(1,4)], by = c("genome" = "genome")) %>%
  ggplot() + 
  geom_col(aes(x = reorder(genome,pan_gene_counts), y = pan_gene_counts , fill = subspecies), position = "stack") + theme_cowplot() +
  theme(axis.text.x = element_text(size = 12,family  = "arial", angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = cbPalette) + labs(y = "Number of pan-genes", x = "Genome", subtitle = "PF05132/IPR007811" , fill = "Subpopulation"  )

# dev.off()

## numbers

magic18_interpro_clusters %>% 
  filter(V5 == "PF05132" | V12 == "IPR007811") %>% 
  filter(!(genome %in% c("Nipponbare(MSU)" , "Nipponbare(RAPDB)"))) %>%
  group_by(taxa_size) %>% summarise(pan_gene_counts = n_distinct(Pan.id)) %>% 
  arrange(taxa_size)



