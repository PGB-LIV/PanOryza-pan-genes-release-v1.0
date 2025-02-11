# PanOryza-pan-genes-release-v1.0


This repository hosts the code for recreating analyses in the PanOryza manuscript. Code for GET_PANGENES available from: <https://github.com/Ensembl/plant-scripts/blob/master/pangenes/> The input files( .fasta and .gff format) for running GET_PANGENES are available from zenodo (10.5281/zenodo.14772953). Else, the output files for Os4530.POR.1 (version 1.0) are also available at the zenodo repository and can be used for various downstream analyses of the pan-genes using the code available here.

### To reproduce the entire analyses starting with the GET_PANGENES result, prepare various tables and intermediate files to recreate manuscript figures.
Output of get_pangenes using RPRP (MAGIC-16 accessions) as input gives out the following set of files: 

1) .cluster_list   ---> parsed in tabular format using function [parse_clusters](scripts/parse_clusters.R) --> output table named as "df_merged"
2) .matrix_genes.tr.tab --> read directly as table named "pangene_list"
3) .matrix.tr.tab 
4) Individual clusters inside folder 'oryzasativanipponbaremerged' --> *.cds.faa files of clusters used to calculate and summarise clusters and individual protein lengths
   clusters sequence summary can be created in R using [create_cluster_sum](scripts/create_cluster_sum.R). NOTE: There are several ways to do this using a Linux terminal.
   The resulting clusters sequence summary can be further parsed into a dataframe using [read_parse_clusters_summary](scripts/read_parse_clusters_summary.R)

Additional "cluster_merged" named table used at various places, created by combining "pangene_list" and "df_merged"

Interproscan tabular results for magic18 sequences were merged with the cluster files above. Recommended to load the workspace core_workspace.RData in R/Rstudio that will also load these Interproscan results for pan-genes. 
Else, [core_files.R](scripts/core_files.R) can be used to input all these files needed for downstream analysis. core_workspace.RData can be obtained from ????? zenodo ??? 

To repoduce the figure-wise analysis, please refer to the [scripts folder](scripts/)

- [Figure 1](scripts/Figure_1.R)

- [Figure 2](scripts/Figure_2.R)

- [Figure 3](scripts/Figure_3.R)

- [Figure 4](scripts/Figure_4.R)

- [Figure 5 Shiny app](heatmap_app/)
