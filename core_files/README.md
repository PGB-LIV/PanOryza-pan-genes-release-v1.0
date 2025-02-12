The folder contains intermediate files that were either obtained from various databases, collaborators or created during the analyses. The other files can be created or retrieved as follows :

-  Clusters sequence summary : using [create_cluster_sum](scripts/create_cluster_sum.R) and [read_parse_clusters_summary](scripts/read_parse_clusters_summary.R)
-  magic18_models.tab : contains gene and protein ids of all the 18 genomes used in the analyses. These can be obtained from ensembl plants (https://plants.ensembl.org/info/data/ftp/index.html). 
-  oryza_sativa_from_homologies.tsv : contains homology relationships for species oryza_sativa and can be filtered from Compara.111.protein_default.homologies.tsv available at ensembl compara (https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-58/tsv/ensembl-compara/homologies) . 
-  magic18_interpro_clusters : contains results of InterproScan-v5.62-94.0 run on all genomes along with the pan-gene identifiers.
-  peptide_mapping.tsv : The results are available at https://peptideatlas.org/builds/rice/ . 
