###### read and parse clusters summary 

read_parse_clusters_summary <- function(df){ ## df includes path to the summary df
  df <- read.delim(df)
  
  # separate out header ID and rename
  df <- df %>% separate_wider_delim(., cols = Sequence_ID, delim = " " , names_sep = "") %>% 
    dplyr::rename(transcript = "Sequence_ID1" , merged.or.gene_id = "Sequence_ID2" , region = "Sequence_ID3", source = "Sequence_ID4" )
  
  #remove brackets
  df$source <- gsub(pattern = "\\[|\\]", replacement = ""  , x = df$source)
  
  #remove prefixes to gene nand transcript names
  df$transcript <- str_replace_all(df$transcript , "transcript:", "")
  df$merged.or.gene_id <- str_replace_all(df$merged.or.gene_id , "gene:", "")
  
  df$genome <- df$source
  
  #separate nipponbare and rest
  df_nip <- df[df$source == "oryza_sativa_nipponbaremerged",]
  df_rest <- df[!(df$source == "oryza_sativa_nipponbaremerged"),]
  
  df_nip$genome <- ifelse(grepl("^LOC", df_nip$transcript), "Nipponbare(MSU)", 
                          ifelse(grepl("OsNip", df_nip$transcript), "Nipponbare(Gramene)",
                                 ifelse(grepl("Os", df_nip$transcript), "Nipponbare(RAPDB)", NA)))
  
  
  df_rest <- df_rest %>% 
    mutate(., genome = recode(genome, "oryza_sativa_gobolsailbalam" = "Gobol Sail Balam", 
                              "oryza_sativa_liuxu_chr" = "Liu Xu" , 
                              "oryza_sativa_khaoyaiguang" = "Khao Yai Guang",
                              "oryza_sativa_larhamugad" = "Larha Mugad",
                              "oryza_sativa_mh63" = "Minghui 63",
                              "oryza_sativa_natelboro" = "Natel Boro",
                              "oryza_sativa_ir64" = "IR64" ,
                              "oryza_sativa_arc" = "ARC 10497" ,
                              "oryza_sativa_ZS97" = "Zhenshan 97" ,
                              "oryza_sativa_pr106" = "PR106" ,
                              "oryza_sativa_lima" = "Lima" ,
                              "oryza_sativa_ketannangka" = "Ketan Nangka" ,
                              "oryza_sativa_chaomeo" = "Chao Meo" ,
                              "oryza_sativa_azucena" = "Azucena" ,
                              "oryza_sativa_n22" = "N22"))
  
  df_final <- bind_rows(df_nip, df_rest)
  
  df_final <- df_final %>% mutate(Pan.id = str_replace_all(File, ".cds.faa", ""))
  
  return(df_final)
  
}


# usage : 
# clusters_summary_path <- # REPLACE WITH THE PATH TO "clusters_sequence_summary.tsv"
# clusters_summary <- read_parse_clusters_summary(clusters_summary_path)
