
## if output of get_pangenes 

clusters <- readLines(paste0(clusters_dir,"/oryzasativanipponbaremerged.cluster_list")) # clusters_dir is the path to clusters

parse_clusters <- function(clusters){
  Matrix <- do.call(rbind,strsplit(clusters, " ", fixed = TRUE))
  
  cluster_name = list()
  size = list()
  taxa = list()
  cdnafile = list()
  cdsfile = list()
  pepfile = list()
  df = data.frame()
  
  for (i in 1:length(Matrix[,1])) {
    if(Matrix[[i,1]] %in% "cluster") {
      cluster_name[[i]] = Matrix[[i,2]]
      size[[i]] = Matrix[[i,3]]
      taxa[[i]] = Matrix[[i,4]]
      cdnafile[[i]] = Matrix[[i,7]]
      cdsfile[[i]] = Matrix[[i,9]]
      pepfile[[i]] = Matrix[[i,11]]
    }
  }
  
  members  = list()
  Csize = list()
  
  Matrix1 <- gsub(pattern  = 'size=', replacement = '', Matrix)
  
  
  for (i in 1:length(Matrix1[,1])) {
    if(Matrix1[[i,1]] %in% "cluster") {
      Csize[[i]] = Matrix1[i,3]
      members[[i]] = Matrix1[as.integer(i+1):(as.integer(Csize[[i]]) + i) , 2]
    }
  }
  
  df1 = do.call(rbind, Map(data.frame, A=cluster_name, B=size , C = taxa , cdna = cdnafile , cds = cdsfile , pep = pepfile))
  df1 <- df1 %>% separate(B, c('First', 'Cluster_size'))
  df1 <- df1 %>% separate(C, c('First', 'taxa_size'))
  df1 <- df1[,-3]
  df1$gene_id <- df1$A
  df1$gene_id <- str_replace_all(df1$gene_id, 'gene:', '')
  
  df2 = do.call(rbind, Map(data.frame, A=cluster_name, B=size , C = taxa , cdna = cdnafile , cds = cdsfile , pep = pepfile, csize = Csize, members = members))
  df3 = aggregate(members ~ A + B + C + cdna + cds + pep + csize , df2, FUN = paste, collapse=",") # taxa collapsed together with , in between. 
  df_merged <- merge(df1, df3, by.x = "A" , by.y = "A")
  df_merged <- df_merged[,c("A","gene_id", "Cluster_size", "taxa_size","cdna.y", "cds.y", "pep.y","members" )]
  df_merged$taxa_size <- as.numeric(df_merged$taxa_size)
  df_merged$Cluster_size <- as.numeric(df_merged$Cluster_size)
  return(df_merged)
}