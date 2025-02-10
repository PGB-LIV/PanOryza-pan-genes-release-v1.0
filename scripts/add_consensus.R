## add consensus column of yes / no bases on inputs in 'colname' and the pan-genes dataframe grouped by group_name(pan_gene, usually Pan.id or oryzasativanipponbaremerged)

add_consensus <- function(df, group_name, colname) {
  
  df_list <- split(df, df[[group_name]])
  
  df_list2 <- lapply(df_list, function(x){ 
    split(x, x$accession)
  } )
  
  df_list3 <- unlist(df_list2, recursive = FALSE)
  
  df_list3 <- lapply(df_list3, function(x) { 
    x %>% mutate(.,consensus = ifelse(any(x[[colname]] == "yes"), "yes", "no"))
  })
  
  return(bind_rows(df_list3))
  
}