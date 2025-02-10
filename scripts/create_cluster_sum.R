suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(progress))

setwd("Os4530.POR.1/oryzasativanipponbaremerged/") # Os4530.POR.1 (version 1) can be obtained from zenodo and extracted from archive

destination_folder <- "Os4530.POR.1/oryzasativanipponbaremerged"

# List of FASTA files
fasta_files <- as.character(list.files(destination_folder, pattern = "*.cds.faa"))


pb <- progress_bar$new(
  format = "  working [:bar] :percent eta: :eta",
  total = length(fasta_files), clear = FALSE, width= 60)

calculate_protein_lengths <- function(fasta_files) {
  data <- data.frame(File = character(), Sequence_ID = character(), Length = integer(), stringsAsFactors = FALSE)
  
  for (fasta_file in fasta_files) {
    pb$tick()
    sequences <- readAAStringSet(fasta_file)
    lengths <- width(sequences)
    ids <- names(sequences)
    
    file_data <- data.frame(File = basename(fasta_file), Sequence_ID = ids, Length = lengths)
    data <- bind_rows(data, file_data)
  }
  
  return(data)
}


# Calculate lengths and tabulate results
clusters_sequence_summary <- calculate_protein_lengths(fasta_files) ## Takes 1-2 hours


write.table(clusters_sequence_summary, "clusters_sequence_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
