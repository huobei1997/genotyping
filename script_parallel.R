#BiocManager::install("sangerseqR")
#install.packages("pbmcapply")

library(sangerseqR)
library(parallel)
library(pbmcapply)  

args <- commandArgs(trailingOnly=T)
folder_path <- "/home/ubuntu/project/test/1_Huang28086_N"
output <- "/home/ubuntu/project/test/results.xlsx"

# Get all ab1 files
ab1_files <- list.files(folder_path, pattern = "\\.ab1$", full.names = TRUE)

# Split the file list into several batches
batch_size <- 1  
file_batches <- split(ab1_files, ceiling(seq_along(ab1_files)/batch_size))

# Determine the number of available cores
num_cores <- detectCores() - 1
num_cores <- max(1, num_cores)  

# Define the batch processing function
process_batch <- function(batch_files) {
  
  batch_results <- data.frame(Filename = character(), 
                             Sequence1 = character(), 
                             Sequence2 = character(),
                             stringsAsFactors = FALSE)
  
  for (file in batch_files) {
    tryCatch({
      abif_data <- read.abif(file)
      hetcalls <- makeBaseCalls(sangerseq(abif_data))
      seq1 <- as.character(hetcalls@primarySeq)
      seq2 <- as.character(hetcalls@secondarySeq)
      
      batch_results <- rbind(batch_results, 
                           data.frame(Filename = basename(file), 
                                     Sequence1 = seq1, 
                                     Sequence2 = seq2, 
                                     stringsAsFactors = FALSE))
    }, error = function(e) {
      message("Error processing file:", file, "-", e$message)
    })
  }
  
  return(batch_results)
}

# Use pbmclapply to show a progress bar
result_list <- pbmcapply::pbmclapply(file_batches, process_batch, mc.cores = num_cores)
results <- do.call(rbind, result_list)

# Define prefix and suffix
prefix <- "TTAGCATTTG"
suffix <- "AAAGCAGTTA"

# Check the sequences and find the middle bases
for (i in 1:nrow(results)) {
  seq1 <- results$Sequence1[i]
  seq2 <- results$Sequence2[i]
  
  # Find matches for prefixes and suffixes
  start_index1 <- regexpr(prefix, seq1)
  end_index1 <- regexpr(suffix, seq1)
  
  start_index2 <- regexpr(prefix, seq2)
  end_index2 <- regexpr(suffix, seq2)
  
  # Extract middle bases
  if (start_index1 > 0 && end_index1 > 0) {
    middle_base1 <- substr(seq1, start_index1 + nchar(prefix), end_index1 - 1)
    results$middle_base1[i] <- middle_base1
  } else {
    results$middle_base1[i] <- "missing"  
  }
  
  if (start_index2 > 0 && end_index2 > 0) {
    middle_base2 <- substr(seq2, start_index2 + nchar(prefix), end_index2 - 1)
    results$middle_base2[i] <- middle_base2
  } else {
    results$middle_base2[i] <- "missing"  
  }
}

# Judge
for (i in 1:nrow(results)) {
  middle1 <- results$middle_base1[i]
  middle2 <- results$middle_base2[i]
  
  if (middle1 == "A") {
    results$judge[i] <- "b6 control"
  } else if (middle1 == "C" && middle2 == "C") {
    results$judge[i] <- "Hom"
  } else if (middle1 == "C" && middle2 == "A") {
    results$judge[i] <- "Het"
  } else {
    results$judge[i] <- "others"
  }
}

writexl::write_xlsx(results,output)

  
  
  
  
  
