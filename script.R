library(sangerseqR)

# Setting the folder path
folder_path <- "/Users/huobei/Desktop/241128_Huang26581/"

# Access to all ab1 documents
ab1_files <- list.files(folder_path, pattern = "\\.ab1$", full.names = TRUE)

results <- data.frame(Filename = character(), Sequence1 = character(), Sequence2 = character(),stringsAsFactors = FALSE)

# Read ab1 files and extract sequences
for (file in ab1_files) {
  abif_data <- read.abif(file)
  hetcalls <- makeBaseCalls(sangerseq(abif_data))  
  seq1 <- as.character(hetcalls@primarySeq)  
  seq2 <- as.character(hetcalls@secondarySeq)
  results <- rbind(results, data.frame(Filename = basename(file), Sequence1 = seq1,Sequence2 = seq2, stringsAsFactors = FALSE))
}

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

writexl::write_xlsx(results,"/Users/huobei/Desktop/results.xlsx")
