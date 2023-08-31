#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# THIS SCRIPT WILL PUT ALL ALLELES WITHIN EACH GROUP IN NUMERIC/ALPHABETIC ORDER

# set-up:
library(tidyverse)

# pull arguments from the command line by default
fasta <- read.delim(args[1], header = FALSE)
file_basename <- str_remove(args[1], ".fasta")

main <- compiler::cmpfun(function(sequences, basename) {

  # how many rows of allele groups are there?
  paste("There are", length(sequences[grepl("|", sequences[, 1], fixed = TRUE), 1]),
        "rows to examine in this fasta.", sep = " ")

  # correcting allele order in unsorted allele group rows
  rows_fixed <- NA
  for (i in 1:nrow(sequences)) {

    if (!grepl(">", sequences[i, 1])) {
      next
    }
    if (!grepl("|", sequences[i, 1], fixed = TRUE)) {
      next
    }
    sub <- sequences[i, 1]

    # PRE-SORTING PROCESSING:
    sub <- str_replace_all(sub, fixed(">"), "") %>%
      str_replace_all(fixed("|"), " ") %>%
      strsplit(split = " ") %>%
      unlist()
    sub <- t(t(sub))

    # SORTING:
    sub_sorted <- sort(sub[, 1])

    # POST-SORTING PROCESSING:
    sub_sorted <- paste(sub_sorted, collapse = "|", sep = "")
    sub_sorted <- paste(">", sub_sorted, sep = "")
    if (sub_sorted == sequences[i, 1]) {
      next
    }

    # INSERTING BACK INTO FASTA:
    sequences[i, 1] <- sub_sorted
    rows_fixed <- c(na.omit(rows_fixed), i)

  } # end for loop

  paste("This script has fixed", length(rows_fixed),
        "rows in the inserted fasta.", sep = " ")

  # exporting now-sorted fasta
  write.table(sequences,
              paste(file_basename, ".sorted.fasta", sep = ""),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

})

# run main
main(fasta, file_basename)
