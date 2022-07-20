#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# THIS SCRIPT WILL PUT ALL ALLELES WITHIN EACH GROUP IN NUMERIC/ALPHABETICAL ORDER

# set-up:
library(tidyr)
library(stringr)
library(tidyverse)
filepath = "/Users/nicholasminor/Documents/dholk_experiments/27453/3b_miSeq-trimming" # replace this filepath with wherever you have stored the fastas
setwd(filepath)
fasta <- read.delim("ipd-mhc-nhp-2022-07-11_cleaned.miseq.trimmed.deduplicated.fasta", header = F)
# fasta <- read.delim(args[1], header = F)

file_basename <- "ipd-mhc-nhp-2022-07-11_cleaned.miseq.trimmed.deduplicated"
# file_basename <- str_remove(args[1], ".fasta")

# how many rows of allele groups are there?
paste("There are", length(fasta[grepl("|", fasta[,1], fixed=T),1]), "rows to examine in this fasta.", sep = " ")

# correcting allele order in unsorted allele group rows
rows_fixed <- NA
for (i in 1:nrow(fasta)){
  if (grepl(">", fasta[i,1])){
    if (grepl("|", fasta[i,1], fixed=T)){
      sub <- fasta[i,1]
      
      # PRE-SORTING PROCESSING:
      sub <- str_replace_all(sub, fixed(">"), "") %>%
      str_replace_all(fixed("|"), " ") %>%
      strsplit(split = " ") %>%
      unlist()
      sub <- t(t(sub))
      
      # SORTING:
      sub_sorted <- sort(sub[,1])
      
      # POST-SORTING PROCESSING:
      sub_sorted <- paste(sub_sorted, collapse = "|", sep = "")
      sub_sorted <- paste(">", sub_sorted, sep = "")
      if (sub_sorted!=fasta[i,1]){
        
      # INSERTING BACK INTO FASTA:
        fasta[i,1] <- sub_sorted
        rows_fixed <- c(na.omit(rows_fixed), i)
        
      } # if-else statement 3
    }  # end if-else statement 2
  }  # end if-else statement 1
} # end for loop

paste("This script has fixed", length(rows_fixed), "rows in the inserted fasta.", sep = " ")

# exporting now-sorted fasta
write.table(fasta, 
            paste(file_basename, ".sorted.fasta", sep = ""),
            sep = "\t", quote = F, col.names = F, row.names = F)

