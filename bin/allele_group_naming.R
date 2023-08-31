#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# THIS SCRIPT WILL RENAME THE ALLELE GROUPS FOR EACH CLASS OF ALLELES IN THE IPD REFERENCE FASTA

# set-up:
library(tidyr)
library(stringr)

fasta <- read.delim(args[1], header = F)
file_basename <- str_remove(args[1], ".fasta")

main <- compiler::cmpfun(function(fasta, filebasename){
  
  if (grepl("-nhp-", file_basename)){
    
    animal_filename <- "nhp"
    
  } else if (TRUE %in% grepl(">Mamu-", fasta[,1])){
    
    animal <- ">Mamu-"
    animal_filename <- "mamu"
    
  } else if (TRUE %in% grepl(">Mafa-", fasta[,1])){
    
    animal <- ">Mafa-"
    animal_filename <- "mafa"
    
  } else if (TRUE %in% grepl(">Mane-", fasta[,1])){
    
    animal <- ">Mane-"
    animal_filename <- "mane"
    
  }
  
  
  if (animal_filename=="nhp"
      | length(unique(substring(fasta[grepl("|", fasta[,1], fixed = T),1], 2,5))) > 1
  ){
    
    print("Error: Allele groups cannot be assigned for multi-species reference files without introducing further ambiguity.")
    
    errline <- data.frame("V1"="Error: Allele groups cannot be assigned for multi-species reference files without introducing further ambiguity.")
    fasta <- rbind(errline, fasta)
    
    write.table(fasta, 
                paste(file_basename, ".NOT_RENAMED.fasta", sep = ""),
                sep = "\t", quote = F, col.names = F, row.names = F)
    
  } else {
    
    fasta[grepl("|", fasta[,1], fixed = T),] <- str_remove_all(fasta[grepl("|", fasta[,1], fixed = T),], "Mamu-") %>%
      str_remove_all("Mafa-") %>%
      str_remove_all("Mane-") %>%
      str_replace_all(fixed("|"), ",")
    
    groups <- as.data.frame(fasta[grepl(",", fasta[,1], fixed = T),])
    rows_renamed <- NA
    rows_to_rename <- nrow(groups)
    total_rows <- nrow(fasta)
    
    
    # A1 ALLELES ####
    # ---------
    name_rows <- groups[grepl("-A1_", groups[,1]) | grepl(">A1_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 6)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    
    # A2 ALLELES ####
    # ---------
    name_rows <- groups[grepl("-A2_", groups[,1]) | grepl(">A2_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 5)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # A3 ALLELES ####
    # ---------
    name_rows <- groups[grepl("-A3_", groups[,1]) | grepl(">A3_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 5)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # A4 ALLELES ####
    # ---------
    name_rows <- groups[grepl("-A4_", groups[,1]) | grepl(">A4_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 5)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # A5 ALLELES ####
    # ---------
    name_rows <- groups[grepl("-A5_", groups[,1]) | grepl(">A5_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 5)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # A6 ALLELES ####
    # ---------
    name_rows <- groups[grepl("-A6_", groups[,1]) | grepl(">A6_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 5)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # A7 ALLELES ####
    # ---------
    name_rows <- groups[grepl("-A7_", groups[,1]) | grepl(">A7_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows)>0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 5)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # AG ALLELES ####
    # ---------
    name_rows <- groups[grepl("-AG", groups[,1]) | grepl(">AG", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 6)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # B02Ps_01 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">B02Ps", groups[,1]) | grepl("-B02Ps", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 8)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # B10Ps_01 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">B10Ps_01", groups[,1]) | grepl("-B10Ps_01", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 8)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # B11L_01 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">B11L_01", groups[,1]) | grepl("-B11L_01", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 7)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # B14Ps_175 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">B14Ps_175", groups[,1]) | grepl("-B14Ps_175", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 9)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # B16_01 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">B16_01", groups[,1]) | grepl("-B16_01", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 6)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # B17_01 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">B17_01", groups[,1]) | grepl("-B17_01", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 6)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # B21Ps_01 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">B21Ps_01", groups[,1]) | grepl("-B21Ps_01", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 8)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # B ALLELES ####
    # ---------
    name_rows <- groups[grepl(">B_", groups[,1]) | grepl("-B_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 5)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # DPA1 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">DPA1", groups[,1]) | grepl("-DPA1_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 7)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # DPB1 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">DPB1", groups[,1]) | grepl("-DPB1_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 7)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # DQA1 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">DQA1", groups[,1]) | grepl("-DQA1_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 7)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # DQB1 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">DQB1", groups[,1]) | grepl("-DQB1_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 7)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # DRB1 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">DRB1", groups[,1]) | grepl("-DRB1_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 7)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # DRB3 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">DRB3", groups[,1]) | grepl("-DRB3_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 7)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # DRB_W0 ALLELES ####
    # ---------
    name_rows <- groups[grepl(">DRB_", groups[,1]) | grepl("-DRB_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 8)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # E ALLELES ####
    # ---------
    name_rows <- groups[grepl(">E_", groups[,1]) | grepl("-E_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 4)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # F ALLELES ####
    # ---------
    name_rows <- groups[grepl(">F_", groups[,1]) | grepl("-F_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 4)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # G ALLELES ####
    # ---------
    name_rows <- groups[grepl(">G_", groups[,1]) | grepl("-G_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 4)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # I ALLELES ####
    # ---------
    name_rows <- groups[grepl(">I_", groups[,1]) | grepl("-I_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 4)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    # J ALLELES ####
    # ---------
    name_rows <- groups[grepl(">J_", groups[,1]) | grepl("-J_", groups[,1]),] %>%
      str_replace_all(fixed(">"), "")
    
    if (length(name_rows) > 0){
      name_rows <- data.frame("group" = NA, "alleles" = name_rows)
      
      for (i in 1:nrow(name_rows)){
        sub <- name_rows[i,"alleles"]
        group <- substr(sub, 1, 4)
        name_rows[i, "group"] <- group
      }
      
      for (i in unique(name_rows$group)){
        sub <- name_rows[name_rows$group==i,]
        if (nrow(sub)>1){
          for (j in 1:nrow(sub)){
            sub[j,"group"] <- paste(sub[j,"group"], "g", j, sep = "")
          } # end for loop 2
        } else {
          sub[,"group"] <- paste(sub[,"group"], "g", sep = "")
        } # end if-else statement
        for (k in unique(sub$alleles)){
          name_rows[name_rows$alleles==k,"group"] <- sub[sub$alleles==k, "group"]
        } # end for loop 3
      } # end for loop 1
      
      final_names <- name_rows
      
      for (i in 1:nrow(name_rows)){
        name_rows[i,2] <- paste(">", name_rows[i,2], sep = "")
      }
      
      for (i in 1:nrow(final_names)){
        final_names[i,] <- paste(final_names[i,], collapse = "|", sep = "")
        final_names[i,] <- paste(animal, final_names[i,], sep = "")
      }
      final_names$alleles <- NULL
      colnames(final_names) <- NULL
      
      for (i in 1:nrow(final_names)){
        fasta[fasta[,1]==name_rows$alleles[i],] <- final_names[i,1]
        rows_renamed <- c(na.omit(rows_renamed), i)
      }
      progress <- paste(round((100 * (length(rows_renamed)/rows_to_rename)), 2), 
                        "% finished renaming the IPD reference database.", sep = "")
    }
    
    
    # SANITY CHECK(S) ####
    # progress
    groups <- as.data.frame(fasta[grepl("|", fasta[,1], fixed = T),]) # visually check groups to make sure that not groups are without a name
    groups_final <- data.frame("group" = str_replace_all(groups[,1], fixed(">"), ""))
    groups_final$group <- str_replace(groups_final$group, fixed("|"), " ")
    groups_final <- separate(groups_final, col = group, into = c("group", "alleles"), sep = " ")
    
    
    # FINAL EXPORTING ####
    # --------------
    write.table(fasta, 
                paste(file_basename, ".renamed.fasta", sep = ""),
                sep = "\t", quote = F, col.names = F, row.names = F)
    write.csv(groups_final, 
              paste("IPD_", animal_filename, "_reference_allele_groups.csv", sep = ""), 
              row.names = F)
  }
  
})

main(fasta, file_basename)
