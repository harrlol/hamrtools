library(data.table)
library(tidyr)
library(dplyr)

args=commandArgs(trailingOnly=TRUE)

bed2modtbl <- function(bed) {
  out <- bed%>%
    mutate(start=bp, end =bp, seqname=as.character(chr), score=".", Name=as.character(pred.mod), strand=as.character(strand))%>%
    select(seqname, start, end, Name, score, strand)
  return(out)
}

comb.intersect <- function(bed1, bed2) {
  bed1.mods = bed2modtbl(bed1)
  bed2.mods = bed2modtbl(bed2)
  consensus <- intersect(bed1.mods, bed2.mods)
  return(consensus)
}

comb.union2 <- function(bedlist) {
  #takes the first bed to be the basis of future extensions
  out <- data.frame(bedlist[[1]])
  #for each bed after the first one, find the difference between it and what's in out already, and add to out
  for (bed in bedlist[-1]) {
    temp <- data.frame(bed)
    new <- setdiff(temp, out)
    out <- rbind(out, new)
  }
  return(out)
}

findConsensus <- function(in_dir, out_dir) {
  file_names <- list.files(path = in_dir)
  file_names <- grep(".mods.txt", file_names,value=TRUE)
  processed_variables <- c()
  for (file_name in file_names) {
    if (!(file_name %in% processed_variables)) {
      # Extract the common part of the variable name
      common_part <- sub("_[^_]*$", "", file_name)
      # Find variables with the same common part
      variables_to_process <- grep(paste0("^", common_part, "_"), file_names, value = TRUE)
      
      # Add the selected files to processed variables
      processed_variables <- c(processed_variables, variables_to_process)
      
      # If only 1 rep, then that rep is the consensus
      if (length(variables_to_process)==1) {
        consensus <- bed2modtbl(fread(file.path(in_dir, variables_to_process)))
        write.table(consensus, paste0(out_dir, "/", common_part, ".bed"), sep='\t', row.names=F, col.names=F, quote=F)
      } else {
        # If >1 rep, find combinations
        # Generate a matrix of 2 combinations of the files
        variable_combinations <- combn(variables_to_process, 2)
        
        # Initiate list to contain all intersections
        intersect_list <- list()
        for (i in 1:ncol(variable_combinations)) {
          var1 <- assign(variable_combinations[1, i], fread(file.path(in_dir, variable_combinations[1, i]),stringsAsFactors = TRUE)) 
          var2 <- assign(variable_combinations[2, i], fread(file.path(in_dir, variable_combinations[2, i]),stringsAsFactors = TRUE))
          output <- comb.intersect(var1, var2)
          output_name <- paste0(variable_combinations[1, i], "_", variable_combinations[2, i])
          assign(output_name, output, envir = .GlobalEnv)
          intersect_list[[output_name]] <- output
        }
        
        # If only 1 intersect (only 2 reps), then that intersect is the consensus
        if (length(intersect_list)==1) {
          write.table(intersect_list, paste0(out_dir, "/", common_part, ".bed"), sep='\t', row.names=F, col.names=F, quote=F)
        } else {
          # If >2 rep, then union is taken 
          # Apply union to the list of intersections
          consensus <- comb.union2(intersect_list)
          write.table(consensus, paste0(out_dir, "/", common_part, ".bed"), sep='\t', row.names=F, col.names=F, quote=F)
        }
      }
    }
  }
}

findConsensus(args[1], args[2])
