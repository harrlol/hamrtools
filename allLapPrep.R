library(data.table)
library(tidyr)
library(dplyr)

args=commandArgs(trailingOnly=TRUE)

lap.clean <- function(bed, type) {
  if (type=="gene") {
    out <- select(bed, seq=V7, pos=V8, mod=V10, gene=V4, score=V11, strand=V12)
  } else if (type=="UTR") {
    out <- select(bed, seq=V7, pos=V8, mod=V10, gene=V5, feature=V6, score=V11, strand=V12)
  } else if (type=="CDS") {
    out <- select(bed, seq=V6, pos=V7, mod=V9, gene=V5, score=V10, strand=V11)
  } else if (type=="mRNA") {
    out <- select(bed, seq=V7, pos=V8, mod=V10, gene=V4, score=V11, strand=V12)
  } else {stop("invalid type")}
  return(out)
}

isoUnsensitive <- function (modtbl) {
  modtbl%>%
    group_by(across(c(-gene)))%>%
    summarize(gene = unique(unlist(strsplit(as.character(gene),".", fixed = TRUE))[1]))
}

allLapPrep <- function(in_dir) {
  file_names <- list.files(path = in_dir)
  
  # Initialize long df
  longdf <- NULL
  for (file_name in file_names) {
    
    # First split the file 
    finfo <- unlist(strsplit(tools::file_path_sans_ext(file_name), split = "_", fixed=TRUE))
    
    # Obtain the type of overlap of current file (UTR, mRNA, CDS, gene)
    lap_type <- finfo[length(finfo)]
    
    # Obtain the sequencing technique (mRNA seq, GMUCT, NAD, etc.)
    seq_tech <- finfo[length(finfo) - 1]
    
    # Obtain the genotype of the sample group
    genotype <- paste(setdiff(finfo, c(lap_type, seq_tech)), collapse = "_")
    
    # Set file path for testing and reading
    fpath <- file.path(in_dir, file_name)
    
    # Some files can have no overlaps, consider only those that have
    if (file.size(fpath)!=0) {
      # Import the file as a variable and pipe through the cleaning steps 
      temp <- assign(file_name, fread(file.path(in_dir, file_name), stringsAsFactors = TRUE))
      temp_clean <- lap.clean(temp, lap_type)
      
      # UTR needs a bit more work because 5' and 3' isoform cleaning
      if (lap_type == "UTR") {
        # Add experimental information alongside hamr predictions and bind to long df
        to_add1 <- data.frame(isoUnsensitive(temp_clean))%>%
          filter(feature=="three_prime_UTR")%>%
          mutate(genotype=genotype) %>%
          mutate(seq_tech=seq_tech) %>%
          mutate(lap_type="3'UTR") %>%
          select(-feature)
        to_add2 <- data.frame(isoUnsensitive(temp_clean))%>%
          filter(feature=="five_prime_UTR")%>%
          mutate(genotype=genotype) %>%
          mutate(seq_tech=seq_tech) %>%
          mutate(lap_type="5'UTR") %>%
          select(-feature)
        longdf <- rbind(longdf, to_add1, to_add2)
      } else if (lap_type == "CDS") {
        to_add <- data.frame(isoUnsensitive(temp_clean))%>%
          mutate(genotype=genotype) %>%
          mutate(seq_tech=seq_tech) %>%
          mutate(lap_type=lap_type)
        longdf <- rbind(longdf, to_add)
      } else {
        to_add <- data.frame(temp_clean)%>%
          mutate(genotype=genotype) %>%
          mutate(seq_tech=seq_tech) %>%
          mutate(lap_type=lap_type)
        longdf <- rbind(longdf, to_add)
      }
    }
  }
  return(longdf)
}

project <- allLapPrep(args[1])

write.table(project, paste0(args[2], "/mod_long.csv"), sep='\t', row.names=F, col.names=T, quote=F)