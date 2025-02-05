# LOAD LIBRARIES  --------------------------------------------------
library(dplyr)
library(magrittr)
library(tibble)
library(seqinr)
library(optparse)
library(purrr)
library(stringr)
library(readr)
library(stringdist)
library(DECIPHER)
library(stringi)
library(tidyr)
library(msa)



# FUNCTIONS  --------------------------------------------------
# Tests whether a sequence is productive (no frameshift and stop codons) 

`%notin%` <- Negate(`%in%`)

STOP_CODONS = c("TAA", "TAG", "TGA")



is_productive <- function(seqlist, stopcodons = c( "TAA", "TAG", "TGA")  ){
  stopcodons <- paste(stopcodons, collapse = "|")
  frameshift <- nchar(seqlist)%%3 # UCA sont des multiples de 3. A vérifier avec echidna 
  
  #sequence_codons <- stri_sub(seqlist, seq(1, stri_length(seqlist),by=3), length=3)
  sequence_codons <- gsub("(.{3})", "\\1 ", seqlist)
  
  #productive <- ifelse(frameshift != 0 |  grepl(pattern = stopcodons, x = seqlist) , FALSE, TRUE)
  productive <- ifelse(  grepl(pattern = stopcodons, x = sequence_codons) , FALSE, TRUE)
  
  return (productive)
}

Mode <- function(x, na.rm = FALSE) {
  # it takes two areguments x, and na.rm (https://www.educative.io/answers/what-is-the-mode-method-in-r)
  if(na.rm){ #if na.rm is false it means no need to remove NA values
    x = x[!is.na(x)]
  }
  valx <- unique(x)
  return(valx[which.max(tabulate(match(x, valx)))])
}
