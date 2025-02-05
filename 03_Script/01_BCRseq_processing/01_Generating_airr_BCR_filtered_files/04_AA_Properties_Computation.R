
## LOADING LIBRARIES 

library(optparse)
library(alakazam)
library(dplyr)
library(ggplot2)
library(readr)

## PARAMETERIZE OPTIONS TO LAUNCH THE SCRIPT FOR EACH SAMPLE PASSED AS ARG -

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="bcr file in AIRR format with germline_alignment_d_mask info", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## MANAGE NULL ARG ---------------------------------------------------------

if (is.null(opt$file) | is.null(opt$out) ){
  print_help(opt_parser)
  stop("2 arguments must be supplied (input and output file).n", call.=FALSE)
}


## LOAD DATA 

bcr.df <- read_tsv(opt$file)


## AA PROPERTIES COMPUTATION 

# Calculate the properties of amino acid sequences
bcr.df <- aminoAcidProperties(bcr.df, seq="junction", nt=TRUE, trim=TRUE, 
                                label="cdr3")

## SAVE THE OUTPUT FILE 

write_tsv(bcr.df, file = opt$out)