
# LOADING LIBRARIES -------------------------------------------------------

library(optparse)
library(alakazam)
library(shazam)
library(dplyr)
library(ggplot2)
library(readr)

# PARAMETERIZE OPTIONS TO LAUNCH THE SCRIPT FOR EACH SAMPLE PASSED AS ARG -

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="bcr file in AIRR format with germline_alignment_d_mask info", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# MANAGE NULL ARG ---------------------------------------------------------

if (is.null(opt$file) | is.null(opt$out) ){
  print_help(opt_parser)
  stop("2 arguments must be supplied (input and output file).n", call.=FALSE)
}


## LOAD DATA 

bcr.df <- read_tsv(opt$file)


## MUTATION ANALYSIS 

# Calculate the counts of mutations per region (i.e cdr1/2/3 and fwr1/2/3/4)
bcr.df <- observedMutations(bcr.df, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment", 
                            regionDefinition=IMGT_VDJ_BY_REGIONS,
                            frequency=FALSE, 
                            nproc=1)

# Calculate the counts of mutations over the entire sequence
bcr.df$tot_mut_count <- rowSums(bcr.df[,grep("mu_count", names(bcr.df))])

# Calculate the frequencies of mutations per region (i.e cdr1/2/3 and fwr1/2/3/4)
bcr.df <- observedMutations(bcr.df, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment", 
                            regionDefinition=IMGT_VDJ_BY_REGIONS,
                            frequency=TRUE,
                            combine=FALSE,
                            nproc=1)

# Calculate the frequencies of mutations over the entire sequence
bcr.df <- observedMutations(bcr.df, sequenceColumn="sequence_alignment",
                            germlineColumn="germline_alignment", 
                            regionDefinition=IMGT_VDJ_BY_REGIONS,
                            frequency=TRUE,
                            combine=TRUE,
                            nproc=1)

# SAVE THE OUTPUT FILE 

write_tsv(bcr.df, file = opt$out)
