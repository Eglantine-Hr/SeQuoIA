# ####################################################################
# This script launch the compilation of report
# ####################################################################

library( knitr)
library( rmarkdown)
library( funr)
library(optparse)



# Define working folder (contains R/Rmd file for current sample, parent contains global project files)
if( exists( "snakemake"))
{
  WORKING_DIR = snakemake@scriptdir
} else
{
  WORKING_DIR = dirname( sys.script())
}


### Load parameters
# Define an environment that will contain parameters
paramsEnv = new.env();

# Parameterize option to lauch the script for each sample passed as arg
option_list = list(
  make_option(c("-I", "--INPUT_TABLE"), type="character", default=NULL, 
              help="dataframe summarizing all sequences scores", metavar="character"),
  make_option(c("-V", "--VDJ_FOLDER"), type="character", default=NULL, 
              help="VDJ summary", metavar="character"),
  make_option(c("-M", "--METADATA_PATH"), type="character", default=NULL, 
              help="METADATA_PATH", metavar="character"),
  make_option(c("-O", "--OUTPUT_PATH"), type="character", default=NULL, 
              help="final report", metavar="character")#,

); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#INPUT_TABLE = paramsEnv[["INPUT_TABLE"]]
#OUTPUT_PATH =  paramsEnv[["OutputFile"]]  

ALL_SCORES_TABLE = opt$INPUT_TABLE
OUTPUT_PATH = opt$OUTPUT_PATH
VDJ_FOLDER = opt$VDJ_FOLDER
METADATA_PATH = opt$METADATA_PATH



#ALL_CLONES_TABLE  = file.path(  dirname(ALL_SCORES_TABLE), "Clonal_Composition_All.txt") 
ALL_CLONES_TABLE <- list.files(path = VDJ_FOLDER,  pattern = "Clonal_Composition_All.txt", full.names = T)[1]
print(VDJ_FOLDER)
print(ALL_CLONES_TABLE)

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }

  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

### Render the report using previously built environment (use link to Rmd file)
rmarkdown::render( input = file.path(WORKING_DIR, "Analysis_All_Samples.Rmd"),
                   output_dir = dirname(OUTPUT_PATH),
                   output_file  = I( "Analysis_All_Samples.html"),
                   envir = paramsEnv,
                   quiet = FALSE)

# Remove temporary 'Rmd' copy
#if(! file.remove( rmdCopyName)) warning( paste0( "Temporary file '", rmdCopyName, "' could not be removed..."));
