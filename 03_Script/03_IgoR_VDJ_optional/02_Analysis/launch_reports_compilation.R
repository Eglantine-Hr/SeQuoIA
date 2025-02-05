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
  make_option(c("-P", "--PGEN_TABLE"), type="character", default=NULL, 
              help="dataframe summarizing all sequences scores", metavar="character"),
              c("-V", "--VDJ_COMBIS_TABLE"), type="character", default=NULL, 
              help="dataframe summarizing VDJ segments names", metavar="character"),
              c("-S", "--SELECT_SCORE"), type="character", default=NULL, 
              help="dataframe summarizing all sequences scores", metavar="character"),
  make_option(c("-O", "--OUTPUT_REPORT"), type="character", default=NULL, 
              help="final report", metavar="character")#,

); 



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#INPUT_TABLE = paramsEnv[["INPUT_TABLE"]]
#OUTPUT_PATH =  paramsEnv[["OutputFile"]]  

ALL_SCORES_TABLE = opt$INPUT_TABLE
OUTPUT_REPORT = opt$OUTPUT_REPORT
PGEN_TABLE = opt$PGEN_TABLE
VDJ_COMBIS_TABLE = opt$ALL_CLONES_TABLE




### Render the report using previously built environment (use link to Rmd file)
rmarkdown::render( input = file.path(WORKING_DIR, "Analysis_All_Samples.Rmd"),
                   output_dir = dirname(OUTPUT_REPORT),
                   output_file  = I( "Analysis_All_Samples.html"),
                   envir = paramsEnv,
                   quiet = FALSE)

# Remove temporary 'Rmd' copy
#if(! file.remove( rmdCopyName)) warning( paste0( "Temporary file '", rmdCopyName, "' could not be removed..."));
