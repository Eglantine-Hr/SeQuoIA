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

# Parameterize option to lauch the script for each sample passed as arg CELL_METADATA_FOLDER
option_list = list(
  make_option(c("-i", "--InputPath"), type="character", default=NULL, 
              help="input objects repertory path", metavar="character"),           
  make_option(c("-p", "--NPROC"), type="numeric", default=NULL, 
              help="max number processors used", metavar="numeric"),
  make_option(c("-c", "--CPU_USAGE"), type="numeric", default=NULL, 
              help="max number processors used", metavar="numeric"),
  make_option(c("-M", "--CELL_METADATA_FOLDER"), type="character", default="test", 
              help="pasth to scRNAseq metadata", metavar="character"), 
  make_option(c("-A", "--airrfile"), type="character", default=NULL, 
              help="Airr file with pairs of chains cells ids", metavar="character")#,

); 

#  make_option(c("-O", "--OutputDir"), type="character", default=NULL, 
#              help="Output repertory name", metavar="character"
#  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#mcmapply parallelization 
print(opt$NPROC)
print(parallelly::availableCores())

max_cores = min( parallel::detectCores() , opt$NPROC  ) 

#available_cores = as.numeric(parallelly::availableCores())

available_cores = floor(max_cores - (as.numeric(opt$CPU_USAGE)/100)*parallel::detectCores()) #% cpu usage retrieved from bash top command 

print(paste0("available cores :    ", available_cores))
n_cores = floor(available_cores/(1+ max_cores -available_cores))
n_cores = floor(max(1, n_cores, na.rm = T ))
#n_cores = floor(max(1, n_cores/5, na.rm = T ))
#n_cores = 1
print(paste0("final cores to use  :    ", n_cores ))

# Load file defining analysis parameters
analysisParamsFilePath = file.path( WORKING_DIR, "analysisParams.R");

if( file.exists(analysisParamsFilePath) ){ source( analysisParamsFilePath, local = paramsEnv)
} else { warning("The file 'analysisParamsFilePath.R' containing analysis-specific parameters is missing.") }

### Prevent overwriting
# Boolean to activate overwriting security (change at your own risk !)
preventOverwrite = FALSE

# Create output file name
reportOutputFilename = paste0( paramsEnv[["ANALYSIS_STEP_NAME"]], ".html");


print(file.path( paramsEnv[["PATH_ANALYSIS_OUTPUT"]], reportOutputFilename) )
alreadyExists = file.exists( file.path( paramsEnv[["PATH_ANALYSIS_OUTPUT"]], reportOutputFilename) )




if( alreadyExists  && preventOverwrite){ stop( paste( c( "Report file already exists:", reportOutputFilename))) }

GCTREE_FOLDER = paramsEnv[["INPUT_OBJECTS_PATH"]]
CONTIG_AIRR_FILE =  paramsEnv[["AIRR_FILE"]]  
print(GCTREE_FOLDER)
print(CONTIG_AIRR_FILE)





### Prevent result mixing when parallel rendering
# Create a copy of original rmd file with unique name.
# Prevents 'render' to overwrite 'md' file created from name of 'Rmd' file,
# and scrambled results when parallel-rendering 'Rmd' files with identical name.
# Specifying 'tmp' folders for 'intermediates_dir' and 'knit_root_dir' prevent
# usage of WORKING_DIR in reports.
# Final output 'html' name provided to render anyway (not based on 'Rmd' name).

# Create a unique temporary file name
rmdCopyName = tempfile( pattern = "tempRmdCopy_", tmpdir = WORKING_DIR, fileext = ".Rmd");
# Copy 'Rmd' file to it
stopifnot( all( file.copy( from = file.path( WORKING_DIR, paste0(  paramsEnv[["ANALYSIS_STEP_NAME"]], ".Rmd")),
                           to   = rmdCopyName)));


### Render the report using previously built environment (use link to Rmd file)
print(GCTREE_FOLDER)
rmarkdown::render( input = rmdCopyName,
                   output_dir = paramsEnv[["PATH_ANALYSIS_OUTPUT"]],
                   output_file  = I( reportOutputFilename),
                   envir = paramsEnv,
                   quiet = FALSE)

# Remove temporary 'Rmd' copy
if(! file.remove( rmdCopyName)) warning( paste0( "Temporary file '", rmdCopyName, "' could not be removed..."));
