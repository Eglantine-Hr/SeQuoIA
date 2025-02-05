###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
###############################################################################


#### General

GLOBAL_DESCRIPTION = "CENTURY PROJECT - DLBCL SCRNASEQ ANALYSIS"

SCIENTIFIC_GROUP = "Collaboration"
SCIENTIFIC_PROJECT_NAME = "LYMPHOMA_ATLAS"
DEV_REP = "04_BIOINFO_DEV/baaklinis"
EXPERIMENT_PROJECT_NAME = "06_FL_Citeseq_setup_experiment"
SAMPLE_NAME = c("Control","Control_TestBarcode","DilutionFactor1","DilutionFactor4")



#### Input / Output

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_PROJECT = file.path( "/mnt/NAS7", 
                          SCIENTIFIC_GROUP,
                          SCIENTIFIC_PROJECT_NAME,
                          DEV_REP)

PATH_PROJECT_RAWDATA = file.path(PATH_PROJECT,EXPERIMENT_PROJECT_NAME)

PATH_PROJECT_OUTPUT = file.path( PATH_PROJECT, EXPERIMENT_PROJECT_NAME, "04_Output")

# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_",     
                            EXPERIMENT_PROJECT_NAME, "_")


#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



