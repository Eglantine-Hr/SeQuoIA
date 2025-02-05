######################################################
#####   This file contains analysis parameters  ######
#####     parameters that should not be changed ######
######################################################

# ANALYSIS PARAM ------------------------------------

ANALYSIS_STEP_NAME = "03_Tree_Selection" #paramsEnv[["ANALYSIS_STEP_NAME"]]
INPUT_OBJECTS_PATH = opt$Input
PATH_ANALYSIS_OUTPUT = opt$Input
AIRR_FILE = opt$airrfile
NPROC = opt$NPROC
#AVAILABLE_CORES = opt$AVAILABLE_CORES

CELL_METADATA_FOLDER = opt$CELL_METADATA_FOLDER

print(CELL_METADATA_FOLDER)
print(list.files(path = CELL_METADATA_FOLDER,  pattern = "etadata|ETADATA|etaData", full.names =T)[1])
