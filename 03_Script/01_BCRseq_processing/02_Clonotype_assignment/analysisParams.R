######################################################
#####   This file contains analysis parameters  ######
#####     parameters that should not be changed ######
######################################################

# ANALYSIS PARAM ------------------------------------

ANALYSIS_STEP_NAME = "02_Clonotype_Assignment"
PATH_ANALYSIS_OUTPUT = opt$OutputDir
SAMPLE_NAME = opt$sample
BCR_INPUT_FILE = opt$file
SEURAT_INPUT_FILE = opt$input
AIRR_OUTPUT_FILE = opt$AirrFormatOutput
SEURAT_OUTPUT_FILE = opt$SeuratCompatibleOutput
METADATA = opt$Metadata
METADATA_FOLDER = opt$METADATA_FOLDER

DONOR_COLUMN = opt$DONOR_COLUMN
METADATA_FOLDER = opt$METADATA_FOLDER

# CLONOTYPE SIMILARITY COMPUTATION
DISTANCE_METHOD = "normlevenshtein"
DISTANCE_MISMATCH = 1
DISTANCE_GAP = 1
DISTANCE_CASE_SENSITIVITY = TRUE
HC_METHOD = "complete"
LC_METHOD = "complete"




MODE = "VJ_CDR3"

#NP 
HC_CUTREE = 0.15
LC_CUTREE = 0.15

HC_CUTREE = 0.2
LC_CUTREE = 0.2

