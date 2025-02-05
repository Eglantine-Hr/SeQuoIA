# ####################################################################
# This script launch the compilation of report
# ####################################################################

library( knitr)
library( rmarkdown)
library( funr)
library(optparse)
library(readr)
`%notin%` <- Negate(`%in%`)


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
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="bcr file in AIRR format", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="B seurat object", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample name", metavar="character"),
  make_option(c("-S", "--SeuratCompatibleOutput"), type="character", default=NULL, 
              help="Seurat Compatible output file name", metavar="character"),
  make_option(c("-a", "--AirrFormatOutput"), type="character", default=NULL, 
              help="Airr Format Output file name", metavar="character"),
  make_option(c("-p", "--PatternFile"), type="character", default=NULL, 
              help="File pattern", metavar="character"),
  make_option(c("-m", "--Metadata"), type="character", default=NULL, 
              help="Single cell metadata", metavar="character"),
  make_option(c("-L", "--DONOR_COLUMN"), type="character", default=NULL, 
              help="Donor Column name in metadata dataframe", metavar="character"), 
  make_option(c("-e", "--METADATA_FOLDER"), type="character", default=NULL, 
              help="Single cell metadata", metavar="character"),           
  make_option(c("-O", "--OutputDir"), type="character", default=NULL, 
              help="Output dir name", metavar="character"
  )
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt$sample)
print(opt$PatternFile)
print(WORKING_DIR)

print(opt$DONOR_COLUMN) 
print(opt$METADATA_FOLDER)
#DONOR_COLUMN = opt$DONOR_COLUMN
METADATA_FOLDER = opt$METADATA_FOLDER




####################################################
# Data frame processing 
####################################################

detect_separator = function(filename){
	asAText = readLines(filename, n = 1)
	
	num_semicolons	<- count.fields(textConnection(asAText), sep = ";")
	num_spaces	<- count.fields(textConnection(asAText), sep = " ")
	num_commas	<- count.fields(textConnection(asAText), sep = ",")
	num_tabs	<- count.fields(textConnection(asAText), sep = "\t")
	
	all_seps = c( ";", " ", ",", "\t" )
	all_counts = c(num_semicolons, num_spaces, num_commas, num_tabs)
	final_sep = all_seps[ which.max(all_counts)  ]
	
	return( final_sep)	
}

#METADATA_FOLDER = "/mnt/DOSI/PMLAB/BIOINFO/FL_modeling/01_FL_vs_Physio_data_analysis/00_Rawdata/02_scRNAseq/LN_BLOOD_BM_kim_nature_2022"


print("Parameters fetched") 
#DONOR_COLUMN = "Donor" # temporary 

# MANAGE NULL ARG 

# if (is.null(opt$file) | is.null(opt$out) | is.null(opt$seuratCompatibleOutput) ){
#   print_help(opt_parser)
#   stop("4 arguments must be supplied (input file name, sample name, airrformat output file name, seurat compatible output file name).n", call.=FALSE)
# }

# Load file defining global parameters
globalParamsFilePath = file.path( WORKING_DIR, "../globalParams.R");
if(file.exists(globalParamsFilePath)) 
{
  source( globalParamsFilePath, local = paramsEnv);
  source( globalParamsFilePath)
} else 
{
  warning("The file 'globalParamsFilePath.R' containing global parameters is missing.");
}


# Load file defining analysis parameters
analysisParamsFilePath = file.path( WORKING_DIR, "analysisParams.R");
if(file.exists(analysisParamsFilePath)) 
{
  source( analysisParamsFilePath, local = paramsEnv);
  source( analysisParamsFilePath)
  source(file.path( WORKING_DIR, "00_generalDeps.R") )
} else 
{
  warning("The file 'analysisParamsFilePath.R' containing analysis-specific parameters is missing.");
}

### Prevent overwriting
# Boolean to activate overwriting security (change at your own risk !)
preventOverwrite = TRUE;
print(paramsEnv[["ANALYSIS_STEP_NAME"]])
# Create output file name
reportOutputFilename = paste0( paramsEnv[["SCIENTIFIC_PROJECT_NAME"]], "_",
                               paramsEnv[["ANALYSIS_STEP_NAME"]], ".html");
reportOutputFilename = paste0( paramsEnv[["ANALYSIS_STEP_NAME"]], "_",
                               opt$PatternFile, ".html");

alreadyExists = file.exists( file.path( paramsEnv[["PATH_ANALYSIS_OUTPUT"]], reportOutputFilename));
if( alreadyExists && preventOverwrite) stop( paste( c( "Report file already exists:", reportOutputFilename)));



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
stopifnot( all( file.copy( from = file.path( WORKING_DIR, paste0( "Report_", paramsEnv[["ANALYSIS_STEP_NAME"]], ".Rmd")), to   = rmdCopyName)));


# Split by donors if relevant to speed up clonotype assignation ------------


scRNAseq_metadata_path = list.files(path = METADATA_FOLDER,  pattern = "etadata|ETADATA|etaData", full.names =T)[1]
print(METADATA_FOLDER)
print(paste0( "metadata FOUND:    ", scRNAseq_metadata_path  ) )
metadata_df = read.table(file = scRNAseq_metadata_path, header = T, sep =  detect_separator(scRNAseq_metadata_path) ) #detect_separator(scRNAseq_metadata_path)
#metadata_df = read.table(file = "/mnt/DOSI/PMLAB/BIOINFO/FL_modeling/01_FL_vs_Physio_data_analysis/00_Rawdata/02_scRNAseq/LN_BLOOD_BM_kim_nature_2022/metadata_kim_et_al_2022_EH.txt", header = T, sep =  " " )
print(paste0( "BCR INPUT FOUND:    "  ) )
print(BCR_INPUT_FILE)
bcr.df.all <- read_tsv(file = BCR_INPUT_FILE)
#bcr.df.all = bcr.df.all %>% group_by(sequence) %>% slice(1) %>% ungroup() # TO REMOVE AFTER, FOR YEAP ONLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
k = 1

CLONE_Specification <-  str_replace(DONOR_COLUMN, "None", "All") 
print("START LOOP")
if(DONOR_COLUMN == "None" | length(na.omit(scRNAseq_metadata_path)) < 1 ){  
	bcr.df = bcr.df.all
	other_donor_barcodes = metadata_df[ which( metadata_df[ DONOR_COLUMN ]  !=  current_donor  ) , ]$barcode 
	rmarkdown::render( input = rmdCopyName,
                   output_dir = paramsEnv[["PATH_ANALYSIS_OUTPUT"]],
                   output_file  = I( reportOutputFilename  ),
                   #envir = paramsEnv,
                   quiet = FALSE)
       
       
       bcr.df.clono_list = bcr.df.clono
       bcr.df_list = bcr.df
       
       # Not mandatory, if whole clono approach 
       #if(    any( grepl("[Dd]onor", names(metadata_df  ) )     ){
       #	metadata_df_barcodes = metadata_df[, which( names(metadata_df  ) %in% c("cell_id", "barcode" )   ) ][,1]  
       #	metadata_df_donors = metadata_df[, which( grepl("[Dd]onor", names(metadata_df  )  ) ][,1]  
       #	metadata_df = cbind( metadata_df_barcodes,metadata_df_donors ) 
       #	
       #	if()
       #	bcr.df.all$cell_id
       #	colnames(metadata_df) = c("barcode", "donor") 
       #}
       
        
}else{  
	print(paste0("DONOR_COLUMN   =   ", DONOR_COLUMN))
	if(DONOR_COLUMN %in% names(metadata_df) == F ){ DONOR_COLUMN = names(metadata_df)[ which(  names(metadata_df) %in% c("donor", "Donor", "patient", "Patient", "mouse", "MULTI_ID" )) ][1]    }  
	
	#grepl( "donor|patient|mouse|MULTI_ID",  names(metadata_df) )   ) ][1]
	all_donors =  unique( as.vector( metadata_df[, which(  names(metadata_df) ==  DONOR_COLUMN   )   ] ) )
	print(all_donors)
	
	### Initialization 
	bcr.df.clono_list = list() 
	bcr.df_list = list() 
	
	print(names(bcr.df.all))
	print(BCR_INPUT_FILE)

	print("FOR LOOP DONOR") 
	#print(names(metadata_df))
	write.table(bcr.df.all, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/bcr.df.all.txt" )
	for(  current_donor in all_donors     ){
		print(current_donor)
	
		### Substing inputs -----
		if("barcode" %in% names(metadata_df) == F & "cell_id" %in% names(metadata_df) == F ){  
					 metadata_df$barcode = rownames(metadata_df) } # barcodes should be in specific columns as rownames of metadata  
		 
		
		##other_donor_barcodes = metadata_df[ which( metadata_df[ DONOR_COLUMN ]  !=  current_donor  ) , ]$barcode # non assigned BCR seq + current donor BCR are taken
		other_donor_barcodes = unique(metadata_df[ which( metadata_df[ DONOR_COLUMN ]  !=  current_donor  ) , ][,which(names(metadata_df) %in% c("barcode", "cell_id") )])
		#print(current_donor_barcodes[1:20])
		print(head( bcr.df.all$barcode))
		if("barcode" %in% names(bcr.df.all)  ){
			 bcr.df = bcr.df.all[ which(    bcr.df.all$barcode %notin%  other_donor_barcodes )  , ]
			 bcr.df.all$cell_id = bcr.df.all$barcode
		}else{ bcr.df = bcr.df.all[ which(    bcr.df.all$cell_id %notin%  other_donor_barcodes )  , ] } 
		print(dim(bcr.df.all)) 
		print(head(rownames(bcr.df.all)))
		print(dim(bcr.df)) 
		print(DONOR_COLUMN)
		print(nrow(bcr.df) )
		if(nrow(bcr.df) == 0){next}
		
		### Render the report using previously built environment (use link to Rmd file)------------------------------------------
		rmarkdown::render( input = rmdCopyName,
                   output_dir = paramsEnv[["PATH_ANALYSIS_OUTPUT"]],
                   output_file  = I(    paste0( current_donor , "_", reportOutputFilename)  ),
                   #envir = paramsEnv,
                   quiet = FALSE)
               
               
               ## SAVE THE OUTPUT FILE -------------------
               bcr.df.clono$cloneID = paste0( bcr.df.clono$cloneID, ".", k, "_", current_donor ) #, "_", current_donor
               write.table(bcr.df.clono, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/bcr.df.clono.txt" )
                bcr.df.clono = as.data.frame( bcr.df.clono )
               bcr.df.clono = bcr.df.clono[,which(  grepl(".[xy]" , names(bcr.df.clono)) == F   )  ] 
		bcr.df.clono_list[[k]] = as.data.frame( bcr.df.clono )
		bcr.df_list[[k]] = as.data.frame( bcr.df )
		k = k +1 
		print("AFTER MERGE")
		print(dim(bcr.df.clono))
		print(dim( bcr.df))
		
		#write.table(bcr.df.clono, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/bcr.df.clono.txt" ) 
		#write.table(bcr.df, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/bcr.df.txt" )

	
	}
	print("END FOR LOO --------------")
	print(dim( bcr.df))
	print(SEURAT_OUTPUT_FILE)
	
	saveRDS(bcr.df.clono_list, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/bcr.df.clono_list.rds" )
	saveRDS(bcr.df_list, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/bcr.df_list.rds" )
	print( dim(do.call("rbind", bcr.df.clono_list)) )
	print(AIRR_OUTPUT_FILE)
	print( dim(do.call("rbind", bcr.df_list )) )
	#bcr.df.clono_list = as.data.frame( do.call("rbind", bcr.df.clono_list) ) %>% group_by(cell_id) %>% slice(1)
	  
	
	#contig without donor are ascribed to largest clone in terms of real cells 
	cell_ids_metadata = unique(as.vector(metadata_df[,which(names(metadata_df) %in% c("barcode", "cell_id") )]    ))
	bcr.df.clono_list = as.data.frame( do.call("rbind", bcr.df.clono_list) ) %>%
				dplyr::mutate( isInMetadata = ifelse(cell_id %in% cell_ids_metadata, T, F )  ) %>%
				##dplyr::mutate( isInMetadata = ifelse(cell_id %in% metadata_df$barcode , T, F )  ) %>%
				group_by(cloneID) %>% dplyr::mutate(Sublclone_Size= sum(isInMetadata == T  ) ) %>% ungroup() #%>% 
				##group_by( heavy_sequence_id  ) %>% dplyr::mutate(n_clones_H = n_distinct(cloneID ) ) %>% #%>% dplyr::filter(n_clones_H < 2 ) %>% ungroup() %>%
				##group_by( light_sequence_id  ) %>% dplyr::mutate(n_clones_L = n_distinct(cloneID ) ) #%>% dplyr::filter(n_clones_L < 2 ) %>% ungroup()
				#group_by(cloneID) %>% dplyr::mutate(Sublclone_Size= sum(isInMetadata == T  ), 
				#				     Sublclone_Size = ifelse(grepl( "Not_CLONAL", cloneID ), 0,   Sublclone_Size     ) ) %>%   
				#ungroup() %>% group_by(cell_id) %>% dplyr::arrange( desc(Sublclone_Size) ) %>% slice(1)
	
	bcr.df_list = as.data.frame( do.call("rbind", bcr.df_list )  ) #%>% group_by(sequence_id)  %>% dplyr::arrange(cloneID) %>% slice(1)
	
	
	
	
	print(names(bcr.df.clono_list))
	#bcr.df_list_clonal  = bcr.df_list %>% dplyr::filter(  grepl("NOT_CLONAL" , cloneID ) ==F ) 
	#bcr.df_list = bcr.df_list %>% dplyr::filter(  grepl("NOT_CLONAL" , cloneID )  &    sequence_id %notin%  bcr.df_list_clonal$sequence_id  )
	#bcr.df_list = rbind(bcr.df_list, bcr.df_list_clonal ) 
	
	write.table(bcr.df.clono_list, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/bcr.df.clono_list.txt" )
	bcr.df.clono_list_clonal  = bcr.df.clono_list %>% dplyr::filter(  grepl("NOT_CLONAL" , cloneID ) ==F ) %>%
				    #group_by(cloneID) %>% dplyr::mutate(Sublclone_Size= sum(isInMetadata == T  ) ) %>% ungroup() %>% 
				    group_by(cell_id) %>% 
				    dplyr::mutate(n_clones = n_distinct(cloneID[Sublclone_Size>0])  )  %>% 
				    dplyr::arrange( desc(Sublclone_Size) ) %>% slice(1) %>% ungroup() %>%
				    dplyr::rowwise() %>% dplyr::filter(  ( n_clones > 1 & Sublclone_Size > 15 ) ==F  ) %>% dplyr::select( -  c(n_clones ) ) 
				    #dplyr::filter(n_clones_H < 2 & n_clones_L < 2  )
	bcr.df.clono_list = bcr.df.clono_list %>% dplyr::filter(  grepl("NOT_CLONAL" , cloneID )  &   
					 heavy_sequence_id %notin%  bcr.df.clono_list_clonal$heavy_sequence_id & light_sequence_id %notin%  bcr.df.clono_list_clonal$light_sequence_id )
	bcr.df.clono_list = rbind(bcr.df.clono_list, bcr.df.clono_list_clonal ) 
	
	print("END LOOP") 
} 




write_tsv(  bcr.df.clono_list, file = SEURAT_OUTPUT_FILE) # will serve in Tree selection (fitlered contig)
write_tsv(   bcr.df_list    , file = AIRR_OUTPUT_FILE)    





# Remove temporary 'Rmd' copy
if(! file.remove( rmdCopyName)) warning( paste0( "Temporary file '", rmdCopyName, "' could not be removed..."));
