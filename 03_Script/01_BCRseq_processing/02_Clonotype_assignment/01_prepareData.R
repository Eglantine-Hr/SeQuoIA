#####################################################################
#######   The aim of this script is to load single cell BCR    ######
#######                   dataset                              ######
#####################################################################


# SET UP CLUSTER PARALLELISATION ------------------------------------
## @knitr clusterSetup


n_cores<-detectCores()-1
#n_cores = 20
cl<-makeCluster(n_cores)

registerDoParallel(cl)

# Add a combined function for parallelization
comb <- function(...) {
  mapply(rbind, ..., SIMPLIFY=FALSE)
}


# PREPARE THE DATA --------------------------------------------------------

## @knitr prepareData 

# File and params infos
print(current_donor) 
print(dim(bcr.df.all)) 
print(dim(bcr.df)) 
print(head(current_donor_barcodes))
print(head(bcr.df.all$barcode))

# Loading BCR dataset
#bcr.df <- read_tsv(file = BCR_INPUT_FILE)
bcr.df$Absolute_Sequence = gsub("[^0-9A-Za-z///']", "", bcr.df$sequence_alignment  )
bcr.df$Original_Sequence = bcr.df$sequence  
bcr.df$IMGT_Sequence =  bcr.df$sequence_alignment  
if("contig_id" %notin% names(bcr.df)){bcr.df$contig_id = bcr.df$sequence_id }

# Loading contig ID data 
Contig_IDs = seqinr::read.fasta(file = SEURAT_INPUT_FILE, as.string = TRUE, set.attributes = FALSE, forceDNAtolower = FALSE) 
Contig_IDs = as.data.frame( t( data.frame( Contig_IDs, check.names = F) ) ) %>% dplyr::rename(Original_Sequence = V1) # %>% rownames_to_column("contig_id")
Contig_IDs$sequence_id  =as.character(rownames(Contig_IDs)) 
print( colnames(Contig_IDs) ) 


# Linking contig ids to cell barcodes
metadata_files <- list.files(path= dirname(SEURAT_INPUT_FILE), pattern= paste0("*", METADATA, "*"), full.names=TRUE, recursive=FALSE)
print(paste0("metadata files", metadata_files) )

if ( length(metadata_files) > 0 ){
	metadataDF = read.table(file = metadata_files[1], header = T, sep = ",")
	#metadataDF = read.csv(file = metadata_files[1], header = T)
	print(colnames(metadataDF))
	unique_colnames = colnames(metadataDF)[colnames(metadataDF) %notin% colnames(bcr.df) ]
	print(unique_colnames)
	if("barcode" %in% names(metadataDF) ){ metadataDF$cell_id = metadataDF$barcode }
	print("BARCODES:")
	print(head(metadataDF$cell_id))
	metadataDF$sequence_id = as.character(metadataDF$contig_id) # to be changed according to metadata ????

	
	metadataDF = metadataDF %>% dplyr::select( c(unique_colnames,"sequence_id", "cell_id" ) ) #%>% dplyr::filter(contig_id %in% Contig_IDs$contig_id ) 
	bcr.df = dplyr::left_join(bcr.df,  metadataDF  ) #[, c(unique_colnames,"contig_id" ) ] #, by = "contig_id"
	
	Contig_IDs = dplyr::left_join(Contig_IDs,  metadataDF , by = "sequence_id" ) 
	print( unique(Contig_IDs$Chain)) 
}else{Contig_IDs$cell_id = Contig_IDs$sequence_id }

# Ascribe contig ID to each row 
#bcr.df = bcr.df[which(bcr.df$Original_Sequence %in% Contig_IDs$Original_Sequence),]
Contig_IDs = Contig_IDs %>% dplyr::select(-Original_Sequence)
#bcr.df = dplyr::left_join(Contig_IDs,  unique_sequences , by = "Original_Sequence" ) 
#bcr.df = dplyr::left_join(bcr.df,  Contig_IDs , by = "sequence_id" )


print(paste0("BCR_INPUT_FILE : ",  BCR_INPUT_FILE))
print(paste0("SEURAT_INPUT_FILE : ", SEURAT_INPUT_FILE))
print(paste0("METADATA_FOLDER : ", METADATA_FOLDER))
print(paste0("metadata_files : ", metadata_files))


# Split table in heavy and light chains
if("Chain" %in% names(bcr.df) ){ bcr.df$chains = recode( bcr.df$Chain, "IGH" = "Heavy" , "IGK" = "Light",  "IGL" = "Light") 
}else{bcr.df$chains = recode( bcr.df$locus, "IGH" = "Heavy" , "IGK" = "Light",  "IGL" = "Light")  } 

#print(names(bcr.df)  )
print(head(bcr.df$locus))
print(head(metadataDF$Chain))

bcr.df.split <- bcr.df %>% 
  #dplyr::mutate(chains = dplyr::recode(locus, "IGH" = "Heavy" , "IGK" = "Light",  "IGL" = "Light")) %>% # Rename IGH/IGL/IGK in chains column as heavy or light chain #locus
  #dplyr::mutate(chains  = ifelse(locus == "IGH" , "Heavy","Light" )) %>% 
  #				     #locus %in% c("IGK", "IGL" ) ~ "Light" )  ) %>%  
  group_split(chains, .keep = FALSE) # create 1 list of 2 df : 1 for heavy chain and one for light chain

# Add heavy prefix to colnames and extract vj genes without allele info 
Heavy_Chain_df<-bcr.df.split[[1]]  %>%  group_by(cell_id) %>% arrange(desc("umis")) %>% slice(1) %>% 
  setNames(paste0('heavy_', names(.))) %>% 
  dplyr::rename(cell_id = heavy_cell_id) %>% 
  mutate(heavy_vj_call_short = paste0(gsub("\\*.*", "", heavy_v_call),"_",gsub("\\*.*", "", heavy_j_call)))
print(paste0("dim heavy chain data:  ", dim(Heavy_Chain_df)))

# Add light prefix to colnames and extract vj genes without allele info 
Light_Chain_df<-bcr.df.split[[2]] %>%   group_by(cell_id) %>% arrange(desc("umis")) %>% slice(1) %>% #dplyr::filter(umis > 6 & reads > 50 ) %>% 
  setNames(paste0('light_', names(.))) %>% 
  dplyr::rename(cell_id = light_cell_id)%>% 
  mutate(light_vj_call_short = paste0(gsub("\\*.*", "", light_v_call),"_",gsub("\\*.*", "", light_j_call)))
print(paste0("dim light chain data:  ", dim(Light_Chain_df)))


if( exists( "Light_Chain_df") == F ){
Light_Chain_df<-bcr.df.split[[1]] %>% slice(1) %>%   group_by(cell_id) %>% arrange(desc("umis")) %>% slice(1) %>% #dplyr::filter(umis > 6 & reads > 50 ) %>% 
  setNames(paste0('light_', names(.))) %>% 
  dplyr::rename(cell_id = light_cell_id)%>% 
  dplyr::mutate(light_vj_call_short = paste0(gsub("\\*.*", "", light_v_call),"_",gsub("\\*.*", "", light_j_call)),
  cell_id = 'virtual' )
print(paste0("dim light chain data:  ", dim(Light_Chain_df)))
}


#Memory optimization 
rm(bcr.df.split, Contig_IDs, metadataDF)
gc()

