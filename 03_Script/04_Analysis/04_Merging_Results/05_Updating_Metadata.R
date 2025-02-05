#########################################################################
#######      The aim of this script is to update         	     ######
#######   			10x scRNAseq metadata		     ######
#########################################################################


# Udating metadata with selection scores ---------------------------------------------------------
## @knitr updateMetadata


`%notin%` <- Negate(`%in%`)

metadata_df = list.files(path = METADATA_PATH,  pattern = "etadata|ETADATA|etaData", full.names =T)[1]

#Selecting columns that are not yet in metadata + contig id for merging 
Columns_to_Add = All_Scores_Df_Repeated_Seqs[which(names(All_Scores_Df_Repeated_Seqs) %notin% names(metadata_df )[which(names(metadata_df ) != "barcode")]  , ]   
#%>% dplyr::select("contig_id", "Chain", "MutSinceParent", "Global_Select_Score","Cumulated_Select_Score", "Clonotype")

#Pivot from long to wide (chain prefix) 
Columns_to_Add = Columns_to_Add %>%
			pivot_wider(names_from = Chain, values_from = c("SeqName", "MutSinceParent", "Global_Select_Score","Cumulated_Select_Score", "Clonotype"), names_sep = '.')


#updating metadata and saving to output folder 
if("barcode" %notin% names(barcode)  ){metadata_df$barcode = rownames(metadata_df)}
metadata_df = merge(metadata_df, Columns_to_Add, by = "barcode", all = FALSE, all.x = TRUE, all.y = FALSE)  
write.table(metadata_df, file = file.path(dirname(OUTPUT_PATH), "scRNAseq_metadata_updated.txt", header = T) 
