## @knitr LocalSelection_Groups

################################################
# Pheno based analysis (from scRNAseq data)  

################################################

scRNAseq_metadata_path = NA
ggplotly_list = list()
j = 0

# Checking if current dataset features scRNAseq metadata (all samples)
Metadata_Dir  =  file.path( sub("/03_Script.*", "", WORKING_DIR), "00_Rawdata/02_scRNAseq", DATASET)            
if(dir.exists( Metadata_Dir  ) ){
  
  
  scRNAseq_metadata_path = list.files(path = Metadata_Dir, 
                                 pattern = "etadata|ETADATA|etaData", full.names =T)[1] 
  print(Metadata_Dir)
  print("METADATA FILE FOUND: ")
  print(scRNAseq_metadata_path)
  if(length(na.omit(scRNAseq_metadata_path)) > 0 ){
    
    scRNAseq_metadata = read.table(file = scRNAseq_metadata_path, header = T, sep = detect_separator(scRNAseq_metadata_path) )
    

    if("barcode" %in% names(scRNAseq_metadata) == F & "cell_id" %in% names(scRNAseq_metadata) == F ){scRNAseq_metadata = scRNAseq_metadata %>% rownames_to_column(var = "barcode") }

    contig_annot = read.csv( list.files(file.path( dirname(dirname(Metadata_Dir)), "01_VDJ_BCR_contigs", DATASET, SAMPLE), pattern = "*annotations.csv" , full.names = TRUE)[1] )

    
    # Check if specified phenotypes are in metadata 
    print(paste0( "the following phenotypes were not found in metadata: ", PHENOS_OF_INTEREST[-which(PHENOS_OF_INTEREST %in% names(scRNAseq_metadata) )]  ))
    PHENOS_OF_INTEREST = PHENOS_OF_INTEREST[which(PHENOS_OF_INTEREST %in% names(scRNAseq_metadata) )]
    
    #getting rows corresponding to patient if merged object
    scRNAseq_metadata_sample= scRNAseq_metadata %>%  dplyr::filter(if_any( matches("ident$|ample$|atient$|maxID$|MULTI_ID$"), ~ grepl(SAMPLE, .)  )) 
    if(nrow(scRNAseq_metadata_sample)==0){scRNAseq_metadata_sample = scRNAseq_metadata %>% filter_all(any_vars(grepl(SAMPLE, .)))} #misleading HTO_secondID => priority to cellranger cols
    if( nrow(scRNAseq_metadata_sample) <2 ){scRNAseq_metadata_sample = scRNAseq_metadata } #no sample column found
    

    
    for (current_pheno in PHENOS_OF_INTEREST){ 
      #fetching cellular barcodes associated to pheno 
      colNums <- match( c( current_pheno),names(scRNAseq_metadata_sample))
      

      all_current_pheno_values = unlist( unique(scRNAseq_metadata_sample %>% dplyr::select(all_of(colNums)) )[,1]  )
      #print(all_current_pheno_values)
      if (   all(is.na( as.numeric(all_current_pheno_values)  )   ) == F    ){
        next # continuous numerical variable, impossible to plot all values 
      }
      
      for (current_pheno_value in all_current_pheno_values[!is.na(all_current_pheno_values)]){
      
      	cellular_barcodes =  scRNAseq_metadata_sample[which(scRNAseq_metadata_sample[,colNums] ==current_pheno_value  ),] 
      	cellular_barcodes = unique(as.vector( cellular_barcodes[,which(  cellular_barcodes %in% c("barcode", "cell_id") ) ] )) 

     
      	
      	 #if several libraries have been merged, a prefix or suffix must have been added. Here we get only the barcode sequence:
      	 cellular_barcodes_nuc = as.vector(str_match(cellular_barcodes, "[ATGC]{10,30}")) #length of barcode ~ 10-20nt 

      	 #get corresponding contig ID in original cell ranger data 
      	 contig_annot$cell_column = contig_annot[,which(names(contig_annot) %in% c("barcode", "cell_id") ) ][1:nrow(contig_annot)]
      	 contig_of_pheno = contig_annot %>% dplyr::rowwise() %>% dplyr::filter(cell_column %in% cellular_barcodes) #%>% dplyr::select(contig_id) )


      	 
      	 # subsetting Codon data frame
      	 subset_score_DF = Codon_Selection_Scores[which(Codon_Selection_Scores$SeqID %in%  as.character(contig_of_pheno$contig_id)  ),]

      	 
      	 
      	 # ploting and save 
      	 write.table(subset_score_DF, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "subset_score_DF.txt"), sep = "\t",  row.names = F, col.names =T) 

      	 current_pheno_value = str_replace_all(current_pheno_value, "/", "-") # Preventing nomeclatures featuring "/" symbols from generating bugs in automated file saving 
      	 if(nrow(subset_score_DF) > 1){
      	 
      	 selection_plot = plot_condon_scores(subset_score_DF, BCR_Regions, lines = T) +ggtitle( paste0("Selection scores of Seq. among  ", current_pheno_value) )
      	 ggsave( filename =  paste0( OUTPUT_PATH ,  "/", SIMU_TYPE,    "_local_scores_", current_pheno   , "_", current_pheno_value,  ".png"),
              plot =  selection_plot, width = 13, height = max(3, 3*length(unique(subset_score_DF$Chain)))  ) 
          #print( htmltools::tagList( list( ggplotly(selection_plot  )  )    ) )
          if(grepl("[Aa]nno|[Ss]ubset|[Pp]heno|[Tt]ype|[Tt]ime|LFC|BLI|[Cc]ycle|is", current_pheno ) | current_pheno_value == TRUE ){ #TO CHANGE 
          	j = j +1
          	temp_plot <- ggplotly(selection_plot , tooltip="text_plotly",  width = 1000, height = 1000 )
          	ggplotly_list[[j]] <-temp_plot
          	print(selection_plot)
          }
          
          
          
          
          #ggplotly_list[[j]] <- plotly_build(selection_plot)    
        }
      	 
              
      
      }
      
      
      
      
    }# end pheno loop  
    
    
  } #end if file 
  
} #end if directory


htmltools::tagList( ggplotly_list   ) 


#Ancestor
NCA_selection_plot = plot_condon_scores(Codon_Selection_Scores[which(Codon_Selection_Scores$SeqName =="Ancestor"),], BCR_Regions, lines = T) +ggtitle( paste0("NCA pattern") )

if(nrow(Codon_Selection_Scores[which(Codon_Selection_Scores$SeqName =="Ancestor"),]) > 0 ){ htmltools::tagList( ggplotly(NCA_selection_plot)   ) }
