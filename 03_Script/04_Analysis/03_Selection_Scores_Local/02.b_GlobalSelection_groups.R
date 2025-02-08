

## @knitr GLOBAL_Score_Computation_Comparisons

Selected_Trees = list.files(file.path( PATH_PHYLO, current_clono), pattern = "*Selected_Trees.txt", full.names =T, recursive = T)[1]
Selected_Trees = read.table(Selected_Trees, header = T, sep= " ")




################################################
# Cleaning data  
################################################
Codon_Selection_Scores = Codon_Selection_Scores %>% dplyr::mutate(SelectScore = ifelse(is.infinite(abs(SelectScore)), NA, SelectScore ) )


################################################
# Computing Instantaneous Positive, Negative, Total Selection Scores 
################################################

#Positive and Negative 

ColumnsToGroupby <- c("SeqID", "SeqName", "Chain", "ParentNode", "MutSinceParent", "NmutSinceNCA","Isotype", "Dist_ToGermline", "Abundance", "MutSinceParent_NS", "Sampled_Ancestor") # "MutSinceParent_NS","NmutSinceNCA_NS", 
Sequence_Selection_Scores = left_join( distinct(Codon_Selection_Scores[,which(names(Codon_Selection_Scores) %in% ColumnsToGroupby)] ), 
					ComputeSelectScores(Codon_Selection_Scores) %>% distinct() ) %>% 	#Immediate Selection Scores 
					dplyr::group_by(SeqName, Chain) %>% dplyr::mutate(		#Normalize to Number of Nucleotides involved
							SelectScore_Inst_Pos 	= SelectScore_Inst_Pos, #/log2(1+MutSinceParent)
							SelectScore_Inst_Neg 	= SelectScore_Inst_Neg  

					) 




################################################
# Cumulated selection score 
################################################  

Lineage_Cumulated_Scores   = dplyr::left_join(Sequence_Selection_Scores, Selected_Trees %>% dplyr::select(Chain, Tree_Nk)  ) 


Lineage_Cumulated_Scores = Lineage_Cumulated_Scores  %>% dplyr::rowwise() %>%  dplyr::mutate(NodeParents = as.character(mcmapply(getParents, Tree_Nk, SeqName, SIMPLIFY = FALSE) ) )

 

Lineage_Cumulated_Scores = 	ComputeCumulatedScores(Lineage_Cumulated_Scores)
Sequence_Selection_Scores = left_join( distinct(Sequence_Selection_Scores), distinct(Lineage_Cumulated_Scores) ) 	 %>% dplyr::mutate(	Cell_Group  = 'Total')
Sequence_Selection_Scores[sapply(Sequence_Selection_Scores, is.infinite)] <- NA

Time_Df =  read.csv( text =  paste( c( "Clono_Chain", "SeqName", "Time"  ) , collapse = ",")  ) 

################################################
# Ploting results 
################################################
#to plot trees
 SeqName_Features = Sequence_Selection_Scores %>% ungroup()  %>% distinct()
 
 

# Linking sequences to phenos of interest 
if(length(na.omit(scRNAseq_metadata_path)) > 0 ){ #if metadata found 
  tokeep = names(scRNAseq_metadata_sample)[which( names(scRNAseq_metadata_sample) %in% names(contig_annot)==F | names(scRNAseq_metadata_sample) %in% c("barcode", "cell_id"))  ]
  Selected_PhenoMetadata = inner_join(contig_annot , scRNAseq_metadata_sample[tokeep  ]) %>% dplyr::rename( SeqID= "contig_id") # check with other cell ranger versions cell_id (by)
  print(dim(Selected_PhenoMetadata))
  


  
  
  

  All_Select_Scores =  grep("SelectScore", colnames(Sequence_Selection_Scores), value = T)
  cat("\n\n")
     
  for (current_pheno in   PHENOS_OF_INTEREST  ){
  
	cat("## ", current_pheno, "{.tabset}", "\n\n") 
	
	
	#Fetch right metadata  -------------
	colNums <- match( c(current_pheno, "SeqID", "barcode", "cell_id") , names(Selected_PhenoMetadata))
	colNums = colNums[which(is.na(colNums) ==F)]
	 
    
    	current_pheno_to_add = Selected_PhenoMetadata %>% dplyr::select(colNums) %>%
    				dplyr::filter(SeqID %in% SeqName_Features$SeqID ) 
    	current_pheno_to_add = left_join(Sequence_Selection_Scores,  current_pheno_to_add  %>% distinct() )
    	
    	SeqName_Features = left_join(SeqName_Features, Selected_PhenoMetadata %>% dplyr::select(colNums) %>% distinct() ) %>% dplyr::mutate(PhenoName = as.factor(as.character(get(current_pheno )) ))  %>%
    	  			mutate_at(vars(contains("SelectScore")), funs(as.numeric)) 
    				

    	#Store Time Data (for aa analysis) -------------

    	if( grepl("[Tt]ime",current_pheno )  ){  
    		Time_Df =  SeqName_Features  %>% ungroup() %>% dplyr::mutate(Clono_Chain = paste0(current_clono , "_", Chain) , Time = current_pheno ) %>% 
    				dplyr::select(SeqName, Clono_Chain, Time) %>% dplyr::filter( is.na(Time) ==F ) %>% as.data.frame() 
     	}
    
 
    				
	#Ploting within pheno of interest comparisons for current score  -------------
	if(  all( is.na(  SeqName_Features$PhenoName  )) ==F &  length( unique(SeqName_Features$PhenoName  ) ) > 1     ){
	
	
		for(current_score in All_Select_Scores ){
	
		cat("### ", current_score, "\n\n")
	
		
		
		#if( nrow(distinct( SeqName_Features[,  which(names(SeqName_Features) %in% c(current_score, current_pheno )  )])) <2 | all(is.na(SeqName_Features$PhenoName)) ){next}
		SeqCounts = SeqName_Features %>% group_by(Chain, PhenoName) %>% dplyr::summarize( n = n()  ) 
	
		current_plot  = ggplot(data = SeqName_Features %>% dplyr::filter(grepl("X[0-9]", SeqName)==F) %>% dplyr::mutate(Score = round( as.numeric( get(current_score) ),5   ) ), 
				aes(x = PhenoName, y =Score )) + 
    				#geom_violin(aes(fill = as.factor(PhenoName)), alpha = 0.2, color = "transparent") +
    				geom_jitter(aes(color = as.numeric(MutSinceParent)), width = 0.1, size = 0.75) +  facet_grid(Chain ~ .) +
    				theme_bw() +  theme(legend.position = "none") + 
    				 theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle(as.character(current_score))
    		if( 1 %in% c(SeqCounts$n) ==F ){current_plot = current_plot + geom_violin(aes(fill = as.factor(PhenoName)), alpha = 0.2, color = "transparent") }
    		 
    				 
    				 
	
		print(current_plot)
	
		 
	
		cat("\n\n")
	
		}
	
		
	
	}
	
	
  }
	
		
}	



################################################
# Saving  node features 
################################################
SeqName_Features_to_Save = SeqName_Features %>% dplyr::mutate(across(!SeqName, as.character)) %>% dplyr::filter(grepl("X[0-9]", SeqName)==F)

SeqName_Features_to_Save[is.na(SeqName_Features_to_Save)] <- "NaN" #for python bug '<' not supported between instances of 'str' and 'float' (colouring tree nodes with ete3) 

SeqName_Features_to_Save = SeqName_Features_to_Save[,which(names(SeqName_Features) %in% c(PHENOS_OF_INTEREST, All_Select_Scores, "SeqID" , "SeqName" ,"Chain" , "MutSinceParent", "ParentNode", "barcode" , "NmutSinceNCA", "Isotype", "Dist_ToGermline", "Abundance", "Sampled_Ancestor", "MutSinceParent_NS" )) ]  %>% #, "MutSinceParent_NS", "NmutSinceNCA_NS"
mutate_at(vars(contains("SelectScore")), funs(as.numeric)) 


SeqName_Features_to_Save = dplyr::left_join(SeqName_Features_to_Save, distinct(Real_Nodes_List)) %>% dplyr::filter(Inferred ==F) %>% dplyr::select( -c(Inferred) )


SeqName_Features_to_Save$Sample =  basename(PATH_PHYLO)
SeqName_Features_to_Save$Clonotype =  current_clono
write.table( distinct(SeqName_Features_to_Save ) ,     paste0( PATH_PHYLO , "/", current_clono, "/", SIMU_TYPE,    "_NodeFeatures_treePlot",   ".txt") ,  sep = " ") #%>% dplyr::select(-SeqID)
print(paste0( PATH_PHYLO , "/", current_clono, "/", SIMU_TYPE,    "_NodeFeatures_treePlot",   ".txt"))



################################################
# Saving  aa changes data 
################################################
PATTERN_TABLE = PATTERN_TABLE %>% dplyr::select(Regex, Size, Properties) #, Color
beforemut_df =PATTERN_TABLE  %>% dplyr::rename( aa.beforeMut = Regex, properties.beforeMut =  Properties,  size.beforeMut = Size  )   #, color.beforeMut = Color
aftermut_df =PATTERN_TABLE %>% dplyr::rename( aa.afterMut = Regex, properties.afterMut =  Properties,  size.afterMut = Size )   #, color.afterMut = Color 


Codon_Selection_Scores$Sample = basename(dirname( dirname( OUTPUT_PATH) ) )

AA_Properties_toSave = Codon_Selection_Scores %>%  dplyr::filter(Inferred ==F) %>% 
                          dplyr::group_by(SeqName, Chain  ) %>% dplyr::mutate(
                          				aa.beforeMut = substr(Parent_ProtSeq, CodonIndex, CodonIndex)  , 
                                                      aa.afterMut  = substr(ProtSeq, CodonIndex, CodonIndex) ,
                                                      aa.NCA = substr(NCAProt, CodonIndex, CodonIndex) ,
                                                      aa.Germline = substr(GermlineProt, CodonIndex, CodonIndex),
                                                      is.mutated = ifelse(aa.beforeMut != aa.afterMut, T, F ),
                                                      InheritedFrom = case_when(is.mutated == T ~ "NaN",
                                                      			   aa.beforeMut == aa.NCA & aa.beforeMut != aa.Germline  ~ "fromNCA",
                        							    aa.beforeMut != aa.NCA & aa.beforeMut != aa.Germline ~ "afterNCA",
                        							    aa.beforeMut == aa.Germline ~ "Germline", #, InheritedFromGermline == T ~ "Germline"
                        							    TRUE ~ "undefined" ), 
                        				Dist_Germline_to_NCA = Dist_ToGermline - NmutSinceNCA
                                                      ) %>%  #dplyr::filter(aa.beforeMut != aa.afterMut) %>%
                        dplyr::select(SeqName, Chain, aa.beforeMut, aa.afterMut , SelectScore, Clono_Chain, Sample, Region, Chain, InheritedFrom, 
                        		aa.Germline, aa.NCA , Dist_Germline_to_NCA) %>% distinct()
AA_Properties_toSave = left_join( left_join(AA_Properties_toSave, distinct(beforemut_df) )   ,  distinct(aftermut_df) )      

#write.table(AA_Properties_toSave, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "_mutated_aa_properties0.txt"), sep = "\t",  row.names = F, col.names =T) 


if(nrow(Time_Df)> 0 ){ 
 AA_Properties_toSave = dplyr::left_join( as.data.frame(AA_Properties_toSave) %>% distinct(),  Time_Df %>% distinct(), multiple = "first"  )
 AA_Properties_toSave$Time = paste0(AA_Properties_toSave$Time)
 AA_Properties_toSave = AA_Properties_toSave %>% dplyr::mutate(across(-c(SelectScore), as.character ))
  }



             

write.table(Time_Df, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "Time_Df.txt"), sep = ",",  row.names = F, col.names =T) 



write.table(AA_Properties_toSave %>% dplyr::select(everything()), file = paste0( OUTPUT_PATH, "/", "mutated_aa_properties.txt"), sep = ",",   col.names =T) #row.names = F,




