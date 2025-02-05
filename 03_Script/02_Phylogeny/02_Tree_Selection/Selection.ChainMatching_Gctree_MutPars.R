# Selecting one tree or a heavy/light pair 
# Based on cChain concordance, Mut Parisomony, Gctree criteria 
# -----------------------

## @knitr ChainMatching_Gctree_MutPars






#Updating available chains (in case no trees could pass the filters in the previous steps)
available_chains_updated = unique(Tree_Nk_Table[which(Tree_Nk_Table$toKeep == T),]$Chain)
if ( length(available_chains) > length(available_chains_updated)){
  print(paste0("only ", available_chains_updated[1], " chain trees could pass the filters"))
  available_chains = available_chains_updated
}




while( length(available_chains)>1 ){

  heavy_confidence = 1- (quantileH/(quantileH+quantileL)) #trees we least incertainty are preferred (gaps in nmuts)
   print(paste0("heavy confidence = ", heavy_confidence))
     
   
  # Wide to long table and gathering common columns to heavy and light data. (unique pairs BCR are taken )

  cell_seq_IDs = airr_table  %>% dplyr::select(-contains( c("sequence", "Isotype")     )) %>% na.omit() %>% tidyr::gather(Chain, SeqName, -cell_id) %>%
                      dplyr::mutate(Chain = str_to_title(str_extract(Chain, "[^_]+") ),
                             SeqNameChain = paste0(Chain, SeqName) ) %>% arrange(cell_id) #%>% #%>% group_by(SeqNameChain) %>% slice(1) 
                             
  
  #write.table(cell_seq_IDs, file = file.path(GCTREE_FOLDER, "cell_seq_IDs1.txt") ) 
  
  # Take into consideration cell quality in barcode information ------------------------------------------
   scRNAseq_metadata_path = list.files(path = CELL_METADATA_FOLDER,  pattern = "etadata|ETADATA|etaData", full.names =T)[1]
   metadata_df = read.table(file = scRNAseq_metadata_path, header = T, sep =  detect_separator(scRNAseq_metadata_path) )   
   
   write.table(metadata_df, file = file.path(GCTREE_FOLDER, "metadata_df.txt") )
   
   
   metadata_cell =  names(metadata_df)[which( names(metadata_df) %in% c("barcode", "cell_id")  )]
   metadata_cell = metadata_df[, c(metadata_cell)]
   cell_seq_IDs =  cell_seq_IDs %>% dplyr::mutate( cell_id_inMetadata = ifelse(  cell_id %in% c(rownames(metadata_df), unlist(metadata_cell) ) , 
   					cell_id, paste0(cell_id, Chain) )     ) 
   					
   					
                       
   # Chain pairs in good quality cells 
    cell_seq_IDs = cell_seq_IDs %>%                         
                      group_by(cell_id_inMetadata) %>% dplyr::mutate(BCRcombi = paste0( unique( paste0(Chain, SeqName)  ), collapse = ".")) #paste0(Chain, str_to_title(SeqName))  )
      

   # if no combis available => keep one chain only -----------------------
   paired_chains = cell_seq_IDs$BCRcombi[which( grepl(  "\\.", cell_seq_IDs$BCRcombi  ) ) ]
   if(length(paired_chains)  < 1 ){
   	available_chains = 1
   	if(heavy_confidence <  0.5 ){Tree_Nk_Table = Tree_Nk_Table %>% dplyr::filter( Chain == "Heavy" ) 
   	}else{  Tree_Nk_Table = Tree_Nk_Table %>% dplyr::filter( Chain == "Light" )     } 
   	break #cannot break out of if statement => while loop with 1 iteration 
   }
   
   #---------------------------------------------------------------------
   
   
                 

  
   #Here we select unique bcr combinations instead of cells to avoid repeated values to be taken for correlation (to be addded?)
   #Reasonning: most expanded clones would mask the contribution of unique combinations to the correlation score 
  cell_seq_IDs = cell_seq_IDs %>% ungroup() %>% group_by(Chain, BCRcombi) %>% dplyr::arrange(cell_id) %>%  dplyr::slice(1)
  
  #write.table(cell_seq_IDs, file = file.path(GCTREE_FOLDER, "cell_seq_IDs4.txt") )
  
  cloneSize = as.numeric( gsub(".*CLONOTYPE_n(.+)_Id[0-9].*", "\\1", GCTREE_FOLDER ) ) 
  #if(cloneSize > 200 ){Tree_Nk_Table = Tree_Nk_Table %>% dplyr::filter(Pars_Penalty < 0 | Chain == "Heavy"  ) }



  # Repeat with heavy and light chains 
  chain_combis_tree_unique = left_join(cell_seq_IDs %>% dplyr::select(-SeqNameChain), 
                      Tree_Nk_Table %>%  dplyr::filter(toKeep ==T) %>% #dplyr::filter(  (Tree_Index == 6212 & Chain == "Heavy" ) | (Chain == "Light"  )  ) %>% 
                      dplyr::select(Chain,Tree_Index,Tree_Nk)  )  %>% 
                      dplyr::mutate(across(-c(Tree_Index), as.character))

  


  # mutate columns node position , node parent list 
  print("inferring parents")
  #chain_combis_tree_unique = chain_combis_tree_unique %>% ungroup() %>%  dplyr::filter(  grepl(pattern = SeqName, Tree_Nk)  ) #New: seqnames dropped by Gctree if equal to germline 
  
  #write.table(chain_combis_tree_unique, file = file.path(GCTREE_FOLDER, "chain_combis_tree_unique0.txt") )
  
  chain_combis_tree_unique = chain_combis_tree_unique %>% dplyr::rowwise() %>%  dplyr::filter(  grepl(pattern = paste0(SeqName, ":"), Tree_Nk)  ) #New: seqnames dropped by Gctree if equal to germline 
  chain_combis_tree_unique = chain_combis_tree_unique  %>% dplyr::mutate(NodeParents =mcmapply(getParents, Tree_Nk, SeqName, SIMPLIFY = FALSE ) ) #
  chain_combis_tree_unique$NodeParents = as.character(chain_combis_tree_unique$NodeParents)
  #write.table(chain_combis_tree_unique, file = file.path(GCTREE_FOLDER, "chain_combis_tree_unique1.txt") )
  #chain_combis_tree_unique = chain_combis_tree_unique  %>% dplyr::mutate(NodeDaughters =mcmapply(getDaughters, Tree_Nk, SeqName, SIMPLIFY = FALSE) )
 
  chain_combis_tree_unique  = chain_combis_tree_unique %>% dplyr::mutate(NodePosition = 1+ stringr::str_count( NodeParents, ','), #
                             NodeParents = str_remove_all(NodeParents, ",[0-9]{1,3}"),
                             #NodeParents = paste0(NodeParents, ",",  SeqName ), #NEW: including 
                             NodeParents = str_replace_all(NodeParents, "seq", paste0(Chain, "seq") ) , #%>%   #unambiguous seqName
                             #New
                             #NodeDaughters = str_remove_all(NodeDaughters, ",[0-9]{1,3}"),
                             #NodeDaughters = str_replace_all(NodeDaughters, "seq", paste0(Chain, "seq") ) )
                             terminal_leaves = paste0(str_extract_all(Tree_Nk, pattern = "(?<=[(,])seq[0-9]{1,3}")[[1]], collapse ="|")  , #seq following ( or , in nk 
                             is_terminal =   grepl(terminal_leaves, SeqName) ,
                             LastParent = str_extract(NodeParents , pattern = "([^,*]+$)")  )
                      #dplyr::select(-c(SharedSeq, Tree_Nk)) 

  #if nca corresponds to real sequence, it is renamed as such NEW
  #chain_combis_tree_unique = chain_combis_tree_unique %>% group_by(Chain) %>%
  #	dplyr::mutate(debug_column =  ifelse(NodeParents == "" ,"Ancestor", str_replace(NodeParents, "Ancestor", "Seq0") )  ) %>%
  #	dplyr::mutate(NodeParents = case_when( any( NodeParents == "" ) ~debug_column , TRUE ~ NodeParents )) %>%
  #	dplyr::select(-debug_column)
  write.table(chain_combis_tree_unique, file = file.path(GCTREE_FOLDER, "chain_combis_tree_unique.txt") )

  #Dictionary for matching BCR pairs 
  lookup_table = chain_combis_tree_unique  %>% ungroup() %>% #dplyr::mutate(#SeqName = paste0(Chain, SeqName ), 
  						#		LastParent = str_extract(NodeParents , pattern = "([^,*]+$)")  ) %>% #after last ,
  			#dplyr::select(BCRcombi, LastParent, Chain) %>% spread( key = Chain, value =  LastParent) %>% rename(HeavyParent = Heavy, LightParent = Ligth) %>% 
  			dplyr::select(BCRcombi) %>% #, HeavyParent, LightParent
                 	dplyr::mutate(HeavySeq = sub("\\..*", "", BCRcombi), LightSeq = sub(".*\\.", "", BCRcombi)) %>% distinct()
                                                       

  
  #Computing correlation position + jacquard index parents 
  heavy_table = chain_combis_tree_unique %>% dplyr::filter(Chain == "Heavy") %>% dplyr::select(BCRcombi, Tree_Index, Tree_Nk, NodeParents, NodePosition, is_terminal,LastParent) %>% 
   				dplyr::rename_with(~ paste0("Heavy.", .), -BCRcombi)
  light_table = chain_combis_tree_unique %>% dplyr::filter(Chain == "Light") %>% dplyr::select(BCRcombi, Tree_Index, Tree_Nk, NodeParents, NodePosition, is_terminal,LastParent) %>% 
  				dplyr::rename_with(~ paste0("Light.", .), -BCRcombi)
  
  chain_combis_tree_pair = merge(heavy_table, light_table, by ="BCRcombi", all.x =TRUE, all.y = TRUE)
   
  
  ################TO REMOVZ AFTER
  write.table(light_table, file = file.path(GCTREE_FOLDER, "light_table.txt") , row.names = T)
  write.table(lookup_table, file = file.path(GCTREE_FOLDER, "lookup_table.txt") , row.names = T)


  
  #######################
  
  # filtering out tree pairs with uncorrelated leaf positiions (to save time in the next step)

   
  chain_combis_tree_pair = chain_combis_tree_pair  %>% dplyr::filter(Heavy.is_terminal  | Light.is_terminal )
  

  
  chain_combis_tree_pair = chain_combis_tree_pair %>% group_by(Heavy.Tree_Index, Light.Tree_Index) %>%
                           dplyr::mutate(CorCoeff = stats::cor(Heavy.NodePosition, Light.NodePosition,  method = "spearman") ) 
                           

   hist1 = ggplot(chain_combis_tree_pair, aes(x=CorCoeff)) + geom_histogram(color="darkblue", fill="lightblue", binwidth=0.01) + ggtitle("Position Correlation")                         
   
 
   filter_cor_leaves = min( tail(  sort(chain_combis_tree_pair$CorCoeff) ,  1000000 ) )
   
   if( nrow(chain_combis_tree_pair)  >  150000 ){ filter_cor_leaves =  -0.7 #quantile(  chain_combis_tree_pair$CorCoeff[which(chain_combis_tree_pair$CorCoeff > -0.1 ) ]   , 0.05 )  
   }else if(nrow(chain_combis_tree_pair)  >  50000){filter_cor_leaves  = -0.9  } 

    
   if( all(is.na(  filter_cor_leaves   ))   ){  filter_cor_leaves = -1
   }else{filter_cor_leaves = max( c(-1, filter_cor_leaves), na.rm = T)} #0.1 too stringent for small trees and few data 
   #filter_cor_leaves = max( c(0, filter_cor_leaves), na.rm = T) 
   #keep only positively correlated or threshold to reduce computing time. is.na for small trees with no lineage
   
  




   print("prefiltering on leaf position correlation")
   print(dim(chain_combis_tree_pair))
   
   #Leaves position are positively correlated #keeping terminal leaves for next criteria 
   chain_combis_tree_pair = chain_combis_tree_pair   %>% dplyr::filter( CorCoeff >= filter_cor_leaves | is.na(CorCoeff )    ) # &  (Heavy.is_terminal  | Light.is_terminal )	
   print(dim(chain_combis_tree_pair))

  if(nrow(chain_combis_tree_pair) < 1){
  	#sometimes no consistent trees are available: in this case, the script is stopped and the status will be caught in the bash environement
  	q(save = "no", status = 0, runLast = FALSE)
  	#we are obliged to set exit status to 0, otherwise not supported by snakemake
  	#stop("Node positions are not consistent between chains") 
  }
  
  print(paste0( "correlation leaf position threshold: ,  ", filter_cor_leaves) )

  
  # concordance of lineages 
  print("computing jaccard on lineage")

  
  
  

  

  
  lookup_table = lookup_table[which( grepl(".", lookup_table$BCRcombi )  ),]  %>% dplyr::mutate(across(everything(), as.character)) %>% mutate_all( tolower )  #only pairs
  lookup_table = rbind(lookup_table, c( "heavyancestor.lightancestor", "heavyancestor", "lightancestor" ) )   # adding expected shared ancestor
  
  write.table(lookup_table, file = file.path(GCTREE_FOLDER, "lookup_table.txt") , row.names = T)
  
  
  parent_df = chain_combis_tree_pair %>% dplyr::mutate(HeavySeq = sub("\\..*", "", BCRcombi), LightSeq = sub(".*\\.", "", BCRcombi),
                                                     Heavy.NodeParents = str_replace_all(Heavy.NodeParents , ",", "/" ),  Light.NodeParents = str_replace_all(Light.NodeParents , ",", "/" ), 
                                                     Heavy.LastParent = basename(dirname( Heavy.NodeParents)),  Light.LastParent = basename(dirname( Light.NodeParents) ),
                                                     Heavy.GrandParent =  basename(dirname( dirname(Heavy.NodeParents))),  Light.GrandParent = basename(dirname( dirname(Light.NodeParents)) ) 
                                                     )    %>% 
                dplyr::select(BCRcombi, Heavy.Tree_Index, Light.Tree_Index, Heavy.LastParent, Light.LastParent, Heavy.GrandParent,Light.GrandParent)  %>% distinct() 
  
  
  parent_df = parent_df %>% dplyr::mutate( Light.LastParent = ifelse(Light.LastParent == "Ancestor", "LightAncestor", Light.LastParent ) ,
                                         Heavy.LastParent = ifelse(Heavy.LastParent == "Ancestor", "HeavyAncestor",  Heavy.LastParent ),
                                         Light.GrandParent = ifelse(Light.GrandParent == "Ancestor", "LightAncestor", Light.GrandParent ) ,
                                         Heavy.GrandParent = ifelse(Heavy.GrandParent == "Ancestor", "HeavyAncestor", Heavy.GrandParent ),
                                          ) %>% 
                            dplyr::mutate(across(where(is.character) , tolower  )) %>%  #mutate_all( tolower )
                            dplyr::mutate(HeavySeq = sub("\\..*", "", BCRcombi), LightSeq = sub(".*\\.", "", BCRcombi)) %>% distinct() #NEW
                            
  write.table(parent_df, file = file.path(GCTREE_FOLDER, "parent_df.txt") , row.names = T)                          
                            

    # Fitlering only pairs (now that unpaired data are in parent df) and terminal leaves for concordance criterium 
    
   paired_terminal = chain_combis_tree_pair[which(chain_combis_tree_pair$Heavy.is_terminal & chain_combis_tree_pair$Light.is_terminal   ),]$BCRcombi
   paired_terminal =  paired_terminal[grepl( "\\.", paired_terminal )] 

    if( length(paired_terminal   ) > 0    ){ #bugs for small clones where terminal leaves unpaired 
  	chain_combis_tree_pair = chain_combis_tree_pair   %>% dplyr::filter(    (Heavy.is_terminal| is.na(Heavy.is_terminal ) )  & 
   										(Light.is_terminal | is.na(Light.is_terminal ))  #& grepl("\\." , BCRcombi ) 
   									)
   	       
  	
  }
  rm(paired_terminal)
  gc()   
               
   chain_combis_tree_pair = chain_combis_tree_pair %>% dplyr::filter( (is.na(Heavy.is_terminal ) & is.na(Light.is_terminal ) ) ==F  )  
   if( cloneSize > 10 ){chain_combis_tree_pair = chain_combis_tree_pair %>% dplyr::filter(grepl("\\." , BCRcombi )    ) }     
   #write.table(chain_combis_tree_pair, file = file.path(GCTREE_FOLDER, "chain_combis_tree_pair_prefilter.txt") )
   
   if(cloneSize > 1 ){  n_cores = 40 }
   #if( cloneSize > 100  ){  n_cores = 65 }
   
   chain_combis_tree_pair$Jaccard_Index = mcmapply(compare_lineage, 
  						chain_combis_tree_pair$Heavy.NodeParents, chain_combis_tree_pair$Light.NodeParents,
  						chain_combis_tree_pair$Heavy.Tree_Index, chain_combis_tree_pair$Light.Tree_Index, 
  						#chain_combis_tree_pair$Heavy.NodeDaughters, chain_combis_tree_pair$Light.NodeDaughters,  
  						MoreArgs = list(parent_df, lookup_table) , mc.cores = n_cores  ) #LONG!!!!#, mc.cores = n_cores
  						 
  #write.table(chain_combis_tree_pair, file = file.path(GCTREE_FOLDER, "chain_combis_tree_pair_jaccard.txt") )	
  
   # Remove NAs 
  tree_pair_coeff_filtered = chain_combis_tree_pair %>% ungroup() %>% dplyr::filter(grepl("\\." , BCRcombi ) )
  
  					

  #***************************************************************************
  # NEW: checking isotype order 
     

  CSR_df = airr_table %>% dplyr::select( -contains("sequence_id") ) %>% pivot_longer(cols = -c("cell_id") , 
                                      names_to = c( "Chain", ".value" ), 
                                      names_sep="_" ) %>% dplyr::mutate( Chain = str_to_title(Chain) ) 
  
  
  
  CSR_df = left_join(chain_combis_tree_unique, CSR_df ) %>% dplyr::group_by( Chain ,   Tree_Nk, NodeParents, SeqName) %>% dplyr::summarise( Isotype = paste(unique(Isotype), collapse=",")     ) 

  
  
  
  
  CSR_df = CSR_df %>% dplyr::mutate( #last_parent = stringi::stri_extract_last(NodeParents, regex = "seq[0-9]{1,3}|Ancestor"  ),
   #                                #clono_chain = paste0(Clonotype, Chain),
   				    	last_parent = str_remove(NodeParents, SeqName), 
                                   	last_parent = stringi::stri_extract_last(last_parent, regex = "Ancestor|seq[0-9]{1,3}"  ),
                                   Isotype = str_remove(Isotype, pattern  = "NA,|,NA")        ) 

  
   
  
  
  CSR_df$Parent_Isotype = CSR_df[match(CSR_df$last_parent, CSR_df$SeqName ), ]$Isotype 

  
  CSR_df$CSR_discordance = mcmapply(Check_Isotype_Order, CSR_df$Isotype, CSR_df$Parent_Isotype, 
  				MoreArgs = list( c(  "IgD", "IgM",  "IgG3", "IgG1", "IgA1", "IgG2", "IgG2B", "IgG2A", "IgG4", "IgE", "IgA" , "IgA2" )   ), 
  				# c(  "IgD", "IgM",  "IgG3", "IgG1", "IgG2", "IgG2B", "IgG2A", "IgE", "IgA"  )
  				mc.cores = n_cores   ) 
  

  
  summary_CSR_scores = CSR_df %>% dplyr::group_by( Chain, Tree_Nk ) %>% dplyr::summarise(  n_CSR_discordances = sum(CSR_discordance, na.rm = T)   ) #%>% dplyr::mutate( Tree_Nk = as.character(Tree_Nk)   ) 

   
   summary_CSR_scores_heavy =  summary_CSR_scores %>% dplyr::filter(Chain == "Heavy") %>% dplyr::select( - Chain) %>% dplyr::rename_with(~ paste0("Heavy.", .) )
   summary_CSR_scores_light =  summary_CSR_scores %>% dplyr::filter(Chain == "Light") %>% dplyr::select( - Chain) %>% dplyr::rename_with(~ paste0("Light.", .) )



   #***************************************************************************
  

  write.table(chain_combis_tree_pair, file = file.path(GCTREE_FOLDER, "chain_combis_tree_pair.txt") , row.names = T)
  chain_combis_tree_pair = chain_combis_tree_pair  %>% group_by(Heavy.Tree_Index, Light.Tree_Index) %>%
                           dplyr::summarise(
                           		CorCoeff = mean(CorCoeff, na.rm = T), #stats::cor(Heavy.NodePosition, Light.NodePosition,  method = "spearman") #mean doesn't mattter unique value
                                       MeanJaccard = mean(Jaccard_Index, na.rm = T),
                                       n_lineage_discordances = -trunc(sum(Jaccard_Index, na.rm = T )/200) ) %>% 
                           dplyr::mutate(across(everything(), as.numeric)) 
  


   print("Generating outputs")

   hist1 = hist1 + geom_vline(xintercept = filter_cor_leaves, linetype="dashed", color = "red", size=1)
   hist2 = ggplot(chain_combis_tree_pair, aes(x=MeanJaccard)) + geom_histogram(color="tan4", fill="tan3", binwidth=0.01) + ggtitle("Lineage Matching")  #+ xlim(c(0,1)
   histoplots = cowplot::plot_grid(hist1, hist2, labels=NULL, ncol = 2, nrow = 1)
   
   #keeping track of max concordance score for tree quality filtering 
   max_concordance_score = max(0, chain_combis_tree_pair$MeanJaccard, na.rm = T)
   max_concordance_score = max( chain_combis_tree_pair$MeanJaccard, na.rm = T)
   
   # Apply last filter
   # Computing final score on pairs of trees 
   
   
     tree_pair_coeff_filtered = chain_combis_tree_pair %>% 
    left_join(Tree_Nk_Table_chains[[1]] %>% setNames(paste0('Heavy.', names(.)))  ) %>% 
    left_join(Tree_Nk_Table_chains[[2]] %>% setNames(paste0('Light.', names(.)))  )
  
  tree_pair_coeff_filtered = type.convert(as.data.frame(tree_pair_coeff_filtered) , as.is = TRUE) #some columns were converted to character
  #non parsimonious trees are discarded if alternative  + selection based on gctree scores
  
  
  
  tree_pair_coeff_filtered = left_join( left_join(tree_pair_coeff_filtered , summary_CSR_scores_heavy    ), summary_CSR_scores_light )  %>% ungroup() %>% 
  				dplyr::filter( (is.na(Heavy.Tree_Index) | is.na(Light.Tree_Index) ) == F  ) #remove pairs of NAs
  
  write.table(tree_pair_coeff_filtered , file = file.path(GCTREE_FOLDER, "tree_pair_coeff_filtered.txt") , row.names = T)
  
  
  
   if( all(is.na(tree_pair_coeff_filtered$MeanJaccard) ) ){Selected_Trees = tree_pair_coeff_filtered 
  }else{
  
  
  #if no pair is concordant, only one chain is kept (the one with best confidence ie less mutation gaps) 
  #if(nrow(tree_pair_coeff_concordant) < 1){
  #	if(heavy_confidence < 0.5){  tree_pair_coeff_filtered[,which(grepl("Heavy|MeanJaccard", names(tree_pair_coeff_filtered))   )] <- NA
  #	}else{ tree_pair_coeff_filtered[,which(grepl("Light|MeanJaccard", names(tree_pair_coeff_filtered))   )] <- NA   }
  #}else{ tree_pair_coeff_filtered = tree_pair_coeff_concordant } 
  
  

  

  tree_pair_coeff_filtered = tree_pair_coeff_filtered  %>% 
  		dplyr::filter(n_lineage_discordances == min(n_lineage_discordances, na.rm = T))  %>% dplyr::filter(MeanJaccard >= quantile(MeanJaccard, probs=0.5, na.rm = T)) 
  
  # if none of the trees are concordant => only chain with highest confidence is kept
  ##if( nrow(tree_pair_coeff_filtered) < 1 ){
  ## 	available_chains = 1
  ## 	if(heavy_confidence <  0.5 ){Tree_Nk_Table = Tree_Nk_Table %>% dplyr::filter( Chain == "Heavy" ) 
  ## 	}else{  Tree_Nk_Table = Tree_Nk_Table %>% dplyr::filter( Chain == "Light" )     } 
  ## 	break #cannot break out of if statement => while loop with 1 iteration 
  ## }
  
  
  # Remaining fitlers on concordant trees (CSR, overlap optimization, GcTree) 
  if( heavy_confidence < 0.5){
  	tree_pair_coeff_filtered = tree_pair_coeff_filtered %>% 
  		dplyr::filter(Light.n_CSR_discordances == min(Light.n_CSR_discordances, na.rm = T) ) %>% dplyr::filter(Heavy.n_CSR_discordances == min(Heavy.n_CSR_discordances , na.rm=T)  )  %>% 
  		dplyr::filter(MeanJaccard >= quantile(MeanJaccard, probs = 0.99, na.rm = T) ) %>% 
  		dplyr::filter( Light.Pars_Penalty == min(Light.Pars_Penalty, na.rm = T) ) %>%  dplyr::filter(Heavy.Pars_Penalty == min(Heavy.Pars_Penalty, na.rm = T )  )
  }else{
  	tree_pair_coeff_filtered = tree_pair_coeff_filtered %>% 
  		dplyr::filter(Heavy.n_CSR_discordances == min(Heavy.n_CSR_discordances, na.rm = T) ) %>% dplyr::filter(Light.n_CSR_discordances == min(Light.n_CSR_discordances, na.rm = T )  )  %>% 
  		dplyr::filter(MeanJaccard >= quantile(MeanJaccard, probs = 0.99, na.rm = T) ) %>% 
  		dplyr::filter( Heavy.Pars_Penalty == min(Heavy.Pars_Penalty,na.rm = T) ) %>%  dplyr::filter(Light.Pars_Penalty == min(Light.Pars_Penalty , na.rm = T)  )
  }
  
  tree_pair_coeff_filtered = tree_pair_coeff_filtered %>% 
    			dplyr::mutate(Final_Score = heavy_confidence*( (gctree_contrib*abs_med_div(Heavy.GcTree_ll) ) + 
                                                   mutpars_contrib*abs_med_div(Heavy.Mut_Pars)  +
                                                   as.numeric(Heavy.Pars_Penalty)   ) + 
                                (1-heavy_confidence)*( (gctree_contrib*abs_med_div(Light.GcTree_ll)) + 
                                                        mutpars_contrib*abs_med_div(Light.Mut_Pars)  +
                                                        as.numeric(Light.Pars_Penalty)   ) ) 
  
  
	
  
  tree_pair_coeff_filtered = tree_pair_coeff_filtered %>% dplyr::mutate(
  				 Heavy.n_lineage_discordances = n_lineage_discordances, Light.n_lineage_discordances = n_lineage_discordances ,
  				 Heavy.MeanJaccard = MeanJaccard, Light.MeanJaccard = MeanJaccard ) 
  
  
  Selected_Trees =tree_pair_coeff_filtered %>%  dplyr::arrange( Final_Score  ) %>% dplyr::slice(1) %>% #         dplyr::select(where(is.numeric)) %>% 
    dplyr::select(starts_with(c("Heavy", "Light"))) %>%  #, "n_lineage_discordances", "MeanJaccard", "all_aberrations
    dplyr::mutate(across(everything(), as.character)) #%>% dplyr::mutate(ID =NA) %>% as.data.frame()
    
    
    
  } 
  
  
   
  
   
   
   #reformatting 

    Selected_Trees = melt(setDT(Selected_Trees), measure.vars = 1:ncol(Selected_Trees), variable.name = "Metric")  %>%
              dplyr::mutate(Metric = as.character(Metric) , Chain = substr(Metric,1,  str_locate(Metric, "\\." )-1) , Metric = substr(Metric, str_locate(Metric, "\\." ) +1   ,nchar(Metric)) )
    Selected_Trees = spread(Selected_Trees, Metric, value)
   
   
   #Keeping track of tree quality
   Selected_Trees$MaxConcordanceScore = max( chain_combis_tree_pair$MeanJaccard, na.rm = T)
   #Selected_Trees$n_lineage_discordance = -trunc(max_concordance_score/200) 
   
   
   
   # NEW filtering aberrant tree after concordant criteria if not met (remaining chain still best in terms of concordance) 
  Selected_Trees = Selected_Trees %>% dplyr::mutate(  gap_confidence_chain = ifelse( Chain == "Heavy" , heavy_confidence, 1-heavy_confidence )  )
  number_of_lineages_discordance = max( as.numeric(Selected_Trees$n_lineage_discordance), na.rm = T)
  write.table(Selected_Trees , file = file.path(GCTREE_FOLDER, "Selected_Trees_before_filtering.txt") , row.names = T)
  
  if( number_of_lineages_discordance > 0 ){ write.table("chain_discordance!" , file = file.path(GCTREE_FOLDER, "chain_discordance.txt") )  } 
  

  #max_concordance_score
   
   # Put an end to loop 
   available_chains <- NULL #putting an end to loop while preventing next if

}


if( length(available_chains)==1 ){
  Selected_Trees = Tree_Nk_Table %>%  dplyr::mutate(Final_Score = (gctree_contrib*abs_med_div(GcTree_ll) ) + mutpars_contrib*abs_med_div(Mut_Pars) ) %>%
                    dplyr::arrange(Final_Score) %>% dplyr::slice(1)
  print(paste0( "n selected trees  ", dim(Selected_Trees) ))
  histoplots = ggplot(data = NULL)
}


histoplots

