# Preselection 
# Dropping trees not respecting CSR rules or with unproductive trees.
# Best Mut Pars tree is chosen within each topology 
# -----------------------

## @knitr Productivity_CSR_topology





options(digits=13)


Tree_Nk_Table = read.csv( text =  paste( c( "Chain", "Tree_Nk", "Tree_Index", "GcTree_ll", "Mut_Pars", "Iso_Pars", "Pars_Penalty", "all_aberrations", "Productive", "toKeep" ) , collapse = ",")  )


all_tree_scores = list()
c = 0


for(current_chain in available_chains){
  cat(' \n \n')
  cat('##',current_chain,' \n')
  
  # ----------Fetching data----------- 
  
  
  tree_scores= read.csv(file = file.path(GCTREE_FOLDER, current_chain, "treestat_custom.txt"), header = T, sep = "\t") %>% distinct() #in case
  c =c+1
  all_tree_scores[[c]] <- tree_scores %>% dplyr::mutate(Chain = current_chain)
  
  forest_size = nrow(tree_scores)
  
  # ----------Prefiltering based isotype + thining forest ----------- 
  
    # Keeping least aberrant trees in terms of mutation order 
  max_mut_aberration = min(tree_scores[which(tree_scores$Pars_Penalty == 0  ) ,]$all_aberrations, na.rm = T   )
  tree_scores = tree_scores  %>% dplyr::filter( all_aberrations <=  max_mut_aberration    ) 	 
  
  
  print( paste0("Trees with aberrant isotype switch: " ,nrow(tree_scores[which(tree_scores$Iso_Pars >0),]), "/", forest_size  ) )
  tree_scores = tree_scores %>% dplyr::filter(Iso_Pars == min(Iso_Pars)   ) 
  

  


 # ----------Productivity of parents -----------  
  unproductive_trees = c()
  #aberrant_initialization = c() 
  for (cur_idx in tree_scores$Tree_Index ){
    cur_tree = readChar(file.path(GCTREE_FOLDER, current_chain, paste0("tree_",cur_idx,".nk")) , nchar = 1e6 )
     
    
    parent_nodes = as.vector( str_extract_all(cur_tree,  "(?<=[)])(.*?)(?=[:])", simplify = TRUE ) ) #between ) and :
    
    
    
      if(   file.info(   file.path(GCTREE_FOLDER, current_chain,  paste0("tree_",cur_idx,".fasta") )   )$size ==0   ){
  	#no consistent tree from previous step 
  	q(save = "no", status = 0, runLast = FALSE)
  	}
    
       
    cur_fasta = seqinr::read.fasta(file = file.path(GCTREE_FOLDER, current_chain, paste0("tree_",cur_idx,".fasta")),
                                   as.string = T, set.attributes = F, forceDNAtolower = F) 
    cur_fasta = data.frame(SeqName = names(cur_fasta) , Sequence = paste(cur_fasta) )  #%>% 
    
    
    # Check that sequences that are 1 mut dist from ancestor are directly linked to it (not necessarily ensured by GcTree) 

    # Productivity 
    cur_fasta_filtered  = cur_fasta %>% dplyr::filter( SeqName %in% parent_nodes[parent_nodes !="Ancestor"] ) %>% 
                dplyr::mutate(isProductive  = is_productive(Sequence) )
    
    any_unproductive = F %in% cur_fasta_filtered$isProductive
    if(any_unproductive ){unproductive_trees = c(unproductive_trees, cur_idx )} 
    
    
    
    
    
    # Record all tree data 
    Tree_Nk_Table[nrow(Tree_Nk_Table)+1,] = c(current_chain,  cur_tree, 
                                          as.numeric(tree_scores[which(tree_scores$Tree_Index == cur_idx),]), !any_unproductive, "empty" )
    
  }
  
  ## ---------- confidence level for chains (next script) -----------
  ancestral_seq = as.character( cur_fasta[which(cur_fasta$SeqName == "Ancestor"),]$Sequence   )
  seq_distances_df =   cur_fasta %>% dplyr::filter(  grepl("[Ss]eq|ncestor", SeqName ) ) 
  seq_distances_df$Nmut = stringdist(ancestral_seq, seq_distances_df$Sequence, "osa")
  mutation_gaps = stringdist::stringdistmatrix(seq_distances_df$Sequence, seq_distances_df$Sequence, method = "osa")
  if (current_chain == "Heavy"){quantileH = as.numeric(quantile(mutation_gaps, 0.75))
  } else {quantileL =as.numeric(quantile(mutation_gaps, 0.75))}
  
  
  # ----------1 tree per topology -----------
  
  tree_scores = as.data.frame(tree_scores) %>% filter(!Tree_Index %in% unproductive_trees ) %>% 
  		#dplyr::filter( !Tree_Index   %in% aberrant_initialization ) %>% 
  		group_by(GcTree_ll) %>%
                dplyr::arrange(Mut_Pars) %>% dplyr::slice(1) %>% ungroup() #%>% top_n(n = 1, Tree_Index)
  print(dim(tree_scores))

   


  
  
  cat(' \n \n')
  
 

  
  # ----------Selecting non aberrant trees -----------
  
  
  ##Correspondence between contig Id and Gctree seq names 
  chain_seq_name =  read.csv(file.path(GCTREE_FOLDER, current_chain, "idmap.txt"), header = F) %>% separate_rows(V2, sep=":")
  colnames(chain_seq_name) =  c( paste0(tolower(current_chain), "_SeqName" ), paste0(tolower(current_chain), "_sequence_id" )  )
  airr_table = dplyr::left_join(airr_table, chain_seq_name ) #%>% na.omit()
  
  #Isotype NEW-----
  if( file.info(   file.path(GCTREE_FOLDER, current_chain, 'isotypemap.txt')   )$size > 0 ){
  	isotypes_df = read.table(  file = file.path(GCTREE_FOLDER, current_chain, 'isotypemap.txt'), header = F, sep = "," )  %>% dplyr::rename(sequence_id = V1, Isotype = V2)
  }else{
  	isotypes_df = as.data.frame( chain_seq_name[,2] ) %>% dplyr::mutate(Isotype = NA) 
  	colnames(chain_seq_name) = c("sequence_id", "Isotype")
  } 
  
  

    colnames(isotypes_df) =  c( paste0(tolower(current_chain), "_sequence_id" ), paste0(tolower(current_chain), "_Isotype" )  )

  airr_table = dplyr::left_join(airr_table, isotypes_df )
  
  
  # Saving heavy and light tree indexes to run concordance criterium 
  if(current_chain == "Light"){ Light_Indexes = tree_scores$Tree_Index
  }else if(current_chain == "Heavy" ){Heavy_Indexes = tree_scores$Tree_Index}  

  Tree_Nk_Table = Tree_Nk_Table %>% mutate(toKeep = ifelse(Chain == current_chain, !(Tree_Index %notin% tree_scores$Tree_Index & Chain == current_chain ),  toKeep  ) )

}


Tree_Nk_Table = Tree_Nk_Table %>% mutate_at(vars(c("Tree_Index", "GcTree_ll" , "Mut_Pars", "Iso_Pars", "Pars_Penalty", "all_aberrations" )),as.numeric)
Tree_Nk_Table_chains = Tree_Nk_Table %>% group_split(Chain, .keep = FALSE)
 print(names(Tree_Nk_Table))
 
 write.table(Tree_Nk_Table, file = file.path(GCTREE_FOLDER, "Tree_Nk_Table.txt") )

