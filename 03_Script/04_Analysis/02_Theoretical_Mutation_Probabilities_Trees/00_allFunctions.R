######################################################
#####   This file contains custom functions to  ######
#####   compute theoretical mutation probas     ######
######################################################


print("importing functions")

#### General 
`%notin%` <- Negate(`%in%`)


###################################
# Extracting data from trees
###################################

# Extracts parental nodes from nk file. If unreal node => take node above


get_daughters_from_parents <- function(Newick_String, parent_nodes ){
  nkTree = str_replace_all(nkTree, ":\\d+", "") #removing ":number"
  phylo_table = as.data.frame( cbind(parent_nodes, matrix(nkTree, nrow = length(parent_nodes)) ) %>% set_colnames(c("Parents_Raw", "nk")) ) #build table 
  #Getting children nodes between parenthesis 
  phylo_table = phylo_table %>% dplyr::rowwise() %>% mutate(ChildNodes = gsub( paste0(")", Parents_Raw ,"[,);].*"), "\\1", nk), #getting everythink before parent in nk string
                                                            ChildNodes = gsub( "\\((?>[^()]|(?R))*\\)", "", ChildNodes, perl = T), #https://stackoverflow.com/questions/41749058/r-parse-nested-parentheses
                                                            ChildNodes = substr(ChildNodes,stri_locate_last_regex(ChildNodes, "[(]" )[1] +1  , 10000), 
                                                            DistoRoot =str_count(   gsub( paste0(".*", Parents_Raw), "", Newick_String  )   , "[)]"  )   )
                                                            #ChildNodes = sub("^.+[(]", "", ChildNodes) ) #https://stackoverflow.com/questions/50184888/use-stringr-in-r-to-find-the-remaining-string-after-last-substring
  phylo_table$ParentNode = phylo_table$Parents_Raw #to change whith real nodes afterwards ? 
  return (phylo_table)
}

position_in_tree <- function(Newick_String ){
  allnodes = unlist( str_split( str_replace_all(Newick_String, pattern = ":[0-9]{1,2}|[(;]|Ancestor", replacement = "" ), "[),]" )    )
  position_table = as.data.frame( allnodes[nchar(allnodes) > 0] ) #remove empty elements 
  position_table =  position_table %>% rename_at( 1, ~"NodeName" )   %>% dplyr::rowwise() %>% 
                   mutate( DistoRoot =str_count(   gsub( paste0(".*", NodeName), "", Newick_String  )   , "[)]"  )        ) # counting number of ) after seq in nk   
  return (position_table)
} 







###################################
# RUNNING BCR SIMULATION 
###################################


# Get sequences relations 
get_BCR_relations <- function(Tree_Nk, remove_virtual = T){ 
  if(remove_virtual){ #NEW
  	while( is.na( str_extract(Tree_Nk, "[(](.*?)[)][0-9]")) ==F  ){ #length(str_extract_all(Tree_Nk, "[(](.*?)[)][0-9]")) > 0
  		parenthesis = str_extract(Tree_Nk, "[(](.*?)[)][0-9]")
  		parenthesis =  gsub(".*\\(", "", parenthesis) #treating parenthesis to be recognized as regex 
  		replacement = gsub("\\)[0-9].*","",parenthesis) #removing virtual parent eg. )5
  		Tree_Nk = gsub(paste0("\\(", replacement, "\\)[0-9]{1,5}[:?][0-9]{1,3}"), replacement, Tree_Nk)
	}
  }
  Tree_ape = ape::read.tree(text = Tree_Nk )
  if("Ancestor" %notin% Tree_ape$node.label){ 
    #bug correction for linear trees https://stackoverflow.com/questions/71206191/tip-labels-not-appearing-when-reading-newick-tree-into-r-with-read-tree-and-re
    Tree_Nk_withbifurcation = str_replace(Tree_Nk, pattern = "[)]Ancestor", replacement = ",seqX:100)Ancestor") #adding artificial bifurcation
    Tree_ape = ape::read.tree(text = Tree_Nk_withbifurcation )
    Tree_ape = ape::drop.tip(Tree_ape, tip = c("seqX"), trim.internal = F, collapse.singles = F) #removing dummy sequence
  }
  Tree_Table = tidytree::as_tibble(Tree_ape) %>% as.data.frame() %>% dplyr::rename(SeqName = label) %>% 
               dplyr::mutate(SeqName = ifelse( grepl(pattern = "^[0-9]",SeqName ), paste0("X", SeqName), SeqName ) )
  Tree_Table$ParentName = Tree_Table[match(Tree_Table$parent, Tree_Table$node ), ]$SeqName
  return (Tree_Table)
}





#Checks whether sequence is productive ie not stop codons. !reading frame 
is_productive <- function(seqlist, stopcodons = c( "TAA", "TAG", "TGA")  ){
  stopcodons <- paste(stopcodons, collapse = "|")
  frameshift <- nchar(seqlist)%%3 # check that sequence is a multiple of 3 (for translation into aa)
  
  sequence_codons <- gsub("(.{3})", "\\1 ", seqlist)
  
  productive <- ifelse(  grepl(pattern = stopcodons, x = sequence_codons) , FALSE, TRUE)
  
  return (productive)
}

# Selects right mutation model and apply Shazam 
isMurineModel = "MOUSE|Mouse|mouse|Murine|MURINE|MUS|Mus|mus|usculus|USCULUS"
if( grepl(isMurineModel, DATASET) ){MH_RS5A = readRDS("./01_Reference/Substitution_Models/MH_RS5A.rds")  } # MH_RS5A.rds MH_RS5A.rds are not in the docker environment yet. 
applyMutationModel <- function(Sequence, chain, DATASET ){
  # Getting Right Mutation Model 
  isMurineModel = "MOUSE|Mouse|mouse|Murine|MURINE|MUS|Mus|mus|usculus|USCULUS"
  
  if(grepl(isMurineModel, DATASET)  ){
  	if(chain == "Heavy"){ mutModel = MH_RS5A #MH_RS5A
  	}else{mutModel = HKL_S5F } #MK_RS5NF }
  }else{
  	if(chain == "Heavy"){ mutModel = HH_S5F
  	}else{mutModel = HKL_S5F }
 }
  

  
  # Apply Shazam from Immcantation to mutate sequence + drops unproductive muts 
  Sequence_new = shmulateSeq( Sequence, numMutations = 1 , targetingModel = mutModel )
  z = 0 
  
  while(is_productive( Sequence_new) == F & z < 3 ){ # z < 3
  	z =z+1
  	Sequence_without_STOP = Sequence
  	
  	if(z > 1 ){ 
  		old_mutation = mutated 
  		substr(Sequence_without_STOP, old_mutation, old_mutation) <- "N" 
  	}
  	
  	mutated = as.vector(mapply(function(x, y) which(x != y), strsplit(Sequence_new, ""), strsplit(Sequence, "")))[1]
  	
  	substr(Sequence_without_STOP, mutated, mutated) <- "N"
  	Sequence_new2 = shmulateSeq(Sequence_without_STOP, numMutations=1, targetingModel = mutModel)
  	

  	substr(Sequence_new2, mutated, mutated) <- substr(Sequence, mutated, mutated)
  	
  	
  	if(z > 1 ){ substr(Sequence_new2, old_mutation, old_mutation) <- substr(Sequence, old_mutation, old_mutation) }  
  	
  	Sequence_new = Sequence_new2
  	
  	#for(element in mutated){
  	#	Sequence_new = Sequence
  	#	substr(Sequence_new, element, element) <- substr(Sequence_new2, element, element)
  	#	if(is_productive( Sequence_new)){
  	#	break
  	#	print(is_productive( Sequence_new))
  	#	}
  	#}
  	
  }
  return(Sequence_new)
}


# Get mutation effect (synonymous or non synonymous)



CheckDifferentAA <- function(Sequence, Sequence_new, CondonIndex){
  Sequence_proteic = translateDNA(Sequence)
  Sequence_new_proteic = translateDNA(Sequence_new)
  isNonSyn = T
  if (  is.na(CondonIndex) ){isNonSyn <- F} #here if the scanning through strplit doesn't return anything, the mutation is synonymous 
  if( substr(Sequence_proteic, CondonIndex, CondonIndex  ) == substr(Sequence_new_proteic, CondonIndex, CondonIndex)    ){isNonSyn <- F} 
  isNonSyn = 1
  return(isNonSyn)
}




###################################
#PATTERN ANALYSIS
###################################


# Gets all or nucleotides/codons of interest (eg after start)  involved in pattern 
vectorize_regex_idx <- function(seq, pattern, patternlength, convert_to_nuc = F, after_start=0){
  loc <- str_locate_all(seq, pattern)[[1]] #res <- integer(100)
  if (convert_to_nuc == F){idx <- as.numeric(sapply(loc[,1], function(x) seq(x+after_start, x+after_start+patternlength  ) )) 
  }else{ idx <- as.numeric(sapply(loc[,1], function(x) seq(  ( x*3 -2 )  ,  (x+patternlength)*3   ) ))       }
  return(idx)
}



#-------------------------------------------------------------------------------
# Determines which patterns inherited from NCA are mutated/disrupted 
PatternsFromParent <- function(currentDF, AncestralColumn,  PATTERN_TABLE ){
  # Get mutations to be compared to NCA for pattern analysis
  currentDF$CodonRegex =  paste0( "[( )]" , currentDF$Nuc_Index_Absolute, "[,)]")
  
  ############ Controls. Only cytidines targeted by AID are taken ############
  ## Extracting patterns 
  currentDF$RandomSites = as.character(  map( currentDF[,which(names(currentDF) == AncestralColumn )] , vectorize_regex_idx, pattern= "C", patternlength = 0, convert_to_nuc = F, after_start = 0) )
  currentDF$WRC =    as.character(  map( currentDF[,which(names(currentDF) == AncestralColumn )] , vectorize_regex_idx, pattern= "[AT]{1}[AG]{1}C", patternlength = 0, convert_to_nuc = F, after_start = 2) )
  currentDF$SYC =    as.character( map( currentDF[,which(names(currentDF) == AncestralColumn )] , vectorize_regex_idx, pattern=  "[GC]{1}[TC]{1}C", patternlength = 0, convert_to_nuc = F, after_start = 2) )
  
  # Matching positions with motifs initially in parent node or NCA  
  currentDF$isWRC <- mapply(grepl, pattern= currentDF$CodonRegex , x=currentDF$WRC) # Why is long and not others ? 
  currentDF$isSYC <- mapply(grepl, pattern= currentDF$CodonRegex , x=currentDF$SYC)
  currentDF$isRandom <- mapply(grepl, pattern= currentDF$CodonRegex , x=currentDF$RandomSites)
  
  
  ############ Input Patterns ########################
  for(pattern in unique(PATTERN_TABLE$Pattern_Name) ){
    print(pattern)
    
    # Check whether input patterns are in NCA 
    current_pattern = PATTERN_TABLE %>% dplyr::filter(Pattern_Name == pattern)
    if(current_pattern$Regex =="*"){current_pattern$Regex <- "\\*"} #"\\ not supported by read table. corresponds to stop codons
    
    if( grepl("Prot|prot|AA|aa|mino", current_pattern$Seq_Type)  ){currentDF$newColumn= as.character(  map( currentDF$ProtSeq , vectorize_regex_idx, current_pattern$Regex,  current_pattern$Pattern_Length - 1,  T)   ) 
    }else{ currentDF$newColumn= as.character( map( currentDF[,which(names(currentDF) == AncestralColumn )], vectorize_regex_idx, current_pattern$Regex,  current_pattern$Pattern_Length - 1,  F)  ) }
     
    # Matching positions with motifs initially in NCA 
    currentDF  = currentDF %>% dplyr::rowwise() %>% dplyr::mutate(isPattern = grepl(pattern = CodonRegex, x = newColumn )   )

    #Renaming temporary columns with right pattern 
    currentDF = currentDF %>% dplyr::mutate(newColumn = ifelse( grepl("numeric", newColumn), NA, newColumn )  )  %>% rename_at( c("newColumn") , list( ~paste0(  c(pattern) ) )    ) # getting rid of numeric 0 and renaming new column to pattern 
    currentDF = currentDF %>% rename_at( c("isPattern") , list( ~paste0( "is" , c(pattern) ) ) )
    
  } #END of pattern loop
  
  # Summarizes in which pattern nucleotide is located
  ## Source code for this part: https://stackoverflow.com/questions/63249304/create-new-column-that-takes-column-name-based-on-value
  currentDF$Site <- paste(apply(currentDF[, grep("^is(?!.*WRC|NNC|SYC|NonSyn|Random)",  names(currentDF),  value = TRUE, perl = TRUE)]     ,1  ,   function(x) names(x)[which(x==T)] )) #bug is only NA
  currentDF = currentDF  %>%  dplyr::mutate( AIDSpots = case_when(( isRandom == T & isWRC ==F & isSYC ==F) ~ "NNC",
                                                      ( isWRC == T  ) ~ "WRC",
                                                      ( isSYC == T & isWRC == F) ~ "SYC", 
                                                      TRUE ~ "other") , Site  = ifelse( grepl("is", Site), Site,  AIDSpots  ) ) 
  
  return(currentDF)
} #END of the function



# Determines which patterns inherited from NCA are mutated/disrupted => only AID spots now => other patterns after NCA 
PatternsFromParent <- function(currentDF, AncestralColumn,  PATTERN_TABLE ){
  # Get mutations to be compared to NCA for pattern analysis
  currentDF$CodonRegex =  paste0( "[( )]" , currentDF$Nuc_Index_Absolute, "[,)]")
  
  ############ Controls. Only cytidines targeted by AID are taken ############
  ## Extracting patterns 
  currentDF$RandomSites = as.character(  map( currentDF[,which(names(currentDF) == AncestralColumn )] , vectorize_regex_idx, pattern= "C", patternlength = 0, convert_to_nuc = F, after_start = 0) )
  currentDF$WRC =    as.character(  map( currentDF[,which(names(currentDF) == AncestralColumn )] , vectorize_regex_idx, pattern= "[AT]{1}[AG]{1}C", patternlength = 0, convert_to_nuc = F, after_start = 2) )
  currentDF$SYC =    as.character( map( currentDF[,which(names(currentDF) == AncestralColumn )] , vectorize_regex_idx, pattern=  "[GC]{1}[TC]{1}C", patternlength = 0, convert_to_nuc = F, after_start = 2) )
  
  # Matching positions with motifs initially in parent node or NCA  
  currentDF$isWRC <- mapply(grepl, pattern= currentDF$CodonRegex , x=currentDF$WRC) # Why is long and not others ? 
  currentDF$isSYC <- mapply(grepl, pattern= currentDF$CodonRegex , x=currentDF$SYC)
  currentDF$isRandom <- mapply(grepl, pattern= currentDF$CodonRegex , x=currentDF$RandomSites)
  

  currentDF = currentDF  %>%  dplyr::mutate( Site = case_when(( isRandom == T & isWRC ==F & isSYC ==F) ~ "NNC",
                                                      ( isWRC == T  ) ~ "WRC",
                                                      ( isSYC == T & isWRC == F) ~ "SYC", 
                                                      TRUE ~ "other")  ) 
  
  return(currentDF)
} #END of the function


#-------------------------------------------------------------------------------

####################################################
# Determining mutations from experimental data
####################################################



get_indexes_and_effects_vector <- function( data_df,  seq_col, uca_col, read_frame){ #
  
  # get the cells that correspond to the requirement (root or leaf)

  line_index_set = seq(1:nrow(data_df) )

  # for each cell, look at the parent diff indexes and get their effect
  all_index_set = vector()
  effect_set = vector()
  seqID_set = vector()
  for( line_index in line_index_set){
    
    # get the BCR sequence of the cell
    cell_sequence = toupper(data_df[ line_index, seq_col])
    uca_sequence =  substr( toupper(data_df[ line_index, uca_col])  , 1 , nchar(cell_sequence)  ) #only compare seq of same length
    
    # get the parent diff indexes of the cell
    index_set = as.numeric( mapply(function(x, y) which(x != y), strsplit(cell_sequence, ""), strsplit(uca_sequence, "")  ) ) 

    

    
    if ( is.na(index_set[1]) == FALSE){
      for( index in index_set[which(index_set!=nchar(uca_sequence))]   ){  # bug when index is last nuc
        # get the index of the SNV codon start. In the following computation :
        #   add to the first codon start index to get the SNV codon start index
        codon_start_index = read_frame + 3 * floor( (index - read_frame) / 3)
        
        # get the codon sequences in cell BCR and UCA BCR
        cell_codon = substr( cell_sequence, start = codon_start_index, stop = codon_start_index + 2)

        uca_codon = substr( uca_sequence, start = codon_start_index, stop = codon_start_index + 2)

        
        
        cell_aa = GENETIC_CODE[[ cell_codon]]
        uca_aa = GENETIC_CODE[[ uca_codon]]
        
        if( is.na( cell_aa)){
          cat("\n Error : abnormal codon on cell BCR : ", cell_codon, "=>", cell_aa)
        }else if( is.na( uca_aa)){
          cat("\n Error : abnormal codon on UCA BCR : ", uca_codon, "=>", uca_aa)
        }else{
          effect = "non-syn"
          if( cell_aa == uca_aa){
            effect = "syn"
          }
          all_index_set = append( all_index_set, index)
          effect_set = append( effect_set, effect)
          seqID_set = append( seqID_set,   data_df[ line_index, "sequence_id"]  )
        }
      } 
    }
    
  }
  
  return( list( indexes = all_index_set, effects = effect_set, sequence_id = seqID_set ))  
}


print("functions imported")



