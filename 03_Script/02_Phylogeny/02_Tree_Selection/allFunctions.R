## @knitr load_functions

# #####################################################################################
# #####################################################################################
# #####################################################################################
########             Generate a datatable summarizing values                   ######## 
########    For environments (parameters), all values/variables are shown      ######## 
# #####################################################################################
# #####################################################################################
# #####################################################################################




############################################################################################
# Tool functions  
############################################################################################
`%notin%` <- Negate(`%in%`)



abs_med_div <- function(vector ){
  return(abs(vector/max(vector, na.rm = T)  )) 
}

Mode <- function(x, na.rm = FALSE) {
  # it takes two areguments x, and na.rm (https://www.educative.io/answers/what-is-the-mode-method-in-r)
  if(na.rm){ #if na.rm is false it means no need to remove NA values
    x = x[!is.na(x)]
  }
  valx <- unique(x)
  return(valx[which.max(tabulate(match(x, valx)))])
}


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


############################################################################################
# Testing BCRseq productivity   
############################################################################################


is_productive <- function(seqlist, stopcodons = c( "TAA", "TAG", "TGA")  ){
  stopcodons <- paste(stopcodons, collapse = "|")
  frameshift <- nchar(seqlist)%%3 # UCA sont des multiples de 3. A vérifier avec echidna 
  
  #sequence_codons <- stri_sub(seqlist, seq(1, stri_length(seqlist),by=3), length=3)
  sequence_codons <- gsub("(.{3})", "\\1 ", seqlist)
  
  #productive <- ifelse(frameshift != 0 |  grepl(pattern = stopcodons, x = seqlist) , FALSE, TRUE)
  productive <- ifelse(  grepl(pattern = stopcodons, x = sequence_codons) , FALSE, TRUE)
  
  return (productive)
}


############################################################################################
# Determines a sequence lineage based on newick string. Virtual sequences could be included or not  
############################################################################################


#bug correction for linear trees which generate bugs with ape package:

getParents <- function(Tree_Nk , SeqOfInterest, real = T ){
  Tree_ape = ape::read.tree(text = Tree_Nk )
  if("Ancestor" %notin% Tree_ape$node.label){ #bug correction for linear trees https://stackoverflow.com/questions/71206191/tip-labels-not-appearing-when-reading-newick-tree-into-r-with-read-tree-and-re
    Tree_Nk_withbifurcation = str_replace(Tree_Nk, pattern = "[)]Ancestor", replacement = ",seqX:100)Ancestor") #adding artificial bifurcation
    Tree_ape = ape::read.tree(text = Tree_Nk_withbifurcation )
    Tree_ape = ape::drop.tip(Tree_ape, tip = c("seqX"), trim.internal = F, collapse.singles = F) #removing dummy sequence
  }
  Tree_Table = tidytree::as_tibble(Tree_ape) %>% as.data.frame() 
  Root = Tree_Table[which(Tree_Table$label == "Ancestor"),]$node
  End = Tree_Table[which(Tree_Table$label == SeqOfInterest),]$node
  Node_Path = ape::nodepath(Tree_ape, Root, End) 
  
  Node_Path = Tree_Table[match(Node_Path[-length(Node_Path)], Tree_Table$node ), ]$label #removing last element and matching with label
  if(real == T){Node_Path = Node_Path[grepl(pattern = "Ancestor|seq", Node_Path)] }
  Node_Path = c(Node_Path, SeqOfInterest) #NEW include current node into lineage 
  Node_Path = Node_Path[Node_Path != ""] #IgPhymL
  return( paste(Node_Path, collapse = ",") )
} 


getDaughters <- function(Tree_Nk , SeqOfInterest, real = F ){
    Tree_ape = ape::read.tree(text = Tree_Nk )
    if("Ancestor" %notin% Tree_ape$node.label){ 
        Tree_Nk_withbifurcation = str_replace(Tree_Nk, pattern = "[)]Ancestor", replacement = ",seqX:100)Ancestor") #adding artificial bifurcation
        Tree_ape = ape::read.tree(text = Tree_Nk_withbifurcation )
        Tree_ape = ape::drop.tip(Tree_ape, tip = c("seqX"), trim.internal = F, collapse.singles = F) #removing dummy sequence
    }
    
    Tree_Table = tidytree::as_tibble(Tree_ape) %>% as.data.frame()
    node_number = Tree_Table[which(Tree_Table$label == SeqOfInterest),]$node
    daugthers = Tree_Table[which(Tree_Table$parent == node_number),]$label
    
    if(real == T){daugthers = daugthers[grepl(pattern = "Ancestor|seq", daugthers)] }
    
    
    #return( paste(daugthers, collapse = ",")  )
    return( daugthers  )
}




##########################################################################
#This functions checks whether isotype order is respected along lineages 
##########################################################################






Check_Isotype_Order <- function(Isotypes, ParentIsotypes, CSR_order){
  #Isotypes = "IgG2,IgM"
  #ParentIsotypes = "IgG2"
  
  Isotypes = unlist(str_split(Isotypes, pattern = ',') )
  ParentIsotypes = unlist( str_split(ParentIsotypes, pattern = ',') )
  
  Isotypes = Isotypes[order(match(Isotypes, CSR_order  ))]
  
  #subseting before 1st occurence of child
  CSR_order_before_child = CSR_order[1:min(which(CSR_order == Isotypes[1] ), length(CSR_order) , na.rm = T)]
  ParentIsotypes[which(ParentIsotypes %in% CSR_order_before_child )]
  
  ParentIsotypes = ParentIsotypes[order(match(ParentIsotypes, CSR_order  ))] 
  
  uniqueIso = unique(  c( ParentIsotypes,  Isotypes ) )
  
  
  if(  identical(  uniqueIso,     uniqueIso[order(match(uniqueIso, CSR_order  ))]       )       | is.na(ParentIsotypes )   |   is.na(Isotypes )      ){
     isotype_disconcordance = 0
  }else{ isotype_disconcordance = 1  } 
  
  return( isotype_disconcordance ) 
  
}



############################################################################################
# This function penalizes pairs of trees that are not concordant + computes degree of overlap of lineages   
############################################################################################









compare_lineage <-  function(HeavyParents, LightParents, Hindex, Lindex, parent_df, lookup_table){ 
  
  max_penalty = -200
  
  # Selecting tree pair
  parent_df = parent_df %>% dplyr::filter(Heavy.Tree_Index ==Hindex & Light.Tree_Index == Lindex ) 

  
  # Getting and formating lineage
  HeavyParents = unlist(str_split(tolower(HeavyParents), pattern = ',') )#to avoid ambiguities between annotations 
  LightParents = unlist( str_split(tolower(LightParents), pattern = ',') )
  HeavyParents = rev(HeavyParents) #walk trought lineages from terminal leaf to root
  LightParents = rev(LightParents)
  HeavyParents[HeavyParents == "ancestor"] <- "heavyancestor"
  LightParents[LightParents == "ancestor"] <- "lightancestor"
  
  HeavyParents_all_combis = c()
  LightParents_all_combis = c()
  
  unpaired_penalties = 0

    for (current_rank in 1:max( length(LightParents) , length(HeavyParents)   ) ){
      
      # Heavy Tree concordance -----------------------------------
      
      if( current_rank < length(HeavyParents) ){ #ancestor case is trivial + heavy and light lineages not necessarily same length 
        
        print(HeavyParents[ current_rank  ])
        
        #Getting all possible combinations involving current heavy and light sequence respectively 
        HeavyCombi = unique(tolower(lookup_table[which( tolower(lookup_table$HeavySeq) == HeavyParents[ current_rank  ] ),]$BCRcombi))
        HeavyCombi = HeavyCombi[which( grepl(pattern = "\\.", HeavyCombi)  )] # only paired chains are taken into account 
        
        # Getting parents, sisters and daughters 
        #heavy_tree_parent = HeavyParents[ current_rank  + 1 ]
        #heavy_tree_sisters = parent_df[which(parent_df$Heavy.LastParent %in% heavy_tree_parent  ),]$HeavySeq
        heavy_tree_daugthers = parent_df[which(parent_df$Heavy.LastParent %in%  HeavyParents[ current_rank   ]  ),]$LightSeq
        
        
        # Checking that the sisters in one tree are linked in other tree (same parent, grandparent or daughter of other sister). 
        light_tree_parent = parent_df[which(parent_df$BCRcombi %in% HeavyCombi ), ]$Light.LastParent
        #light_tree_grandparent = parent_df[which(parent_df$LightSeq %in% light_tree_parent ), ]$Light.LastParent
        #light_tree_sisters = parent_df[which(parent_df$Light.LastParent %in% light_tree_parent  ),]$LightSeq
        
        
        # if  heavy parent is unpaired, not taken into account for penalty, but trees where the descendence is clustered downstream of this unpaired sequence are preferred 
        #light_tree_parent_combi =   parent_df[which(parent_df$LightSeq %in% light_tree_parent  ),]$BCRcombi
        #if( length(light_tree_parent_combi) == 0 & all(parent_df[which(parent_df$BCRcombi %in% HeavyCombi ), ]$Light.GrandParent  %in%  c(".", "ancestor") )   ){
        #  unpaired_penalty = -1
        #  unpaired_penalties = unpaired_penalties + unpaired_penalty
        #}else{unpaired_penalty = 0 }
        
        # if  heavy parent is unpaired, not taken into account for penalty, but trees where the descendence is clustered downstream of this unpaired sequence are preferred 
    	heavy_tree_daugthers_light_parent_unpaired = unique(parent_df[which(parent_df$LightSeq %in%  heavy_tree_daugthers  ),]$Light.LastParent)
    	heavy_tree_daugthers_light_parent_unpaired = heavy_tree_daugthers_light_parent_unpaired[which(heavy_tree_daugthers_light_parent_unpaired %notin% c(".", "lightancestor", parent_df$LightSeq ))] 
    	if( length( heavy_tree_daugthers_light_parent_unpaired ) > 0   ){
    		 unpaired_penalty = -1
     		 unpaired_penalties = unpaired_penalties + unpaired_penalty
      	
     		 heavy_tree_daugthers = parent_df[which(parent_df$Heavy.LastParent %in%  HeavyParents[ current_rank   ] & 
                                               parent_df$Light.LastParent %notin% heavy_tree_daugthers_light_parent_unpaired & 
                                               parent_df$Heavy.GrandParent %notin% heavy_tree_daugthers_light_parent_unpaired ),]$LightSeq 
      
    	}else{unpaired_penalty = 0 }
    	
    	#if heavy tree parent is unpaired, replace by grand parent 
    	for( current_parent in light_tree_parent ){
   	    if(length( parent_df[which(parent_df$LightSeq %in% current_parent ), ]$Light.LastParent)  ==0   ){
    	    	light_tree_parent[light_tree_parent == current_parent] <-  unique(parent_df[which(parent_df$Light.LastParent %in% current_parent ), ]$Light.GrandParent) 
    	    	unpaired_penalties = unpaired_penalties -1 
      	    }
    	}
    	light_tree_sisters = parent_df[which(parent_df$Light.LastParent %in% light_tree_parent  ),]$LightSeq
        
        
        
        if(any(light_tree_parent %notin% c(".") )   &  length(light_tree_parent) > 0){ # trivial as all lineages descending from ancestor    #
          
          current_light_chains = str_extract(HeavyCombi, "lightseq[0-9]{1,5}")
          paired_light_chains  = current_light_chains
          light_tree_sisters = light_tree_sisters[which(light_tree_sisters %notin% paired_light_chains )]
          
          
          
          light_tree_direct_descendance = c(paired_light_chains)
          light_tree_sister_descendance = c(light_tree_sisters)
          
          while(length(paired_light_chains) > 0 | length(light_tree_sisters)  > 0 ){
            
            # direct descendance
            paired_light_chains = parent_df[which(parent_df$Light.LastParent %in% paired_light_chains  ),]$LightSeq
            light_tree_direct_descendance = c(light_tree_direct_descendance, paired_light_chains)
            
            # sister descendance
            light_tree_sisters = parent_df[which(parent_df$Light.LastParent %in% light_tree_sisters  ),]$LightSeq
            light_tree_sister_descendance = c(light_tree_sister_descendance, light_tree_sisters)
            
            
          }
          if(  any(heavy_tree_daugthers %notin% c(light_tree_direct_descendance, light_tree_sister_descendance) )  ){return(max_penalty)} # & unpaired_penalty == 0 
          
          
          
          # Detecting gaps in lineages ---------------------------------
          # Order will be checked later 
          
          if( (length(intersect(current_light_chains, light_tree_direct_descendance) )== 0   &  length(intersect(current_light_chains, light_tree_parent) )== 0  ) &
              length(intersect(light_tree_parent, light_tree_direct_descendance)) > 0    ){ return(max_penalty) }
          
          
          
        }
        
        HeavyCombi = HeavyCombi[order(match( sub("\\..*", "", HeavyCombi)  , LightParents   ))]   
        HeavyParents_all_combis = c(HeavyParents_all_combis, HeavyCombi )
        
      }else if( current_rank == length(HeavyParents) ){ # ancestor case 
          print(HeavyParents[ current_rank  ])
          HeavyCombi = unique(tolower(lookup_table[which( tolower(lookup_table$HeavySeq) == HeavyParents[ current_rank  ] ),]$BCRcombi))
          HeavyParents_all_combis = c(HeavyParents_all_combis, HeavyCombi[order(match( sub("\\..*", "", HeavyCombi)  , LightParents   ))] )
        
          # are the combis directly related ie direct daughter of ancestor ie  1 parent ? --------------
          current_light_chains =  sub(".*\\.", "", HeavyCombi)
          light_tree_parent = unique(parent_df[which(parent_df$BCRcombi %in% HeavyCombi ), ]$Light.LastParent)
        
          #NB: parent has to be paired to be informative, otherwise considered as virtual intermediary 
          # checking that parents are in combi 
          for( current_parent in light_tree_parent ){
              if(length( parent_df[which(parent_df$LightSeq %in% current_parent ), ]$Light.LastParent)  ==0   ){
                #changing parent by paired grand parent 
                light_tree_parent[light_tree_parent == current_parent] <-  unique(parent_df[which(parent_df$Light.LastParent %in% current_parent ), ]$Light.GrandParent) 
              }
          }
          
          # penalty if remaining combis are not sisters or parent from one another 
          if(length(   light_tree_parent[which(light_tree_parent %notin% c(current_light_chains, "lightancestor" , ".") )]     ) > 0 ){   return(max_penalty)   }
          
        }
    
        
      
      
      
      # Light Tree concordance -----------------------------------
      
      if( current_rank < length(LightParents) ){ 
        print(LightParents[ current_rank  ])
        #Getting all possible combinations involving current heavy and light sequence respectively 
        LightCombi = unique(tolower(lookup_table[which( tolower(lookup_table$LightSeq) == LightParents[ current_rank  ] ),]$BCRcombi))
        
        # Getting parents, sisters and daughters 
        light_tree_daugthers = parent_df[which(parent_df$Light.LastParent %in%  LightParents[ current_rank   ]  ),]$HeavySeq
    
        # Checking that the sisters in one tree are linked in other tree (same parent, grandparent or daughter of other sister). 
        heavy_tree_parent = parent_df[which(parent_df$BCRcombi %in% LightCombi ), ]$Heavy.LastParent
        #heavy_tree_sisters = parent_df[which(parent_df$Heavy.LastParent %in% heavy_tree_parent  ),]$HeavySeq
        
        
        # if  heavy parent is unpaired, not taken into account for penalty, but trees where the descendence is clustered downstream of this unpaired sequence are preferred 
        #heavy_tree_parent_combi =   parent_df[which(parent_df$HeavySeq %in% heavy_tree_parent  ),]$BCRcombi
        #if( length(heavy_tree_parent_combi) == 0 & all(parent_df[which(parent_df$BCRcombi %in% LightCombi ), ]$Heavy.GrandParent  %in%  c(".", "ancestor") )   ){
        #  unpaired_penalty = -1
        #  unpaired_penalties = unpaired_penalties + unpaired_penalty
        #}else{unpaired_penalty = 0 }
        
        # if  heavy parent is unpaired, not taken into account for penalty, but trees where the descendence is clustered downstream of this unpaired sequence are preferred 
    	light_tree_daugthers_heavy_parent_unpaired = unique(parent_df[which(parent_df$HeavySeq %in%  light_tree_daugthers  ),]$Heavy.LastParent)
    	light_tree_daugthers_heavy_parent_unpaired = light_tree_daugthers_heavy_parent_unpaired[which(light_tree_daugthers_heavy_parent_unpaired %notin% c(".", "heavyancestor", parent_df$HeavySeq ))] 
    	if( length( light_tree_daugthers_heavy_parent_unpaired ) > 0   ){
    		 unpaired_penalty = -1
     		 unpaired_penalties = unpaired_penalties + unpaired_penalty
      	
     		 light_tree_daugthers = parent_df[which(parent_df$Light.LastParent %in%  LightParents[ current_rank   ] & 
                                               parent_df$Heavy.LastParent %notin% light_tree_daugthers_heavy_parent_unpaired & 
                                               parent_df$Heavy.GrandParent %notin% light_tree_daugthers_heavy_parent_unpaired ),]$HeavySeq 
      
    	}else{unpaired_penalty = 0 }
    	
    	
    	#if heavy tree parent is unpaired, replace by grand parent 
    	for( current_parent in heavy_tree_parent ){
   	    if(length( parent_df[which(parent_df$HeavySeq %in% current_parent ), ]$Heavy.LastParent)  ==0   ){
    	    	heavy_tree_parent[heavy_tree_parent == current_parent] <-  unique(parent_df[which(parent_df$Heavy.LastParent %in% current_parent ), ]$Heavy.GrandParent) 
    	    	unpaired_penalties = unpaired_penalties -1 
      	    }
    	}
    	heavy_tree_sisters = parent_df[which(parent_df$Heavy.LastParent %in% heavy_tree_parent  ),]$HeavySeq
        
        
        if(any(heavy_tree_parent %notin% c(".") )   &  length(heavy_tree_parent) > 0){ # trivial as all lineages descending from ancestor    #
          
          current_heavy_chains = str_extract(LightCombi, "heavyseq[0-9]{1,5}")
          paired_heavy_chains  = current_heavy_chains
          heavy_tree_sisters = heavy_tree_sisters[which(heavy_tree_sisters %notin% paired_heavy_chains )]
          
          
          
          heavy_tree_direct_descendance = c(paired_heavy_chains)
          heavy_tree_sister_descendance = c(heavy_tree_sisters)
          
          while(length(paired_heavy_chains) > 0 | length(heavy_tree_sisters)  > 0 ){
            
            # direct descendance
            paired_heavy_chains = parent_df[which(parent_df$Heavy.LastParent %in% paired_heavy_chains  ),]$HeavySeq
            heavy_tree_direct_descendance = c(heavy_tree_direct_descendance, paired_heavy_chains)
            
            # sister descendance
            heavy_tree_sisters = parent_df[which(parent_df$Heavy.LastParent %in% heavy_tree_sisters  ),]$HeavySeq
            heavy_tree_sister_descendance = c(heavy_tree_sister_descendance, heavy_tree_sisters)
            
            
          }
          if(  any(light_tree_daugthers %notin% c(heavy_tree_direct_descendance, heavy_tree_sister_descendance) )    ){return(max_penalty)} #& unpaired_penalty == 0
          
          
          
          # Detecting gaps in lineages ---------------------------------
          # Order will be checked later 
          
          if( (length(intersect(current_heavy_chains, heavy_tree_direct_descendance) )== 0   &  length(intersect(current_heavy_chains, heavy_tree_parent) )== 0  ) &
              length(intersect(heavy_tree_parent, heavy_tree_direct_descendance)) > 0    ){ return(max_penalty) }
          
          
          
        }
        
        LightCombi = LightCombi[order(match( sub("\\..*", "", LightCombi)  , HeavyParents   ))]   
        LightParents_all_combis = c(LightParents_all_combis, LightCombi )
        
      }else if( current_rank == length(LightParents) ){ # ancestor case 
        print(LightParents[ current_rank  ])
        LightCombi = unique(tolower(lookup_table[which( tolower(lookup_table$LightSeq) == LightParents[ current_rank  ] ),]$BCRcombi))
        LightParents_all_combis = c(LightParents_all_combis, LightCombi[order(match( sub("\\..*", "", LightCombi)  , HeavyParents   ))] )
        
        # are the combis directly related ie direct daughter of ancestor ie  1 parent ? --------------
        current_heavy_chains =  sub("\\..*", "", LightCombi)
        heavy_tree_parent = unique(parent_df[which(parent_df$BCRcombi %in% LightCombi ), ]$Heavy.LastParent)
        
        #NB: parent has to be paired to be informative, otherwise considered as virtual intermediary 
        # checking that parents are in combi 
        for( current_parent in heavy_tree_parent ){
          if(length( parent_df[which(parent_df$HeavySeq %in% current_parent ), ]$Heavy.LastParent)  ==0   ){
            #changing parent by paired grand parent 
            heavy_tree_parent[heavy_tree_parent == current_parent] <-  unique(parent_df[which(parent_df$Heavy.LastParent %in% current_parent ), ]$Heavy.GrandParent) 
          }
        }
        
        # penalty if remaining combis are not sisters or parent from one another 
        if(length(   heavy_tree_parent[which(heavy_tree_parent %notin% c(current_heavy_chains, "heavyancestor" , ".") )]     ) > 0 ){   return(max_penalty)   }
        
      }
      
      
      
    }
    
    # removing unpaired data 
    #HeavyParents_all_combis = HeavyParents_all_combis[which(grepl("\\.", HeavyParents_all_combis))]
    #LightParents_all_combis = LightParents_all_combis[which(grepl("\\.", LightParents_all_combis))]
    ##HeavyParents_all_combis[lengths(HeavyParents_all_combis)!=0]
    ##LightParents_all_combis[lengths(LightParents_all_combis)!=0]
    
    # Last check on order eg. inversions for concordance criteria
    #order1 = as.vector( na.omit( match(HeavyParents_all_combis, LightParents_all_combis) ) )
    #order2 = as.vector( na.omit( match(LightParents_all_combis, HeavyParents_all_combis) ) )
    
    
    ###LightChains_heavylineage = sub("\\..*", "", LightParents_all_combis)
    ###HeavyChains_lightlineage = sub(".*\\.", "", HeavyParents_all_combis)
    
    ##LightChains_heavylineage = match(sub("\\..*", "", LightParents_all_combis), unique(sub("\\..*", "", LightParents_all_combis))    )
    ##HeavyChains_lightlineage = match(sub(".*\\.", "", HeavyParents_all_combis), unique(sub(".*\\.", "", HeavyParents_all_combis))    )
    
    ## 0 if false 1 if true. order 1 = order of heavy parents in light tree. should be in the right order
    ##is_concordant = min(as.numeric(order1 == sort(order1)  &  order2 == sort(order2))) 
    #is_concordant = all(order1 == sort(order1)   ) &  all(order2 == sort(order2)   ) #& all(LightChains_heavylineage == sort(LightChains_heavylineage)) &  all(HeavyChains_lightlineage == sort(HeavyChains_lightlineage) )
    ##is_concordant = min(as.numeric(order1 == 1:length(order1)  &  order2 == 1:length(order2)  )) #more general eg. ancestor ... ancestor
    #if(is_concordant==F){return(max_penalty)}
    
    # Overlap score (jaccard)
    Heavy_Lineage <- HeavyParents_all_combis[-which(grepl("heavyancestor", HeavyParents_all_combis) )]  #[-1] # do not take terminal leaf != parent for overlap, neither ancestor
    Light_Lineage <- LightParents_all_combis[-which(grepl("lightancestor", LightParents_all_combis) )]#[-1]
    
    #take only unique seqs from each tree (node abundance does not have an impact on the overlap score )
    Heavy_Lineage  = unique( str_extract(Heavy_Lineage, "heavyseq[0-9]{1,5}") )
    Light_Lineage  = unique( str_extract(Light_Lineage, "lightseq[0-9]{1,5}") )
    
    
    lenH =  length(Heavy_Lineage) 
    lenL = length(Light_Lineage)
    jaccard_score  = (length(intersect(Heavy_Lineage, Light_Lineage)  )  *(lenH*lenL))/(  length(unique(lookup_table$BCRcombi ))*(lenH+lenL)  )  

  return( jaccard_score  + unpaired_penalties/2 ) 
}

