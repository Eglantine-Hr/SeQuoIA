######################################################
#####   This file contains custom functions to  ######
#####   compute theoretical mutation probas     ######
######################################################


print("importing functions")

#### General 
`%notin%` <- Negate(`%in%`)

Mode <- function(x, na.rm = FALSE) {
    # it takes two areguments x, and na.rm (https://www.educative.io/answers/what-is-the-mode-method-in-r)
    if(na.rm){ #if na.rm is false it means no need to remove NA values
        x = x[!is.na(x)]
    }
    valx <- unique(x)
    return(valx[which.max(tabulate(match(x, valx)))])
}


integer_breaks <- function(n = 5, ...) {# Ploting only int values
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}



####################################################
# Determining mutations from experimental data
####################################################


get_diff_nuc <- function(Seq1, Seq2){
  
  different_nuc = mcmapply(function(x, y) which(x != y), strsplit(Seq1, ""), strsplit(Seq2, "")) 
  return( paste(different_nuc, collapse = ',') )
  
}



####################################################
# patterns analysis 
####################################################

vectorize_regex_idx <- function(seq, pattern, patternlength, convert_to_nuc = F, after_start=0){
  loc <- str_locate_all(seq, pattern)[[1]] #res <- integer(100)
  if (convert_to_nuc == F){idx <- as.numeric(sapply(loc[,1], function(x) seq(x+after_start, x+after_start+patternlength  ) )) 
  }else{ idx <- as.numeric(sapply(loc[,1], function(x) seq(  ( x*3 -2 )  ,  (x+patternlength)*3   ) ))       }
  return(idx)
}


#-------------------------------------------------------------------------------
# Determines which patterns inherited from NCA are mutated/disrupted 
GetPatternsCurrentSequence <- function(currentDF, NucleotideSeqColumn,  PATTERN_TABLE ){
  # Get mutations to be compared to NCA for pattern analysis
  currentDF$CodonRegex =  paste0( "[( )]" , currentDF$Nuc_Index_Absolute, "[,)]")
  

  
  ############ Input Patterns ########################
  for(pattern in unique(PATTERN_TABLE$Pattern_Name) ){
    print(pattern)
    
    # Check whether input patterns are in NCA 
    current_pattern = PATTERN_TABLE %>% dplyr::filter(Pattern_Name == pattern)
    if(current_pattern$Regex =="*"){current_pattern$Regex <- "\\*"} #"\\ not supported by read table. corresponds to stop codons
    
    if( grepl("Prot|prot|AA|aa|mino", current_pattern$Seq_Type)  ){currentDF$newColumn= as.character(  map( currentDF$ProtSeq , vectorize_regex_idx, current_pattern$Regex,  current_pattern$Pattern_Length - 1,  T)   ) 
    }else{ currentDF$newColumn= as.character( map( currentDF[,which(names(currentDF) == NucleotideSeqColumn )], vectorize_regex_idx, current_pattern$Regex,  current_pattern$Pattern_Length - 1,  F)  ) }
     
    # Matching positions with motifs initially in NCA 
    currentDF  = currentDF %>% dplyr::rowwise() %>% dplyr::mutate(isPattern = grepl(pattern = CodonRegex, x = newColumn )   )

    
    #Renaming temporary columns with right pattern 
    currentDF = currentDF %>% dplyr::mutate(newColumn = ifelse( grepl("numeric", newColumn), NA, newColumn )  )  %>% rename_at( c("newColumn") , list( ~paste0(  c(pattern) ) )    ) # getting rid of numeric 0 and renaming new column to pattern 
    currentDF = currentDF %>% rename_at( c("isPattern") , list( ~paste0( "is" , c(pattern) ) ) )
    
  } #END of pattern loop
  
  # Summarizes in which pattern nucleotide is located
  ## Source code for this part: https://stackoverflow.com/questions/63249304/create-new-column-that-takes-column-name-based-on-value
  currentDF$SiteAfterMut <- paste(apply(currentDF[, grep("^is(?!.*WRC|NNC|SYC|NonSyn|Random)",  names(currentDF),  value = TRUE, perl = TRUE)]     ,1  ,   function(x) names(x)[which(x==T)] )) #bug is only NA

  
  return(currentDF)
} #END of the function


GetPatternsProtSeq0 <- function(currentDF,   PATTERN_TABLE ){
	# Get mutations to be compared to user input sequence patterns
	
	#By default, the pattern is the codon
	currentDF = left_join( currentDF %>% dplyr::mutate( Regex = substr(ProtSeq,CodonIndex,CodonIndex) ) ,
				PATTERN_TABLE %>% dplyr::select(Regex, Site, Color) )
	
	
	 ############ Input Patterns ########################
	 
	 # Selecting non trivial patterns (ie aa) 
	 aminoacids= c('G','A','V','L','I','P','M','F','Y','W','S','T','C','N','Q','D','E','R','H','K')
	 patterns_of_interest = unique(PATTERN_TABLE[-which(PATTERN_TABLE$Regex %in%  aminoacids ),]$Pattern_Name)
	 
	 #Searching for remaining patterns 
	 for(pattern in patterns_of_interest ){
	 	print(pattern)
	 	current_pattern = PATTERN_TABLE %>% dplyr::filter(Pattern_Name == pattern)
	 	#localization <- str_locate_all(ProtSeq, as.character(current_pattern$Regex) )[[1]] #start and end
	 	#pattern_indexes <- as.numeric(sapply(loc[,1], function(x) seq(x, x+current_pattern$Pattern_Length -1  ) ))  #start intermediate indexes end 
	 	currentDF$pattern_indexes= map( currentDF$ProtSeq, vectorize_regex_idx, as.character(current_pattern$Regex),  current_pattern$Pattern_Length - 1)
	 	currentDF = currentDF %>% dplyr::mutate( Site = ifelse( CodonIndex %in%  pattern_indexes, as.character(pattern), as.character(Site) ), #Update Site
	 						  ColorFinal = ifelse( CodonIndex %in%  pattern_indexes, as.character(current_pattern$Color), as.character(Color) )   ) 
	 }
	 
	 currentDF$Color = replace_na(currentDF$ColorFinal, "gray70") #if not in input table, the pattern will be depicted in gray 
	 currentDF$Site = factor(as.character(currentDF$Site), levels = c(PATTERN_TABLE$Pattern_Name, "other") )
	 return(currentDF %>% dplyr::select(-c(pattern_indexes, Regex) ) ) 
}


GetPatternsProtSeq <- function(currentDF,   PATTERN_TABLE ){
	# Get mutations to be compared to user input sequence patterns
	#By default, the pattern is the codon
	currentDF = dplyr::left_join( currentDF %>% dplyr::mutate( Regex = substr(ProtSeq,CodonIndex,CodonIndex) ) ,
				PATTERN_TABLE %>% dplyr::select(Regex, Pattern_Name, Color, Shape) )
	
	
	 ############ Input Patterns ########################
	 
	 # Selecting non trivial patterns (ie aa) 
	 aminoacids= c('G','A','V','L','I','P','M','F','Y','W','S','T','C','N','Q','D','E','R','H','K')
	 patterns_of_interest = unique(PATTERN_TABLE[-which(PATTERN_TABLE$Regex %in%  aminoacids ),]$Pattern_Name)
	 
	 #Searching for remaining patterns 
	 for(pattern in patterns_of_interest ){
	 	print(pattern)
	 	current_pattern = PATTERN_TABLE %>% dplyr::filter(Pattern_Name == pattern)

	 	current_pattern_name = as.character(pattern)
	 	current_pattern_length = as.numeric(current_pattern$Pattern_Length)
	 	currentDF =  currentDF %>%  dplyr::rowwise() %>% dplyr::mutate(pattern_indexes = paste( str_locate_all(ProtSeq, as.character(current_pattern$Regex))[[1]][1], collapse = ","))
	 	currentDF =  currentDF %>%  group_by(ProtSeq) %>% 
	 			dplyr::mutate(  pattern_indexes  = ifelse(is.na(as.numeric(pattern_indexes)), 1000, as.numeric(pattern_indexes) )  ,
	 			ProtSite = ifelse(as.numeric(CodonIndex) >= pattern_indexes & as.numeric(CodonIndex) < pattern_indexes + current_pattern_length , 
	 						current_pattern_name, as.character(Pattern_Name)  ),
	 			ColorFinal = ifelse( as.numeric(CodonIndex) >= pattern_indexes & as.numeric(CodonIndex) < pattern_indexes + current_pattern_length, 
	 						as.character(current_pattern$Color), as.character(Color) )   
	 			 )
	 }
	 
	 currentDF$Color = replace_na(currentDF$ColorFinal, "gray70") #if not in input table, the pattern will be depicted in gray 
	 currentDF$ProtSite = factor(as.character(currentDF$ProtSite), levels = c(PATTERN_TABLE$Pattern_Name, "other") )
	 return(currentDF %>% dplyr::select(-c(pattern_indexes) ) ) 
}



####################################################
# Rendering
####################################################

# Assign colors to sites 
current_palette = c( "gray", "gray", "darkred", "darkblue",  "#33A02C",    "#FF7F00", "#6A3D9A",  "#B15928","magenta", "red", "#B2DF8A", "#FB9A99", "#FDBF6F",  "#CAB2D6")
colors_sites = c("darkred", "darkblue", "gray", "gray",  "orange", "#33A02C",       "#B15928","magenta", "red", "#B2DF8A", "#FB9A99", "#FDBF6F",  "#CAB2D6", "#6A3D9A")

#this functions takes as an input a data frame matching selection scores to codon and region annotation (background) 

plot_condon_scores <- function(Codon_Selection_Scores, BCR_Regions, lines = F, Nucleic_or_Proteic_Site = "ProtSite" ){
	#Preprocessing
	
	
	##For geom_rectangle/CDR FR boundaries -----
	min_index = min(Codon_Selection_Scores$SelectScore, na.rm = T)
	max_index = max(Codon_Selection_Scores$SelectScore, na.rm = T)
	
	## For geom_segment/Inherited by NCA or germline -----
	Codon_Selection_Scores = Codon_Selection_Scores %>% dplyr::mutate(SiteofInterest = get(Nucleic_or_Proteic_Site) ) %>% 
				dplyr::mutate(StartSeg = CodonIndex - 0.25, EndSeg = CodonIndex +0.25, 
				text_plotly = paste("SeqName: ", SeqName, 
                     "\nSeqID: ", SeqID, 
                     "\nCodonScore: ", SelectScore, 
                     "\nParent: ", ParentNode,
                     "\nRegion: ", Region,
                     "\nindex.nuc: ", Nuc_Index_Absolute,
                     "\nindex.aa: ", CodonIndex,
                     "\nSeq.aa: ", ProtSite, 
                     "\nNCA.aa: ", NCA.aa, sep="")   )
				
	## Lanmarks for Post translational modifications sites (arrows) -----
	arrow.length <- 0.5
	touchoff.distance <- -0.9 # distance between data and start of arrow
	arrowhead.size <- 0.9 # in millimeters 
	PTM_sites = 	c("Cys", "GlycSites", "Asn")			
	
	selection_plot = ggplot() + 
	geom_jitter(data = Codon_Selection_Scores,  aes(x =as.factor(CodonIndex) , y = SelectScore, color = SiteofInterest, shape = SiteofInterest ), size = 0.75, width= 0.09, height = 0.015) + 
	geom_rect(data = BCR_Regions, aes(xmin = start, xmax = end, ymin = min_index, ymax = max_index, fill=Region), alpha = 0.08, show.legend = F)+ # BCR regions in the background
  geom_hline(yintercept =  1 ) + geom_hline(yintercept = -1) + geom_hline(yintercept = 0 , linetype="dashed")  +   # Visual guide for positive and negative selection
  theme_bw(base_size = 2) +theme(text = element_text(size=8), axis.text.x = element_text(angle=90, hjust=1), legend.position = "bottom", axis.title.x = element_blank() ) +
  scale_fill_manual( values = rep(c("white", "mediumslateblue"), 4)  ) + facet_wrap(~Chain, ncol = 1, scales = "free") +
  scale_color_manual(values = Patterns_User_Colors) + scale_shape_manual(values= Patterns_User_Shapes) + 
   geom_segment(data = Codon_Selection_Scores[which(Codon_Selection_Scores$InheritedFromGermline ==T),], show.legend = F,aes(x=StartSeg,xend=EndSeg,y=0,yend=0),col="black") +#already in germline
   geom_segment(data = Codon_Selection_Scores[which(Codon_Selection_Scores$ProtSite %in% PTM_sites),], aes(x = CodonIndex, y = min(-1,SelectScore) + touchoff.distance,
                               xend = CodonIndex, yend = min(-1,SelectScore) + touchoff.distance + arrow.length),
               		arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "last"), show.legend = F, color = "black")  #adding arrows above PTMs 
   if(F %in% Codon_Selection_Scores$InheritedFromGermline){selection_plot = selection_plot + geom_segment(data = Codon_Selection_Scores[which(Codon_Selection_Scores$InheritedFromGermline ==F),], show.legend = F,aes(x=StartSeg,xend=EndSeg,y=0,yend=0),col="red") } #red landmark if relevant (otherwise aestetic bug in ggplot)
  
	if (lines == T){
	
	selection_plot  <- Codon_Selection_Scores %>%
       	ggplot( aes( x = CodonIndex, y = SelectScore,  color = as.factor(SeqChainNumber))
       ) + geom_line(   alpha = 0.65, linewidth = 0.2 ) + 
       geom_jitter(  aes(x = CodonIndex, y = SelectScore, text = text_plotly ,  fill = SiteofInterest) , size = 1.2, width= 0.07, height = 0.01, show.legend = F, shape = 21, stroke = 0.05 , alpha = 1,
       		inherit.aes = F ) + 
       facet_wrap(~Chain, ncol = 1, scales = "free")  +
       scale_fill_manual(values = Patterns_User_Colors, guide = "none") +
       guides(fill = F) + 
       geom_hline(yintercept =  1 ) + geom_hline(yintercept = -1) + geom_hline(yintercept = 0 , linetype="dashed") +
       scale_shape_manual(values= Patterns_User_Shapes) +
   	geom_segment(data = Codon_Selection_Scores[which(Codon_Selection_Scores$InheritedFromGermline ==T),], 
   	show.legend = F,aes(x=StartSeg,xend=EndSeg,y=0,yend=0),col="black") +#already in germline
   	geom_segment(data = Codon_Selection_Scores[which(Codon_Selection_Scores$ProtSite %in% PTM_sites),], 
   				aes(x = CodonIndex, y = min(-1,SelectScore) + touchoff.distance,
                               xend = CodonIndex, yend = min(-1,SelectScore) + touchoff.distance + arrow.length),
               		arrow = arrow(length = unit(arrowhead.size, "mm"), ends = "last"), show.legend = F, color = "black")  #adding arrows above PTMs 

	if(F %in% Codon_Selection_Scores$InheritedFromGermline){selection_plot = selection_plot + geom_segment(data = Codon_Selection_Scores[which(Codon_Selection_Scores$InheritedFromGermline ==F),], show.legend = F,aes(x=StartSeg,xend=EndSeg,y=0,yend=0),col="red") } #red landmark if relevant (otherwise aestetic bug)
	}
  
	
       

  
	return(selection_plot)
	
}








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

print("functions loaded")




process_and_merge_fastas_to_df <- function( Experimental_Sequences, IdMap_Files, Isotype_Files, Virtual_Nodes_Sequences  ){
  
  Experimental_Sequences_Id_Names = list()
  
  k = 0
  for(chain in c( "Heavy", "Light" ) ){
    
    
    
    
    #Virtual Sequences
    if(  length(Virtual_Nodes_Sequences[ which(grepl(chain, Virtual_Nodes_Sequences) ) ])    ){
      k = k+1
      print(k)
      
      print(chain)
      current_fasta = readDNAStringSet(Experimental_Sequences[ which(grepl(chain, Experimental_Sequences) ) ])
      current_fasta  = data.table(  SeqID = names(current_fasta ), Sequence = paste(current_fasta) ) %>% dplyr::mutate(Chain = chain)
      print("current_fasta")
      NCAseq = current_fasta[which( grepl("Ancestor",current_fasta$SeqID)     ),]$Sequence[1] #NEW0
      
      current_idmap =    read.table(IdMap_Files[ which(grepl(chain, IdMap_Files) ) ], header = F, sep = ',') %>% dplyr::rename( SeqName = V1, SeqID = V2) %>% 
        separate_rows( SeqID,  convert = TRUE, sep = ":") %>% dplyr::mutate(Chain = chain)
      print("current_idmap")
      Experimental_SeqNames = current_idmap[which(current_idmap$SeqID != ""),]$SeqName #NEW
      
      
      if( file.info(Isotype_Files[ which(grepl(chain, Isotype_Files) ) ])$size >0 ){
        Current_Isotypes =  read.table(Isotype_Files[ which(grepl(chain, Isotype_Files) ) ], header = F, sep = ',')  %>% dplyr::rename( SeqID = V1, Isotype = V2 ) %>% dplyr::mutate(Chain = chain)
      }else{
        Current_Isotypes = as.data.frame(current_idmap$SeqID ) %>% dplyr::mutate( Isotype = NA, Chain =   chain ) %>% dplyr::rename( SeqID = 1) %>% dplyr::filter( nchar(SeqID) > 0  )
      }
      

      print("Current_Isotypes")
      current_idmap = left_join(current_idmap, Current_Isotypes)
      
      
      
      Virtual_Seq_current_chain = readDNAStringSet(Virtual_Nodes_Sequences[ which(grepl(chain, Virtual_Nodes_Sequences) ) ])  
      Virtual_Seq_current_chain  = data.table(  SeqID = names(Virtual_Seq_current_chain ), Sequence = paste(Virtual_Seq_current_chain) ) 
      
      
      Node_Sizes = Virtual_Seq_current_chain  %>% dplyr::mutate( SeqName =str_extract(SeqID, "[A-z0-9]{1,15}(?=[ ])") , 
                                                           SeqName =  ifelse(grepl("[A-z]", SeqName) , SeqName, paste0("X", SeqName)), 
                                                           Abundance = as.numeric(  sub(".*abundance=", "", SeqID)  ),
                                                           Sampled_Ancestor = ifelse( any(SeqName == "Ancestor" & Abundance > 0 ), T, F )   
                                                            ) %>%
                                                  dplyr::select(  SeqName, Abundance,Sampled_Ancestor)
                                            
      
      print(Virtual_Seq_current_chain$SeqID )
      
      
      Virtual_Seq_current_chain = Virtual_Seq_current_chain%>%
        dplyr::filter(grepl("abundance=0", SeqID) ) %>% dplyr::mutate(SeqID = str_extract(SeqID, "[A-z0-9]{1,15}(?=[ ])"), #get virtual sequences and process their names (seq+number)
                                                                      Chain = chain, 
                                                                      SeqName =ifelse(grepl("[A-z]", SeqID) , SeqID, paste0("X", SeqID)) , #add same columns for rbind
                                                                      Isotype = NA, Inferred = T) %>% dplyr::filter( SeqName %notin%  Experimental_SeqNames )
      Experimental_Seq_current_chain = dplyr::inner_join(current_fasta, current_idmap) %>% dplyr::mutate(Inferred = F)
      print(names(Experimental_Seq_current_chain))
      print(names(Virtual_Seq_current_chain))
      Experimental_Seq_current_chain = rbind(Experimental_Seq_current_chain, Virtual_Seq_current_chain)   %>% dplyr::mutate(NCAseq = NCAseq)
      Experimental_Sequences_Id_Names[[k]] <- left_join(Experimental_Seq_current_chain , Node_Sizes ) # NEW
      print("rbind ok")                                                         
      
      
    }
    
    
    
    
  }
  print("loop ok") 
  print(k)
  
  if(k >= 1){Experimental_Sequences_Id_Names = do.call("rbind", Experimental_Sequences_Id_Names)
  }else{Experimental_Sequences_Id_Names  = rbind(Experimental_Seq_current_chain, Virtual_Seq_current_chain) }
  
  return(Experimental_Sequences_Id_Names)
}







####################################################
# PseudoTime analysis
####################################################

#for each Seq, get ancestral sequences seqnames in right order from ancestor to downstream nodes 
getParents <- function(Tree_Nk , SeqOfInterest, real = F ){
  Tree_ape = ape::read.tree(text = Tree_Nk )
  if("Ancestor" %notin% Tree_ape$node.label){ #bug correction for linear trees 
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
  return( paste(Node_Path, collapse = ",") )
}

####################################################
# Sequence Selection Scores 
####################################################


#Instantaneous selection score (KL divergence) normalized by number of mutations
 
ComputeSelectScores = function(current_df, typeofSelection = NULL ){
  current_df$isMutated = ifelse( current_df$OneNucNS ==1, 1, 0)
  ColumnsToGroupby <- grep("odon|ndex|Region|Germline|Site|Seg|OneNucNS|isMutated|Nuc|aa|Select", colnames(current_df), value = T, invert = T)
  	# The probability is the proba that at least one of the nucleotides is mutated, Hence: KL = isMutated*log(P(NS mut in codon)) 
  	current_df = current_df %>%  group_by(SeqName, Chain, CodonIndex) %>% slice(1) %>% ungroup()%>%  
  				group_by(Chain, SeqName)  %>% #group_by_at(ColumnsToGroupby) %>%  #%>% dplyr::filter( isMutated ==1 )  group_by(Chain, SeqName) %>%
			dplyr::summarize(
  					SelectScore_Inst_Pos = sum(-isMutated*log2(  SumProbaCodon   ), na.rm = T )/sum(-SumProbaCodon*log2(SumProbaCodon),  na.rm =T), 
  					SelectScore_Inst_Neg =  sum(-(1-isMutated)*log2( (0.2-SumProbaCodon )), na.rm = T  )/sum( -(1-isMutated)*log2(0.2-SumProbaCodon) ,  na.rm =T)
) %>%
  		dplyr::rowwise()  %>% dplyr::mutate(SelectScore_Inst = sum(c(SelectScore_Inst_Pos, SelectScore_Inst_Neg) , na.rm =T)) 
  
  return(distinct(current_df)  )
}



#mean of immediate selection scores along a branch to take into account apst selection events in current sequence's score 


#mean of immediate selection scores along a branch to take into account apst selection events in current sequence's score 
ComputeCumulatedScore = function(current_df, Selected_Trees, Select_Score_toSum = "Select_Score" ){
	current_df$SelectScoreToSum = current_df[, c(Select_Score_toSum)]
	current_df   = left_join(current_df, Selected_Trees %>% dplyr::select(Chain, Tree_Nk)  ) %>% 
				dplyr::mutate(NodeParents =mcmapply(getParents, Tree_Nk, SeqName, SIMPLIFY = FALSE))

  lineage_df = current_df %>% ungroup() %>% dplyr::select(Chain, SeqName, NodeParents) %>%
    group_by(Chain, SeqName) %>% separate_rows( NodeParents,  convert = TRUE) %>% dplyr::rename(OriginalSeq = SeqName, SeqName = NodeParents) %>% 
    dplyr::mutate(SeqName = as.character( ifelse(grepl('seq|Seq|cestor', SeqName ), SeqName, paste0('X', SeqName)) )    ) #resolve 'x'number added pb 
  lineage_df = left_join( lineage_df, current_df[,!grepl("ID",names(current_df))] %>% ungroup() %>% 
                             dplyr::select(Chain, SeqName, SelectScoreToSum, MutSinceParent) )  %>% distinct() %>% #%>% dplyr::mutate(Global_Select_Score = replace_na(Global_Select_Score, 0))
    group_by(Chain,  OriginalSeq) %>% dplyr::summarise(Cumulated_Select_Score = sum(as.numeric(SelectScoreToSum)*as.numeric(MutSinceParent), na.rm = T ) ) %>%
    dplyr::rename( SeqName = OriginalSeq)
  current_df = left_join(current_df, lineage_df  ) %>% dplyr::rowwise()  %>% dplyr::mutate(Cumulated_Select_Score = sum(c(Cumulated_Select_Score, SelectScoreToSum) , na.rm =T)) 
  return(current_df    %>% dplyr::select(SeqID, Chain, Cumulated_Select_Score) )
}




ComputeCumulatedScores = function(current_df = Global_Selection_Scores){
  lineage_df = current_df %>% ungroup() %>% dplyr::select(Chain, SeqName, NodeParents) %>%
    group_by(Chain, SeqName) %>% separate_rows( NodeParents,  convert = TRUE) %>% dplyr::rename(OriginalSeq = SeqName, SeqName = NodeParents) %>% 
    dplyr::mutate(SeqName = ifelse(grepl('seq|Seq|cestor', SeqName ), SeqName, paste0('X', SeqName)) ) #resolve 'x'number added pb 

  lineage_df = left_join( lineage_df, current_df[,!grepl("ID",names(current_df))] %>% ungroup() %>% 
                             dplyr::select(Chain, SeqName, SelectScore_Inst_Pos, SelectScore_Inst_Neg, SelectScore_Inst) )  %>% distinct() %>% #%>% dplyr::mutate(Global_Select_Score = replace_na(Global_Select_Score, 0))
    group_by(Chain,  OriginalSeq) %>% dplyr::summarise(
    				SelectScore_Cum_Pos 	= mean(SelectScore_Inst_Pos, na.rm = T ),
    				SelectScore_Cum_Neg 	= mean(SelectScore_Inst_Neg, na.rm = T ),
    				SelectScore_Cum	= mean(SelectScore_Inst, na.rm = T )  ) %>%
    dplyr::rename( SeqName = OriginalSeq)
  current_df = left_join(current_df, lineage_df  ) %>% rowwise()  %>% dplyr::mutate( 
  				SelectScore_Cum_Pos	= mean(c(SelectScore_Cum_Pos, SelectScore_Inst_Pos) , na.rm =T),
  				SelectScore_Cum_Neg 	= mean(c(SelectScore_Cum_Neg, SelectScore_Inst_Neg) , na.rm =T),
  				SelectScore_Cum 	= mean(c(SelectScore_Cum, SelectScore_Inst) , na.rm =T),
  				  ) # Past + Immediate Select Scores 
  return(current_df %>% dplyr::select(SeqName, Chain, SelectScore_Cum_Pos, SelectScore_Cum_Neg, SelectScore_Cum)  %>% distinct() )
}


