####################################################
##      This script aims at inferring UCA/MCRCA   ##
##           for each clone and annotating        ##
##                    sequence regions            ##       
####################################################

library(optparse)
library(funr)


## PARAMETERIZE OPTIONS TO LAUNCH THE SCRIPT FOR EACH SAMPLE PASSED AS ARG -----------

option_list = list(
  make_option(c("-b", "--BcrData"), type="character", default=NULL, 
              help="TSV object", metavar="character"),
  make_option(c("-c", "--ClonotypeData"), type="character", default=NULL, 
              help="TSV object 2", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample name", metavar="character"),
  make_option(c("-p", "--PatternFile"), type="character", default=NULL, 
              help="File pattern", metavar="character"),
              make_option(c("-v", "--VDJ_SUMARY_PATH"), type="character", default=NULL, 
              help="File pattern", metavar="character"),
  make_option(c("-O", "--OutputFile"), type="character", default=NULL, 
              help="Output dir name", metavar="character"
  )
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

BCRdf = opt$BcrData
ClonotypeDF = opt$ClonotypeData
SAMPLE = opt$sample
PATTERN_FILE = opt$PatternFile
OUTDIR = dirname(opt$OutputFile) #we only give last created file as snakemake output 
VDJ_SUMARY_PATH = opt$VDJ_SUMARY_PATH

source(file.path( dirname( sys.script())  ,  "00_generalDeps.R") )


`%notin%` <- Negate(`%in%`)

## OTHER PARAMETERS

Ranks_To_Keep = 2
Consensus_Threshold = 0.8




#print(ClonotypeDF)
#print(BCRdf)

###############################################################################
#PART1: COMPUTING CLONOTYPES NUMBERS 
###############################################################################

BCRdf = as.data.frame( read_tsv(file = BCRdf )   ) #%>% distinct() #avoiding several identical sequences per cell to compute clonotype size 
ClonotypeDF = as.data.frame( read_tsv(file = ClonotypeDF )   )

print(unique(ClonotypeDF$cloneID))
#break

## Computing Clone Size and Updating Clone name 
MaxPad = nchar(  as.character(  max(table(ClonotypeDF[which(ClonotypeDF$cloneID != "NOT_CLONAL" ),  ]$cloneID))   )  )  #How many 0s to add in name 

#NEW ADDING CLONOTYPE COMPO
print("naming clonotypes") 

if( "heavy_v_gene" %notin% names(ClonotypeDF ) ){ ClonotypeDF =  ClonotypeDF %>% dplyr::mutate(
			heavy_v_gene =  str_remove(heavy_v_call, "\\*[0-9][0-9]" ), 
			#heavy_c_gene =  str_remove(heavy_isotype, "\\*[0-9][0-9]" ),
			light_v_gene =  str_remove(light_v_call, "\\*[0-9][0-9]" )  ) }
			
if( "heavy_c_gene" %notin% names(ClonotypeDF ) & "heavy_isotype" %in% names(ClonotypeDF ) ){ 
		ClonotypeDF =  ClonotypeDF %>% dplyr::mutate(heavy_c_gene =  str_remove(heavy_isotype, "\\*[0-9][0-9]" ))   }
if("heavy_c_gene" %notin% names(ClonotypeDF ) ){ClonotypeDF$heavy_c_gene <- NA } #YEAP				
			
Clonotype_Summary = ClonotypeDF  %>% group_by(cloneID) %>%
	dplyr::summarise(CloneSize = n(),  heavy_v_gene = Mode(na.omit(heavy_v_gene)), light_v_gene = Mode(na.omit(light_v_gene) )) %>%
	dplyr::mutate( N = ifelse( cloneID ==  "NOT_CLONAL" ,CloneSize , str_pad(string = CloneSize, width = MaxPad, pad = "0")  )  )
Clonotype_Summary$CloneName = paste0( sub("TYPE_.*","",Clonotype_Summary$cloneID), "TYPE_n" , Clonotype_Summary$N, "_Id", sub(".*TYPE_","",Clonotype_Summary$cloneID), "_", Clonotype_Summary$heavy_v_gene, ".", Clonotype_Summary$light_v_gene  )





ClonotypeDF = left_join( ClonotypeDF, Clonotype_Summary %>% dplyr::select(!contains("gene"))  ) %>% dplyr::filter( cell_id != 'virtual') # A VERIFIER 
print("1st join ok")

BCR_Clono_DF = ClonotypeDF
BCR_Clono_DF$heavy_Absolute_Sequence = str_replace_all(BCR_Clono_DF$heavy_sequence_alignment, "\\.", "")
BCR_Clono_DF$light_Absolute_Sequence = str_replace_all(BCR_Clono_DF$light_sequence_alignment, "\\.", "")

################################################################################
#PART2: COMPUTING DISTANCES TO UCA 
# In this part, N will not be taken into account (insertions). It will be useful to compute MRCA 
################################################################################


#BCR_Clono_DF = BCR_Clono_DF %>% group_by(germline_alignment) %>%  mutate(countN = str_count(germline_alignment, pattern = "N"), DistToGermline = stringdist( sequence_alignment,germline_alignment  , "osa") -  countN)  


################################################################################
#PART2: Inferring isotypes
# In this part, cells will be classified as "IgM", "IgD", "IgG3", "IgG1", "IgG2", "IgE", "IgA" 
################################################################################

#removing numbers except for IgG https://stackoverflow.com/questions/59225069/a-regex-to-remove-digits-except-for-words-starting-with
BCR_Clono_DF = BCR_Clono_DF %>% mutate(Isotype = str_replace(heavy_c_gene, "IGH", "Ig" ), Isotype = gsub("(IgG\\S+)|\\d+","\\1",Isotype)  )  

################################################################################
#PART3: SAVE & ANNOTATE UCA IMGT AND ABSOLUTE SEQUENCES
################################################################################

VDJ_Usage_All_Clones = read.csv( text =  paste( c("cloneID", "v_call", "d_call", "j_call", "Chain", "Sample", "Clonotype", "Clone_Size"  ) , collapse = ",")  ) 


NCA_Dist_DF = read.csv( text =  paste( c( "sequence_id", "cloneID" , "CloneName", "Chain" , "GermlineAbsolute" ,"DistToGermline_chain", "NCAseqIMGT", "NCAseqAbsolute"  , "DistToNCA", "Absolute_Sequence", "AncestorType", "Isotype" ) , collapse = ",")  ) 


summary_chain = c("heavy", "light") 


print(names(BCR_Clono_DF))


# Filtering out aberrant cells (nmut light - nmut heavy < 5) 
#test = BCR_Clono_DF[which( BCR_Clono_DF$CloneName %in% clones_of_interest ) ,  ] %>% ungroup() 
#print(dim(test))
#test = test %>% dplyr::rowwise() %>% dplyr::filter( ((light_tot_mut_count - heavy_tot_mut_count) < 6 |  heavy_tot_mut_count > 1  | is.na(heavy_tot_mut_count) | is.na(light_tot_mut_count) ) ==F  )
#print(table(test$CloneName))
#print(dim(test))

print(dim(BCR_Clono_DF) )
BCR_Clono_DF = BCR_Clono_DF %>% ungroup() %>% dplyr::rowwise() %>% dplyr::filter( (light_tot_mut_count - heavy_tot_mut_count) < 6 |  heavy_tot_mut_count > 1 ) %>% ungroup()
print(dim(BCR_Clono_DF) )

# only clonal sequeces are selected with at least 2 sequences for phylogeny reconstruction
clones_of_interest = BCR_Clono_DF %>% dplyr::select( CloneSize, CloneName, heavy_sequence, light_sequence) %>% group_by(CloneName) %>% 
		dplyr::mutate(uniqueH = n_distinct(heavy_sequence), uniqueL = n_distinct(light_sequence, na.rm = T) ) %>% 
		dplyr::filter(CloneSize > 1) %>% 
		#dplyr::filter(uniqueH > 2 & (uniqueL > 2 |  uniqueL ==0 ) & CloneSize > 2) %>%
		dplyr::filter(  ( uniqueH > 1 | uniqueH == 0   ) & (uniqueL > 1 |  uniqueL ==0 )    )
		 
clones_of_interest = unique(clones_of_interest$CloneName)[grepl("NOT_CLONAL", unique(clones_of_interest$CloneName) )== F ]
print("REMAINING CLONES ----------------------------------------------------------------------------------------------------")
print(sort(clones_of_interest ))


#BCR_Clono_DF = BCR_Clono_DF %>% group_by(CloneName) %>% dplyr::mutate(uniqueL = n_distinct(light_sequence) )
#print(head(unique(BCR_Clono_DF$light_sequence)))
#print(clones_of_interest)
#print("BEFORE LOOP")


# removing unpaired data for big clones 
BCR_Clono_DF = BCR_Clono_DF %>% dplyr::filter( CloneSize < 20 | ( chains %in% c( "both", "heavy_chain"  )   )  )



BCR_Clono_DF = BCR_Clono_DF %>% ungroup() %>% dplyr::rowwise() %>% dplyr::filter( (light_tot_mut_count - heavy_tot_mut_count) < 5 |  heavy_tot_mut_count > 1 ) %>% ungroup()
print(dim(BCR_Clono_DF[which( BCR_Clono_DF$CloneName %in% clones_of_interest ) ,  ] ))

for (current_clono in  unique(BCR_Clono_DF$CloneName )){
  for(current_chain in summary_chain){
    
    
    
    #Filtering relevant data ---------------
    currentDF =  BCR_Clono_DF[, which( grepl( paste0("[Cc]lone|DistToGermline|Isotype|[Cc]ell|", current_chain ), names(BCR_Clono_DF)    ) )] %>% dplyr::filter(CloneName == current_clono)  
    #print(names(currentDF))
    names(currentDF) =  names(currentDF) %>% gsub(paste0(current_chain, "_")  , "", .) 
    currentDF = currentDF[which(is.na(currentDF$sequence_alignment) == F ),]
    
    
    
    
    #Saving VDJ usage ---------------------------------------
    VDJ_names = currentDF %>% group_by(cloneID) %>% 
    		   dplyr::summarize(v_call = Mode(v_call, na.rm = T), 
    		   		     d_call = Mode(d_call, na.rm = T), j_call = Mode(j_call, na.rm = T)) %>% 
                  dplyr::mutate(across(everything(),~ gsub("[,*].*","", .))) %>%
                  dplyr::mutate( Chain = str_to_title(current_chain),  Sample = as.character(SAMPLE),
                   		  Clonotype = current_clono , Clone_Size = gsub(".*_n(.*)_Id.*", "\\1", Clonotype),
                   		  v_call = ifelse(cloneID == "NOT_CLONAL", NA, v_call), d_call =ifelse(cloneID == "NOT_CLONAL", NA, d_call)  )
    
    VDJ_Usage_All_Clones = rbind(VDJ_Usage_All_Clones, VDJ_names  )
   
    
    if(nrow(currentDF) > 1 & current_clono %in% clones_of_interest){ #recent change bugg fasta part + not interesting othrewise
      print(current_clono)
      OUTDIR2 = file.path( OUTDIR, current_clono) 
      #sometimes, there are spaces/* in the VDj annot, generating bugs in subsequent folder search in bash 
      OUTDIR2 = str_replace_all(str_replace_all(OUTDIR2, "\\*[0-9]{1,2}|[ ]", ""), "\\*[0-9]{1,2}|[ ]", "")
      if ( dir.exists( OUTDIR2  ) == F ){ dir.create( OUTDIR2   , recursive = TRUE, showWarnings = FALSE) } #creating one folder per clonotype
    

      # Differences in the start recovery of the sequence---------------------------------------
      #newstart = max( nchar(gsub("[A-Z].*", "", currentDF$sequence_alignment)) )   #length of largest gap at the beginning. Should be in frame 
      # avoid few outliers to change crop start 
      currentDF = currentDF %>% dplyr::mutate(initial_gap = nchar(gsub("[A-Z].*", "", sequence_alignment))   ) %>% dplyr::filter( initial_gap <= quantile(initial_gap, 0.9)) 
      newstart = max( nchar(gsub("[A-Z].*", "", currentDF$sequence_alignment)) )   #length of largest gap at the beginning. Should be in frame 
      if(newstart >0){newstart = newstart + (3 - newstart%%3 )  }#getting nucleotides in frame only. bug with newstart = 0 #+1 !!!!error in previous version relaunch
      currentDF = currentDF  %>% dplyr::mutate( sequence_alignment = substr(sequence_alignment, newstart, nchar(sequence_alignment) ) , germline_alignment = substr(germline_alignment, newstart , nchar(germline_alignment)  ) )
      
      
      
      # For each clonotype, we're taking most complete UCA (longest and least N)---------------------------------------
      IMGT_Germline_Seq = as.character( currentDF %>% group_by(cloneID, germline_alignment) %>% summarise(n = n()) %>% dplyr::mutate(string_length = nchar(germline_alignment), count_n = str_count(germline_alignment, "N"))  %>%  ungroup() %>%
                                          dplyr::filter(string_length == max(string_length, na.rm = T) ) %>%  #dplyr::arrange(count_n) %>% dplyr::slice(1) %>% dplyr::select(germline_alignment)   )
                                          dplyr::filter(count_n == min(count_n, na.rm = T) ) %>%  dplyr::slice(1) %>% dplyr::select(germline_alignment)   )
      
      
      
      #We retrieve absolute continuous sequences (!= IMGT)---------------------------------------
      Absolute_Germline_Seq = str_replace_all(IMGT_Germline_Seq, "\\.", "")
      
      ## Here we test whether germline (or experimental seq) is productive and extract ORF---------------------------------------
      testing_productivity_df = currentDF %>% dplyr::select(sequence_alignment) %>% dplyr::mutate(sequence_alignment = str_replace_all(sequence_alignment, "\\.|\\-", ""), SeqType ="Exp" )

      testing_productivity_df = rbind(testing_productivity_df, c(Absolute_Germline_Seq, "Germline"))

      
      

      start_list = c()
      Absolute_Germline_Seq_copy = Absolute_Germline_Seq
      for(first_index in 1:5){
      	testing_productivity_df = testing_productivity_df  %>% dplyr::rowwise() %>% dplyr::mutate(isProductive = ifelse(is_productive( sequence_alignment), 1,0 ) )
      	if( is_productive(Absolute_Germline_Seq_copy) ){start_list = c(start_list, sum(testing_productivity_df$isProductive, na.rm=T)  ) 
      	}else{start_list = c(start_list, -1) }
      	if( all(testing_productivity_df$isProductive) ==1 ){break}
      
      	testing_productivity_df = testing_productivity_df %>% dplyr::mutate( sequence_alignment = substr(sequence_alignment, 2, nchar(sequence_alignment)) ) 
      	 
      	Absolute_Germline_Seq_copy = substr(Absolute_Germline_Seq_copy, 2, nchar(Absolute_Germline_Seq_copy))
      	

      } 
      first_index = which.max( start_list  )
      
      
      Absolute_Germline_Seq = substr(Absolute_Germline_Seq, first_index, nchar(Absolute_Germline_Seq) -  nchar(Absolute_Germline_Seq)%%3 ) #1
      
      ## Saving sequences for each clonotype and chain 
      write.table(Absolute_Germline_Seq, file = file.path( OUTDIR2 , paste0( str_to_title(current_chain)  , "_Absolute_Germline.txt")    ), sep = "\t", row.names = F, col.names =F)

      
      # FWR/CDR annotation ---------------------------------------
      BCR_regions_seq = as.data.frame(currentDF) %>% dplyr::filter(nchar(sequence_alignment) == max(nchar(currentDF$sequence_alignment), na.rm = T) ) %>% dplyr::slice(1) %>% dplyr::select( c(fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4) )  
      
      
      ## Counting number of letters + cummulative sum to infer nuc positions involved in each region 
      BCR_regions_Absolute =  Reduce(function(last, current) current+last, x = nchar(  gsub("\\.", "", BCR_regions_seq) ), accumulate = T)  
      BCR_regions_Absolute = rbind( data.table::shift(as.matrix(BCR_regions_Absolute), n=1, fill=0, type="lag")  +1 , BCR_regions_Absolute )
      rownames(BCR_regions_Absolute) = c( "beginning", "end")
      colnames(BCR_regions_Absolute) = colnames(BCR_regions_seq)
      #write.table(BCR_regions_Absolute, file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "_Absolute_BCR_Regions.txt")    ), sep = "\t", row.names = T, col.names =T)


      # IMGT/Absolute equivalence ---------------------------------------
      BCR_regions_seq = as.data.frame( t(BCR_regions_seq))  %>% dplyr::rename(NucSeq = names(.)[1]) %>% rownames_to_column("Region") %>%  dplyr::mutate(NucSeq = strsplit(as.character(NucSeq), "", fixed = TRUE)) %>%
        unnest(NucSeq) %>% mutate(Nuc_Index_IMGT =   row_number() ) 
      BCR_regions_seq_dots = BCR_regions_seq[which(BCR_regions_seq$NucSeq == "."),] %>% mutate(Nuc_Index_Absolute =  NA)
      BCR_regions_seq = BCR_regions_seq[which(BCR_regions_seq$NucSeq != "."),] %>% arrange(Nuc_Index_IMGT) %>% mutate(Nuc_Index_Absolute =  row_number())
      BCR_regions_seq  = rbind(BCR_regions_seq, BCR_regions_seq_dots)  %>% arrange(Nuc_Index_IMGT)
      BCR_regions_seq  = BCR_regions_seq %>% mutate(Sample = SAMPLE, Clono_Chain = paste0(current_clono, "_", str_to_title(current_chain) ))       
      write.table(BCR_regions_seq, file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "_IMGT_Absolute_Regions.txt")    ), sep = "\t", row.names = F)
      
      # Computing MRCA
      currentDF = currentDF   %>%  dplyr::mutate(countN = str_count(germline_alignment, pattern = "N"), DistToGermline_chain = stringdist( sequence_alignment, germline_alignment  , "osa") -  countN)
      
      
      currentDF = currentDF  %>% group_by(cell_id) %>% dplyr::mutate( total_mut_noCDR3 = sum(DistToGermline_chain) - sum(mu_count_cdr3_r) - sum(mu_count_cdr3_s) ) %>% ungroup()
      
      #currentDF$countN = str_count(currentDF$germline_alignment, pattern = "N")
      #currentDF$DistToGermline_chain = stringdist( currentDF$sequence_alignment, currentDF$germline_alignment  , "osa") -  currentDF$countN
      currentDF = currentDF  %>% filter(nchar(sequence_alignment)>10 ) %>% mutate(rank = dense_rank( total_mut_noCDR3), toKeep = ifelse(rank < Ranks_To_Keep + 1 , T, F) ) %>% 
      					dplyr::filter(DistToGermline_chain < 100)  #to remove after?
      
      #defining number of ranks to keep 
      summary_ranks = currentDF %>% dplyr::group_by(rank) %>% dplyr::summarize(distinct_seq = n_distinct(sequence_alignment))
      if( min(summary_ranks[which(summary_ranks$rank == 1),]$distinct_seq ) > 3 ){ Ranks_To_Keep_Auto = 1 
      }else{Ranks_To_Keep_Auto = 2 }
      print("Ranks_To_Keep_Auto = ")
      print(Ranks_To_Keep_Auto)
      
      #Computing NCA
      NCAseq = as.character(ConsensusSequence(DNAStringSet( currentDF[which(currentDF$rank <  Ranks_To_Keep_Auto + 1  ), ]$sequence_alignment ), 
                                              threshold = 0.99, #Consensus_Threshold
                                              ambiguity= FALSE, 
                                              noConsensusChar="+", 
                                              minInformation = 1, 
                                              ignoreNonBases = FALSE,includeTerminalGaps = TRUE) )
      NCAseqIMGT = str_replace_all(NCAseq, "\\-", "\\.")
      write.table(NCAseqIMGT, file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "_IMGT_NCA.txt")    ), sep = "\t", row.names = F, col.names =F)
      
      # Get Absolute NCAseq (no gaps and ORF )---------------------------------------
      NCAseq = str_replace_all(NCAseq, "\\-", "")
      NCAseq = substr(NCAseq, first_index, nchar(Absolute_Germline_Seq) )
     
      not_conserved =stri_locate_all_regex(NCAseq, '[^ATGC]')[[1]][,1] #stri_locate_all_regex(NCAseq, '\\+|R|Y|K|M|S|W|N')[[1]][,1]   
      if( is.na(not_conserved)[1] == F ){for(k in not_conserved ){ substr(NCAseq, k, k) <-  substr(Absolute_Germline_Seq, k , k)  } }  #IMGT_Germline_Seq 
      
      
      
      # if no strict consensus for junction region => allow 50% consensus---------------------------------------
      not_conserved =  stri_locate_all_regex(NCAseq, 'N')[[1]][,1]
      print("NO CONSENSUS IN JUNCTION SEQ")
      if( is.na(not_conserved)[1] == F ){
      #50% consensus will be used to fill in the remaining gaps in the junctional regions 
        NCAseq_50 = as.character(ConsensusSequence(DNAStringSet( currentDF[which(currentDF$toKeep == T   ), ]$sequence_alignment ), 
                                                   threshold = 0.50, #Consensus_Threshold
                                                   ambiguity= FALSE, 
                                                   noConsensusChar="+", 
                                                   minInformation = 1, 
                                                   ignoreNonBases = FALSE,includeTerminalGaps = TRUE) )
        NCAseq_50 = str_replace_all(NCAseq_50, "\\-", "")
        # In some rare cases, the newly generated consensus is not productive. In this case we go choose the least mutated seuquence to resolve the conflict
        NCAseq_1strank = as.character(ConsensusSequence(DNAStringSet( currentDF[which(currentDF$rank == 1   ), ]$sequence_alignment ), 
                                                threshold = 0.50, #Consensus_Threshold
                                                ambiguity= FALSE, 
                                                noConsensusChar="+", 
                                                minInformation = 0.5, 
                                                ignoreNonBases = FALSE,includeTerminalGaps = TRUE) )
        NCAseq_1strank = str_replace_all(NCAseq_1strank, "\\-", "")
        for(k in not_conserved ){
          substr(NCAseq, k, k) <-  substr(NCAseq_50, k , k)
          if(is_productive(NCAseq) == F |  substr(NCAseq_50, k , k) %notin% c("A", "T", "G", "C")  ){substr(NCAseq, k, k) <- substr(NCAseq_1strank, k , k) }
          #if(substr(NCAseq_50, k , k) %notin% c("A", "T", "G", "C")  ){substr(NCAseq, k, k) <- substr(NCAseq_1strank, k , k) } #YEAP
          } 
      }
      
      
      #Replacing N in germline by NCA nucleotides ---------------------------------------
      Absolute_Germline_Insert = unlist(strsplit(Absolute_Germline_Seq,""))
      Absolute_Germline_Insert[which(Absolute_Germline_Insert == "N")] <- unlist(strsplit(NCAseq,""))[which(Absolute_Germline_Insert == "N")]
      Absolute_Germline_Insert  = paste0(Absolute_Germline_Insert,collapse='')
      
      

      
      print(paste0(current_chain, ",", current_clono)) # to follow 

      
      currentDF$Absolute_Sequence = str_replace_all(currentDF$sequence_alignment, "\\.|\\-", "")
      currentDF$Absolute_Sequence = substr(currentDF$Absolute_Sequence ,first_index, nchar(NCAseq) ) #to check ?????

      
      # classifying primary & secondary responses if only one unique seq => take germline insert ---------------------------------------
      if (  min(currentDF$DistToGermline_chain) < 4 | stringdist(NCAseq, Absolute_Germline_Insert, "osa") < 3    ){ 
        AncestorType = "Germline"
        AncestorSequence = Absolute_Germline_Insert
        print(  paste0( current_clono, current_chain, ":    " , AncestorType) )
      }else{AncestorType =  "NCA"
        AncestorSequence = NCAseq } # add nchar(NCAseq) > 0 condition ? 
      
      #AncestorSequence = Absolute_Germline_Insert
      # selecting sequences and generating fasta files
      print(AncestorType)
      
      #remove seq with gaps in mutations 
      #if(quantile(currentDF$DistToNCA, 1) - quantile(currentDF$DistToNCA, 0.995) >= 10  ){  
      #		currentDF = currentDF[which(  currentDF$DistToNCA >= quantile(currentDF$DistToNCA, 0.995)   ),]   }
      
      ## Optimization of cropping threshold 
      #AllSeqSorted = unique(currentDF$Absolute_Sequence[order(nchar(currentDF$Absolute_Sequence) ) ])
      AllSeqSorted = c(NCAseq, unique(currentDF$Absolute_Sequence) )
      AllSeqSorted = AllSeqSorted[order(nchar(AllSeqSorted), AllSeqSorted)]
      Distinct_Vector_Seq = c()
      ##treat rare cases of insertions (crop before 1st)
      #sequence_aligment = Biostrings::AAStringSet(AllSeqSorted)
      #sequence_aligment = as.character(   msa(sequence_aligment, "ClustalOmega")  )
      ##before_insertion = min( replace_na(    str_locate(sequence_aligment, "-")[,1] , 5000)   )
      #before_insertion = min(  str_locate(sequence_aligment, "-|R|Y|S|W|K|M|B|D|H|V|N")[,1] , na.rm = T )
      pc1_short = quantile( nchar(AllSeqSorted) , 0.01 )
      AllSeqSorted = AllSeqSorted[ which(  nchar(AllSeqSorted) >=  pc1_short) ] # removing 1% shortest sequences
      for(index in 1:length(AllSeqSorted)  ){
        #if (index >= before_insertion){ break }
        AllSeqCropped = AllSeqSorted[index:length(AllSeqSorted)]
        AllSeqCropped = substr(AllSeqCropped, 1, nchar(AllSeqCropped[1]) ) #we cropped to min nchar from seqs that are kept
        Distinct_Vector_Seq = c(Distinct_Vector_Seq,  length(unique(AllSeqCropped)) )
      }
      
	
      ## Generating fasta file ---------------------------------------
      # avoid 1 sequence leading to decreased cropindex, especially if 
      #if( max(Distinct_Vector_Seq) - max(Distinct_Vector_Seq[Distinct_Vector_Seq!=max(Distinct_Vector_Seq)] ) < 2 ){
      #			Distinct_Vector_Seq  = Distinct_Vector_Seq[Distinct_Vector_Seq!=max(Distinct_Vector_Seq)] } 
      
      CropIndex = nchar(AllSeqSorted[which.max(Distinct_Vector_Seq)]) #quantile(nchar(currentDF$Absolute_Sequence), 1/( nrow(currentDF) - 1 ))
      #CropIndex = min(CropIndex, before_insertion) #too stringent 
      
      # Filtering out remaining sequences with indels (NEW)
      currentDF$last_nuc_conservation = substr(currentDF$Absolute_Sequence, CropIndex -20, CropIndex )
      currentDF$last_nuc_conservation = stringdist(   currentDF$last_nuc_conservation, substr(NCAseq, CropIndex -20, CropIndex )   , "osa"      )
      write.table(currentDF, file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "currentDF_", CropIndex,  ".txt")    ) )
      
      fastaTable = currentDF %>%  dplyr::filter(last_nuc_conservation < 5  )  # same length, if indel, last nucs are not conserved 
      
      if( nrow(fastaTable) == 0 ){ next } # insertion in all seqs 
      #write.table(fastaTable, file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "fastaTable_TEMPORARY1",  ".txt")    ) )
      fastaTable = fastaTable %>% 
        			dplyr::filter( grepl("-", Absolute_Sequence ) == F  ) %>% 	
      		dplyr::select(sequence_id, Absolute_Sequence) %>%  rbind(  c( "Ancestor", AncestorSequence )  )  %>% arrange(desc(row_number()))  %>%
       	 mutate(Absolute_Sequence = substr(Absolute_Sequence, 1, CropIndex )) %>% dplyr::filter(nchar(Absolute_Sequence)  >= CropIndex) #ancestor at top, same lengths
      write.table(fastaTable, file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "fastaTable_TEMPORARY1",  ".txt")    ) )
        	
            # Saving NCA and germline seq as txt files NEW 
      #write.table(NCAseq, file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "_IMGT_NCA.txt")    ), sep = "\t", row.names = F, col.names =F)
      write.table(NCAseq, file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "_Absolute_NCA.txt")    ), sep = "\t", row.names = F, col.names =F)
      write.table(substr(Absolute_Germline_Insert, 1, CropIndex ), file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "_Absolute_Germline_Insert.txt")    ), 
      				sep = "\t", row.names = F, col.names =F)
      write.table(substr(NCAseq, 1, CropIndex ), file = file.path( OUTDIR2,  paste0( str_to_title(current_chain)  , "_Absolute_NCA.txt")    ), sep = "\t", row.names = F, col.names =F)
        
      if (length(unique(fastaTable[which(fastaTable$sequence_id != "Ancestor"),]$Absolute_Sequence)) <= 1  ){ 
      	AncestorSequence = Absolute_Germline_Insert
      	AncestorType = "Germline"
      	fastaTable =fastaTable %>% dplyr::mutate(   Absolute_Sequence = ifelse( sequence_id == "Ancestor" ,  AncestorSequence, Absolute_Sequence )  ) 
      }
        
      if (length(unique(fastaTable$Absolute_Sequence)) > 2  ){ 
      # add conditions here to build tree here=> at least 2seq and ancestor  [which(fastaTable$sequence_id != "Ancestor"),]
        write.fasta(as.list(fastaTable$Absolute_Sequence), names=fastaTable$sequence_id , 
                    as.string=FALSE, file.out=  file.path( OUTDIR2,  paste0( str_to_title(current_chain), "_Sequences.fasta")    ) )
      }
      
      # Saving Isotypeidmap files ---------------------------------------
      isotype_table = currentDF %>% dplyr::select(sequence_id, Isotype) %>% na.omit()
      write.table(isotype_table , file = file.path( OUTDIR2,  paste0( str_to_title(current_chain), "_isotypemap.txt" )   )  , sep = ",",  col.names = F, row.names = F)  #GCtree format
      
      #Summary ---------------------------------------
      if(nchar(NCAseq) > 0 ){
        currentDF$croppedNCA = substr(matrix(NCAseq, nrow(currentDF)), 1, nchar(currentDF$Absolute_Sequence)  ) #  incomplete sequences are not penalized
        currentDF = data.frame(sequence_id = currentDF$sequence_id, cloneID= currentDF$cloneID, CloneName = currentDF$CloneName, Chain = matrix(str_to_title(current_chain), nrow(currentDF)), GermlineAbsolute =   str_replace_all(currentDF$germline_alignment, "\\.", "") ,  Absolute_Germline_Insert = matrix( Absolute_Germline_Insert, nrow(currentDF)) ,   DistToGermline_chain = currentDF$DistToGermline_chain,  NCAseqIMGT = matrix( NCAseqIMGT, nrow(currentDF)),  NCAseqAbsolute =  matrix(str_replace_all(NCAseq, "\\.", ""), nrow(currentDF)),  DistToNCA = stringdist(currentDF$Absolute_Sequence , currentDF$croppedNCA, "osa") , Absolute_Sequence = currentDF$Absolute_Sequence, AncestorType =  matrix(AncestorType, nrow(currentDF)),  Isotype = currentDF$Isotype  ) 
        NCA_Dist_DF = rbind(NCA_Dist_DF, currentDF) #%>% dplyr::filter(DistToNCA < 60)
        print(dim(NCA_Dist_DF))
      }
      
      
      
      
      
    }
     
    
    
    
    
  }
}

########################################
# Saving BCR contif df updated with NCA 
########################################
print(names(BCRdf))
BCRdf$Chain = "Heavy" #YEAP
if("Chain" %notin% names(BCRdf)){ BCRdf = BCRdf %>% mutate(Chain = ifelse(chain == "IGH", "Heavy", "Light")  ) }
BCRdf = BCRdf %>% distinct()


NCA_Dist_DF = left_join(BCRdf[,-which(names(BCRdf) %in% c("Absolute_Sequence") )], NCA_Dist_DF )  


write.table(NCA_Dist_DF , file = paste0(OUTDIR, "/01_", SAMPLE, "_", PATTERN_FILE, "_contigID_compatible_filtered_bcr_with_mut_aa_prop_clonotypes_NCA.txt"), 
            sep = " ", dec = ".", row.names = F, col.names = TRUE)



###################
#Saving VDJ usage table 
########################


write.table(VDJ_Usage_All_Clones, file = VDJ_SUMARY_PATH, sep = " ", row.names = F, col.names =T)

