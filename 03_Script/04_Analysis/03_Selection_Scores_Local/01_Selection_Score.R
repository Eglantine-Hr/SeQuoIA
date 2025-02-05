
################################################
# Imports 
################################################


if( exists( "snakemake"))
{
  WORKING_DIR = snakemake@scriptdir
} else
{
  WORKING_DIR = dirname( funr::sys.script())
}

paramsEnv = new.env();

# libraries and functions 
source( file.path(WORKING_DIR, "00_generalDeps.R") )
source( file.path(WORKING_DIR, "00_allFunctions.R") )


#  Inputs 
source( file.path(WORKING_DIR, "00_analysisParams.R") )


IdMap_Files = list.files(file.path( PATH_PHYLO, current_clono), pattern = "*idmap.txt", full.names =T, recursive = T)
Isotype_Files = list.files(file.path( PATH_PHYLO, current_clono), pattern = "*isotypemap.txt", full.names =T, recursive = T)
print(Isotype_Files)

mutation_probas_table = list.files(THEORETICAL_PROBAS, pattern = "probabilities_df.txt", full.names = T, recursive = T)
mutation_probas_table = read.table(file = mutation_probas_table, header = T, sep = "\t" )  
print(names(mutation_probas_table))

Experimental_Sequences = list.files(path = file.path(PATH_EXPERIMENTAL_ANNOTATION_SEQ, current_clono) , 
                                    pattern = "*Sequences.fasta", full.names =T) #heavy and light, original with cell ids

Virtual_Nodes_Sequences = list.files(path = file.path(PATH_PHYLO, current_clono) , pattern = "*.fasta", full.names =T, include.dirs = F)


Sequence_Annotation = list.files(path = file.path(PATH_EXPERIMENTAL_ANNOTATION_SEQ, current_clono) , 
                                    pattern = "*Regions.txt", full.names =T) %>%
                      map_df(function(x) read.table(x, sep = "\t", header = T) %>% dplyr::mutate(Chain =   str_extract(basename(x), "[^_]+")  )) %>%
                      dplyr::mutate(CodonIndex_IMGT = 1+ floor( (Nuc_Index_IMGT -1  )/3  ) , CodonIndex = 1+ floor( (Nuc_Index_Absolute -1  )/3  ) ) %>%
                      dplyr::group_by(Chain, CodonIndex) %>% dplyr::slice(1)




################################################
# Preprocessing
################################################


# Merging lineage and sequence data (Pmut table with parent seq + experimental sequences)---------------------------



Experimental_Sequences = process_and_merge_fastas_to_df(Experimental_Sequences, IdMap_Files, Isotype_Files, Virtual_Nodes_Sequences)
Real_Nodes_List = Experimental_Sequences %>% dplyr::select(Chain, SeqName, Inferred)
write.table(Experimental_Sequences, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "Experimental_Sequences0.txt"), sep = "\t",  row.names = F, col.names =T)


 



# Computing Nmut since last parent ---------------------------

## Correspondance between parent and daughter seq 

Experimental_Sequences =  left_join(Experimental_Sequences,  
                 mutation_probas_table %>% dplyr::select(ParentNode, ParentNode_seq , Chain, Daughters) %>%
                 separate_rows( Daughters,  convert = TRUE) %>% dplyr::rename(SeqName = Daughters) %>% distinct()  ) %>% dplyr::filter(  is.na(ParentNode) ==F )

Experimental_Sequences = Experimental_Sequences %>% dplyr::mutate( MutSinceParent = stringdist(Sequence, substr(ParentNode_seq,1, nchar(Sequence)), "osa") )


## Nmut 
#print(names(Experimental_Sequences))

Experimental_Sequences = Experimental_Sequences %>% dplyr::mutate( NmutSinceNCA = stringdist(Sequence, substr(NCAseq,1, nchar(Sequence)), "osa")) #NEW


#######NEW!!! Taking only NS mut 
Experimental_Sequences = Experimental_Sequences %>% dplyr::mutate( 
				MutSinceParent_NS = stringdist(  translateDNA(Sequence), translateDNA(substr(ParentNode_seq,1, nchar(Sequence))  ), "osa") 
				#NmutSinceNCA_NS = stringdist(    translateDNA(Sequence), translateDNA(substr(NCAseq,1, nchar(Sequence))          ), "osa")
				)


################################################
# Retrieving mutations from experimental data 
################################################




# Determining which sites are mutated -----------------------------------------

Experimental_Sequences = Experimental_Sequences %>% rowwise()  %>% dplyr::mutate( Nuc_Index_Absolute = get_diff_nuc( as.character(Sequence), as.character(ParentNode_seq)) )

Experimental_Sequences = Experimental_Sequences %>% separate_rows( Nuc_Index_Absolute,  convert = TRUE) %>%
                          dplyr::mutate(Nuc_Index_Absolute = as.numeric(Nuc_Index_Absolute),isMutated = T) 

# Determining mutation effect----------------------------------------------
Experimental_Sequences = Experimental_Sequences %>% 
  mutate(CondonIndex = 1+ floor( (Nuc_Index_Absolute -1  )/3  ),
         Parent_ProtSeq = translateDNA(ParentNode_seq),
         ProtSeq = translateDNA(Sequence), 
         isNonSyn = ifelse(!identical(  substr(ProtSeq, CondonIndex, CondonIndex),   substr(Parent_ProtSeq, CondonIndex, CondonIndex)  ), 1, -1)) 
         

# Adding unmutated nuc data----------------------------------------------
Experimental_Sequences  = Experimental_Sequences %>% group_by(across(- c(Nuc_Index_Absolute, isMutated,  isNonSyn )   )) %>%
  complete( Nuc_Index_Absolute = full_seq(c(1, max( nchar(ParentNode_seq    ) ) ), period = 1),  
            fill = list(isMutated = F,  isNonSyn = -1      ) ) %>%
  group_by(SeqID, Chain, Nuc_Index_Absolute) %>% arrange(desc(isMutated)) %>% dplyr::slice(1) %>%   # to remove nucs that are mutated (doublons) 
  ungroup() #%>% dplyr::mutate(CondonIndex = 1+ floor( (Nuc_Index_Absolute -1  )/3  ) )

################################################
# Getting patterns
################################################


Experimental_Sequences = Experimental_Sequences %>% dplyr::mutate(CodonIndex = 1+ floor( (Nuc_Index_Absolute -1  )/3  ) )
Experimental_Sequences = GetPatternsProtSeq(Experimental_Sequences %>% dplyr::select(-isMutated),  PATTERN_TABLE ) #avoid interferring "is" pattern search



################################################
# Computing selection score 
################################################




# Basal probability per mutation level -----------------------------------------
mutation_probas_table = mutation_probas_table %>% 
                        dplyr::select(ParentNode, ParentNode_seq, Chain, MutSinceParent, Nuc_Index_Absolute, Site, Pmut, PnonSyn) %>% 
                        dplyr::mutate(CodonIndex = 1+ floor( (Nuc_Index_Absolute -1  )/3  ) ) %>% 
                        dplyr::group_by(ParentNode, ParentNode_seq, Chain, MutSinceParent) %>% 
                        dplyr::mutate(PbasalMut = mean(Pmut), PbasalEffect = mean(PnonSyn, na.rm = T)) %>% #MutSinceParent*(1/nchar(ParentNode_seq))  median(Pmut)  
                        dplyr::group_by(ParentNode, ParentNode_seq, Chain, MutSinceParent, CodonIndex) %>%
                        dplyr::mutate(Pbasal_Codon = sum(PbasalMut*PbasalEffect)) 
                        


write.table(Experimental_Sequences, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "Experimental_Sequences1.txt"), sep = "\t",  row.names = F, col.names =T)
write.table(mutation_probas_table, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "mutation_probas_table.txt"), sep = "\t",  row.names = F, col.names =T)  

# Merging theory and experimental data ----------------------------------------
Experimental_Sequences = left_join(Experimental_Sequences %>% dplyr::select(-CondonIndex), mutation_probas_table  ) #%>% #NEW as compared to previous version
                         #dplyr::mutate(Site = ifelse( grepl("is", SiteAfterMut), SiteAfterMut, Site   ) ) # updating sites in experimental sequences 


write.table(Experimental_Sequences, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "Experimental_Sequences2.txt"), sep = "\t",  row.names = F, col.names =T)


## Compute probas for codons (union nucleotides)--------------------------------
print("NAMES")
print(names(Experimental_Sequences))
Codon_Selection_Scores = Experimental_Sequences %>% dplyr::group_by(SeqID, SeqName, Chain, ParentNode, MutSinceParent, CodonIndex, NmutSinceNCA, Isotype, Sequence, ProtSeq, Parent_ProtSeq, Inferred, Abundance, Sampled_Ancestor, MutSinceParent_NS ) %>%  # , MutSinceParent_NS,NmutSinceNCA_NS # at least one of the nucleotides of the codon leads to NS mut 
        dplyr::summarise(Site = Mode(Site, na.rm = T), 
        		  ProtSite = Mode(ProtSite, na.rm = T), 
                         OneNucNS = max(isNonSyn), #at least one of the nucleotides in the codon has a nS mut 
                         SumProbaCodon = sum(Pmut*PnonSyn, na.rm = T), #P one of the nuc leads to NS mut in Codon I given Mut and Seq 
                         Pbasal_Codon = mean(Pbasal_Codon )) %>% #there is only one Pbasal value, mean jst to avoid line duplication 
        dplyr::mutate( SelectScore = OneNucNS/(   1+ OneNucNS*(    (SumProbaCodon -  Pbasal_Codon)/(SumProbaCodon + Pbasal_Codon)    )  ), 
        		#SelectScore = SelectScore/MutSinceParent, #NEW METHOD TO TEST !!!! NOT KEPT
                       Site = ifelse(is.na(Site), "other", Site  ) )   



################################################
# Germline vs NCA inheritance
################################################



Germline_table = read.csv(text="Chain,GermlineSeq,NCAseq")
Current_Clono_Germlines = list.files(file.path( PATH_EXPERIMENTAL_ANNOTATION_SEQ, current_clono), pattern = "*_Absolute_Germline.txt", full.names =T, recursive = T) 
Current_Clono_NCAs = list.files(file.path( PATH_EXPERIMENTAL_ANNOTATION_SEQ, current_clono), pattern = "*_Absolute_NCA.txt", full.names =T, recursive = T) 
for(k in 1:length(Current_Clono_Germlines)){ 
	current_chain = str_extract(Current_Clono_Germlines[k], "[^/]+(?=\\_Absolute_Germline.txt$)") #Extract chain regex
	current_germline = as.character( read.table(Current_Clono_Germlines[k], header = F) )
	current_NCA = as.character( read.table(Current_Clono_NCAs[k], header = F) )
	Germline_table[nrow(Germline_table)+1,] <- c(current_chain, current_germline, current_NCA)
}

print(names(left_join(Codon_Selection_Scores , Germline_table)))
Codon_Selection_Scores = left_join(Codon_Selection_Scores , Germline_table) %>% 
			   	dplyr::mutate( Dist_ToGermline = stringdist(GermlineSeq, Sequence) -str_count(GermlineSeq, "N"),
			   			GermlineProt = translateDNA(GermlineSeq), 
			   			InheritedFromGermline = ifelse(identical( substr(GermlineProt,CodonIndex,CodonIndex), substr(Parent_ProtSeq,CodonIndex,CodonIndex)  ), T,F),
			   			NCAProt = translateDNA(NCAseq), 
			   			NCA.aa = substr(NCAProt,CodonIndex,CodonIndex),
			   			current.aa  = substr(ProtSeq, CodonIndex, CodonIndex)
			   			) %>% dplyr::filter( is.na(CodonIndex)==F )

write.table(Codon_Selection_Scores, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "_Codon_Selection_Scores.txt"), sep = "\t",  row.names = F, col.names =T) 


################################################
# Rendering
################################################


#adding unique seqname chain_number (phenotype part) + ascribe  colours to given sites
Codon_Selection_Scores = Codon_Selection_Scores %>% dplyr::mutate(SeqChainNumber = paste0( substr(Chain, 1, 1),  str_extract(SeqName, "[0-9]{1,4}") ) )
Codon_Selection_Scores$Site = factor(Codon_Selection_Scores$Site, levels = c("WRC", "SYC", "other", "NNC",  unique(Codon_Selection_Scores$Site[Codon_Selection_Scores$Site %notin% c("WRC", "SYC", "other", "NNC"  ) ])    ) )


# IMGT Annotation + FWR/CDR ----------------------------------------------
#to change with IMGT ?
BCR_Regions =  Sequence_Annotation %>% dplyr::select(Chain, Region, CodonIndex) %>% group_by(Chain, Region) %>%
  dplyr::summarise(start = min(CodonIndex, na.rm = T) , end = max(CodonIndex, na.rm = T) )   #%>% column_to_rownames(var = 'Region')

BCR_Regions$Region = factor(BCR_Regions$Region, levels = c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"))


Codon_Selection_Scores = left_join(Codon_Selection_Scores, Sequence_Annotation) 



# Plotting select scores locally ----------------------------------------------

Patterns_User_scale_manual = Experimental_Sequences %>% dplyr::arrange(factor(ProtSite, levels = c(PATTERN_TABLE$Pattern_Name, "other") ) ) %>% 
			group_by(ProtSite) %>% dplyr::summarise(Color = Mode(Color), Shape = Mode(Shape) ) %>% distinct()
Patterns_User_Shapes = as.vector(Patterns_User_scale_manual$Shape) 
Patterns_User_Colors = as.vector(Patterns_User_scale_manual$Color)

print(Patterns_User_Shapes)
print(Patterns_User_Colors)

#colors_sites =current_palette #rev(current_palette[1:length(unique(Codon_Selection_Scores$Site)) ])
# Final plot 

print("PLOT")
selection_plot = plot_condon_scores(Codon_Selection_Scores[which(Codon_Selection_Scores$Inferred ==F),], BCR_Regions)


ggsave( filename =  paste0( OUTPUT_PATH ,  "/", SIMU_TYPE, "_MutProbabilitiesDistribution.png"),
        plot =  selection_plot, width = 13, height = 3*length(unique(Codon_Selection_Scores$Chain)) ) 



#PART2: GLOBAL SELECTION SCORES ----------------------------------------------------------------------

# Create output file name
reportOutputFilename = paste0( OUTPUT_PATH ,  "/", SIMU_TYPE,    "_GLOBAL_scores_analysis"   ,  ".html")

print(OUTPUT_PATH)


# Create a unique temporary file name
rmdCopyName = tempfile( pattern = "tempRmdCopy_", tmpdir = WORKING_DIR, fileext = ".Rmd");
# Copy 'Rmd' file to it
stopifnot( all( file.copy( from = file.path( WORKING_DIR, paste0( "02_Selection_Score_GLOBAL", ".Rmd")),
                               to   = rmdCopyName)));






### Render the report using previously built environment (use link to Rmd file)
PHENOS_OF_INTEREST = PHENOS_OF_INTEREST[PHENOS_OF_INTEREST %notin% c("SeqName", "barcode")]





print("SARTING RMD COMPILATION -----------------------------------------------------------------------")
rmarkdown::render( input = rmdCopyName,
                   output_dir = OUTPUT_PATH, #paramsEnv[["PATH_ANALYSIS_OUTPUT"]],
                   output_file  = I( reportOutputFilename),
                   envir = paramsEnv,
                   quiet = FALSE)
    
#exit() #for tests


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
                        							    InheritedFromGermline == T ~ "Germline", #aa.beforeMut == aa.Germline ~ "Germline",
                        							    TRUE ~ "undefined" ), 
                        				Dist_Germline_to_NCA = Dist_ToGermline - NmutSinceNCA
                                                      ) %>%                         #dplyr::filter(aa.beforeMut != aa.afterMut) %>%
                        dplyr::select(SeqName, Chain, aa.beforeMut, aa.afterMut , SelectScore, Clono_Chain, Sample, Region, Chain, InheritedFrom, 
                        		aa.Germline, aa.NCA , Dist_Germline_to_NCA) 
AA_Properties_toSave = left_join( left_join(AA_Properties_toSave, beforemut_df )   ,  aftermut_df )                   

#write.table(AA_Properties_toSave, file = paste0( OUTPUT_PATH, "/", SIMU_TYPE, "_mutated_aa_properties.txt"), sep = "\t",  row.names = F, col.names =T) 

