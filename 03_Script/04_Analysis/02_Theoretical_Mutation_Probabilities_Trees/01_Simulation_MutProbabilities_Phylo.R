# ##################################################################################
# This script estimates relative codon mutation probabilities, + rev mutations, starting from NCA 
# ##################################################################################
##################################################################################



par(mar = rep(2, 4))



####################################
#Defining inputs and params 
####################################

# Import inputs 
if( exists( "snakemake"))
{
  WORKING_DIR = snakemake@scriptdir
} else
{
  WORKING_DIR = dirname( funr::sys.script())
}



# Import libraries and functions 
#paramsEnv = new.env();
source( file.path(WORKING_DIR, "00_generalDeps.R") ) #, local = paramsEnv)


### Load parameters
# Define an environment that will contain parameters


# Parameterize option to lauch the script for each sample passed as arg
option_list = list(
  make_option(c("-I", "--inputPath"), type="character", default=NULL, 
              help="input path, clono subfolder", metavar="character"),
  make_option(c("-P", "--PATTERN_TABLE"), type="character", default=NULL, 
              help="table with DNA or AA sequence patterns", metavar="character"),
  make_option(c("-O", "--outputPath"), type="character", default=NULL, 
              help="where to write theoretical probas table", metavar="character"),
  make_option(c("-S", "--SIMU_TYPE"), type="character", default=NULL, 
              help="typoe of simulation: etiher tree based or from ancestor", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

TREE_PATH_CLONO = file.path(opt$inputPath)
PATTERN_TABLE = read.csv(opt$PATTERN_TABLE)
Clono = basename(TREE_PATH_CLONO)
DATASET = basename(   dirname(dirname(TREE_PATH_CLONO)) ) 
print(Clono)
OUTPUT_PATH = file.path(opt$outputPath, Clono)
print(opt$outputPath)
print(OUTPUT_PATH)
if ( !dir.exists( OUTPUT_PATH   )  ){ dir.create( OUTPUT_PATH   , recursive = TRUE, showWarnings = FALSE) }

SIMU_TYPE = opt$SIMU_TYPE

source( file.path(WORKING_DIR, "00_allFunctions.R") ) #, local = paramsEnv)


# Default Parameters (it is not recommended to change them)
N_REPLICATES = 3000  #cf saturation graph 
read_frame = 1 #the reading frame was identified in the 3rd step of BCRseq preprocessing  
#SIMU_TYPE = "Tree_Based" #either Tree_Based or From_Ancestor (not recommended)




#GermlinePath = file.path( dirname(dirname(dirname(dirname(TREE_PATH_CLONO)))), "03_UCA_MRCA_inference_annotation" , Clono )
GermlinePath = str_replace(TREE_PATH_CLONO, "05_BCR_Phylogeny", "03_UCA_MRCA_inference_annotation")
GermlinePath = str_replace(TREE_PATH_CLONO, "02_Phylogeny", "01_BCRseq_processing/03_UCA_MRCA_inference_annotation")

print("GERMLINE PATH") 
print(GermlinePath) 

########################################################################################
# Simulation of mutations (NULL model )
########################################################################################

Clono = basename(TREE_PATH_CLONO)


Selected_Trees = read.table(file.path(TREE_PATH_CLONO, "Selected_Trees.txt"), header =T)

column_names_df = c("ParentNode", "ParentNode_seq", "Clono", "Chain", "replicate",  "MutSinceParent", "Sequence", "Nuc_Index_Absolute", "CondonIndex" ,"isNonSyn", "all_mutated") 
mut_table = read.csv( text =  paste(column_names_df, collapse = ",")  ) 

parent_daugthers_table = read.csv( text =  paste(c("Chain", "ParentNode", "ParentNode_seq", "Daughters"), collapse = ",")   ) #to compute selection score 


ancestral_seq_df = read.csv( text =  paste(c("Chain", "Ancestral_Sequence"), collapse = ",")   )


for(chain in Selected_Trees$Chain){ 
  print(chain)
  nMut_distrib = c()
  current_tree = Selected_Trees[which(Selected_Trees$Chain == chain),]
  
  # Getting BCR relations from newick string 
  multiple_topologies = ifelse( as.numeric(current_tree$n_topo)[1] > 1 | as.numeric(current_tree$Pars_Penalty)[1] != 0 , T, F )
  if( file.info(file.path(TREE_PATH_CLONO, chain, "default.inference.parsimony_forest.p")  )$size > 10000 ){ multiple_topologies = T }
  print("REMOVE VIRTUAL")
  print(multiple_topologies)
  BCR_relations  = get_BCR_relations(current_tree$Tree_Nk, remove_virtual = multiple_topologies) #remove_virtual = multiple_topologies
  if(SIMU_TYPE == "From_Ancestor"){BCR_relations$ParentName = rep("Ancestor", nrow(BCR_relations))  } #old strategy 
  
  #Fetching real and inferred sequences 
  BCR_Sequences = seqinr::read.fasta(file = file.path(TREE_PATH_CLONO, chain,  paste0("tree_", current_tree$Tree_Index , ".fasta")    ), 
                                     as.string = TRUE, set.attributes = FALSE, forceDNAtolower = FALSE) %>% as.data.frame()
  BCR_Sequences = t(BCR_Sequences) %>% as.data.frame() %>% dplyr::rename( Sequence = V1) %>% rownames_to_column("SeqName")
  
  
  #adding germline if needed NEW---------------
  GermlineSeq = as.character(read.table( file.path(GermlinePath,  paste0(chain, "_Absolute_Germline_Insert.txt") ), header = F ))
  print(BCR_Sequences[which(BCR_Sequences$SeqName == "Ancestor"),]$Sequence)
  print(GermlineSeq)
  print(str_count(GermlineSeq, "N"))
  NCA_Germline_distance = stringdist::stringdist(GermlineSeq, BCR_Sequences[which(BCR_Sequences$SeqName == "Ancestor"),]$Sequence, "osa" ) - str_count(GermlineSeq, "N")
  print("Germline to NCA distance = ")
  print(NCA_Germline_distance)
  if( NCA_Germline_distance > 1 & NCA_Germline_distance < 1 ){
  	print(paste0("Germline distance = ", NCA_Germline_distance))
  	BCR_Sequences[nrow(BCR_Sequences)+1,] <- c("Germline", GermlineSeq)
  	BCR_relations[which(BCR_relations$SeqName =="Ancestor"),]$ParentName = "Germline"
  	
  }else{BCR_relations = BCR_relations %>% dplyr::filter(SeqName != "Ancestor" ) } #do not model from germline to ancestor 
  print(BCR_relations)
  
  # Keeping track of ancestral sequences 
   ancestral_seq_df[nrow(ancestral_seq_df)+1,] = c(chain, BCR_Sequences[which(BCR_Sequences$SeqName == "Ancestor"),]$Sequence )
  
  
  for ( ParentNode in unique(BCR_relations$ParentName) ){
    
    
    ParentNode_seq = BCR_Sequences[which(BCR_Sequences$SeqName == ParentNode),]$Sequence
    ChildrenNodes = BCR_relations[which(BCR_relations$ParentName == ParentNode),]$SeqName
    ChildrenNodes = ChildrenNodes[grepl("seq|ncestor", ChildrenNodes)]

    if(length(ChildrenNodes) <1){ next} #no seq directly related to parent usually ancestor 
    
    #Updating lineage table  
    parent_daugthers_table[nrow(parent_daugthers_table)+1,] = c(chain,  ParentNode , ParentNode_seq, paste(ChildrenNodes, collapse = ",") )
    
    #Computing max num of mut to compute 
    maxMut_num = max( stringdist(ParentNode_seq,  BCR_Sequences[which(BCR_Sequences$SeqName %in% ChildrenNodes),]$Sequence  , "osa" ) ) - str_count(ParentNode_seq, "N")
    print(paste0("num of Muts from ", ParentNode, " to ", paste0(ChildrenNodes) , " =         " , maxMut_num ))
    nMut_distrib = c(nMut_distrib, maxMut_num)
    if( min( stringdist(ParentNode_seq,  BCR_Sequences[which(BCR_Sequences$SeqName %in% ChildrenNodes),]$Sequence  , "osa" ) ) - str_count(ParentNode_seq, "N") > 35 ){maxMut_num = 1}
    
    # compute mutation probabilities over several replicates 
    for (replicate in 1:N_REPLICATES){
      #print(replicate) 
      
      rev_muts_potential = c("init") #list of potential rev mut
      Sequence = ParentNode_seq
      all_mutated = c()
      
      
      for (mut in 1:maxMut_num){
        
        
        #Applying mutation model
        Sequence_new = applyMutationModel(Sequence = Sequence, chain = chain, DATASET = TREE_PATH_CLONO)
        
        Nuc_Index_Absolute = mcmapply(function(x, y) which(x != y), strsplit(Sequence_new, ""), strsplit(Sequence, "")) #long!
        
        #here we determine the mutation effect (syn or non syn)
        CondonIndex =  1+ floor( (Nuc_Index_Absolute -1  )/3  )
        isNonSyn = CheckDifferentAA(Sequence, Sequence_new, CondonIndex)
        
        #update tables
        all_mutated = c(all_mutated, Nuc_Index_Absolute)
        mut_table[nrow(mut_table)+1,] = c(ParentNode, ParentNode_seq, Clono,  chain, as.numeric(replicate), mut,  Sequence_new, Nuc_Index_Absolute, CondonIndex,  isNonSyn, paste(all_mutated, collapse = ",") ) 
        Sequence = Sequence_new
        
        
        
      } # end of mutation loop (1-15 for FL)
      
    } # end of replicate loop (1000)
    
  } #end of parent loop (max 30 for FL)

} #end of chain loop (1-2)




# Converting factors to numeric 


mut_table = mut_table %>% dplyr::mutate(across(where(~is.character(.) && grepl("^[0-9]$", .) ) & !c(all_mutated), ~  as.numeric(.x)  ))  #converting character to numeric if relevant




########################################################################################
# Estimating Pmut and Pnon syn
########################################################################################


# Computing cumulative proabilities (summing proba of mutations accumulated until rankmut)
mut_table_cumulative = na.omit(mut_table) %>% group_by(Clono, Chain, replicate) %>% separate_rows( all_mutated,  convert = TRUE) #one row per mutation 
mut_table_cumulative$Nuc_Index_Absolute = mut_table_cumulative$all_mutated


# computing relative mutation probabilities 
Pmut_df = mut_table_cumulative %>% group_by( ParentNode , ParentNode_seq, Clono, Chain, MutSinceParent, Nuc_Index_Absolute ) %>% dplyr::summarise(n = n() ) %>% group_by(Clono, Chain, MutSinceParent) %>% dplyr::mutate(Pmut = n/(N_REPLICATES) )  

## Proba for non synonymous mutations 
Pmut_PnonSyn_df = mut_table_cumulative %>% group_by(ParentNode , ParentNode_seq, Clono, Chain, MutSinceParent, Nuc_Index_Absolute, isNonSyn) %>% dplyr::summarise(PnonSyn = n() )  %>% dplyr::filter(isNonSyn == TRUE) %>% group_by(Clono, Chain, MutSinceParent) %>% dplyr::mutate(PnonSyn = PnonSyn/(N_REPLICATES) ) 
Pmut_PnonSyn_df = full_join(Pmut_df, Pmut_PnonSyn_df)

Pmut_PnonSyn_df[is.na(Pmut_PnonSyn_df)] <- 0 #the mutations at this site never turned out non syn 
Pmut_PnonSyn_df$PnonSyn = Pmut_PnonSyn_df$PnonSyn/Pmut_PnonSyn_df$Pmut


mean_probabilities_df =  Pmut_PnonSyn_df %>% mutate_at(vars(starts_with("Mut")), as.numeric) %>% 
  group_by(ParentNode , ParentNode_seq, Clono, Chain) %>% complete(MutSinceParent, Nuc_Index_Absolute = full_seq(c(1, max( nchar(ParentNode_seq    ) ) ), period = 1), fill = list(n = 0, Pmut = 0, PnonSyn = 0)) # %>% mutate(Clono_Chain = paste0(Clono, "_", Chain)) %>% dplyr::select(-c(Clono, Chain) )





########################################################################################
# Impacted sequence motifs 
########################################################################################

#Translate to aa seq 
mean_probabilities_df$ProtSeq = translateDNA( mean_probabilities_df$ParentNode_seq  ) #long

# Inherited from ancestor 
mean_probabilities_df = PatternsFromParent(mean_probabilities_df, "ParentNode_seq" ,PATTERN_TABLE ) 
mean_probabilities_df = mean_probabilities_df %>% dplyr::mutate( CondonIndex =  1+ floor( (Nuc_Index_Absolute -1  )/3  ) ) 








########################################################################################
# Saving outputs 
########################################################################################



#probabilities of mutation table 
mean_probabilities_df = left_join(  mean_probabilities_df,      parent_daugthers_table ) %>% dplyr::mutate( Site = ifelse(Site == "other", NA, Site) )
#adding ancestral seq and dist to NCA for further analyses 
mean_probabilities_df  = left_join(mean_probabilities_df, ancestral_seq_df) 
write.table(mean_probabilities_df, file = paste0( OUTPUT_PATH,   "/", SIMU_TYPE,  "_mean_probabilities_df.txt"), sep = "\t",  row.names = F, col.names =T) 



p1 = ggplot(mean_probabilities_df[which(mean_probabilities_df$Site != "other" & mean_probabilities_df$MutSinceParent ==1),], aes(Site, Pmut )) +
  geom_violin(width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = ParentNode_seq, color = ParentNode_seq ), #trim = FALSE,
    binaxis='y', stackdir='center' , dotsize =  0.05,  position = position_jitter(width = 0.001, height = 0.0005)) + 
  theme_bw() + ggtitle(paste0('Nmut since Parent= ',1 )) +  facet_wrap(~Chain, ncol = 1, scales = "free") + theme(legend.position = "none")

p2 = ggplot(mean_probabilities_df[which(mean_probabilities_df$Site != "other" & mean_probabilities_df$MutSinceParent ==max(nMut_distrib)),], aes(Site, Pmut )) +
  geom_violin(width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = ParentNode_seq, color = ParentNode_seq ), #trim = FALSE,
    binaxis='y', stackdir='center' , dotsize =  0.05,  position = position_jitter(width = 0.001, height = 0.0005)) + 
  theme_bw() + ggtitle(paste0('Nmut since Parent= ',max(nMut_distrib) )) +  facet_wrap(~Chain, ncol = 1, scales = "free")  + theme(legend.position = "none")

combined_plot = cowplot::plot_grid(p1, p2, labels=NULL, ncol = 2, nrow = 1)


ggsave( filename =  paste0( OUTPUT_PATH ,  "/", SIMU_TYPE,  "_MutProbabilitiesDistribution.png"), plot =  combined_plot, width = 10, height = 5 ) 







