print("parsing params")
# Parameterize option to lauch the script for each sample passed as arg
option_list = list(
  make_option(c("-T", "--THEORETICAL_PROBAS"), type="character", default=NULL, 
              help="main input path, theoretical probas + clono", metavar="character"),
  make_option(c("-P", "--PHENOS_OF_INTEREST"), type="character", default=NULL, 
              help="phenotypes from scRNAseq metadata to compare local patterns", metavar="character"),
  make_option(c("-E", "--PATH_EXPERIMENTAL_ANNOTATION_SEQ"), type="character", default=NULL, 
              help="path to experimental sequences and their annotations", metavar="character"), 
  make_option(c("-I", "--PATH_PHYLO"), type="character", default=NULL, 
              help="idmap file to make correspondance between seqIds and SeqNames", metavar="character"),            
  make_option(c("-O", "--OUTPUT_PATH"), type="character", default=NULL, 
              help="where to write theoretical probas table", metavar="character"),
  make_option(c("-M", "--PATTERN_TABLE"), type="character", default=NULL, 
              help="table with DNA or AA sequence motifs", metavar="character"),
  make_option(c("-S", "--SIMU_TYPE"), type="character", default=NULL, 
              help="tree based or from ancestor", metavar="character")            
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


PATH_EXPERIMENTAL_ANNOTATION_SEQ = opt$PATH_EXPERIMENTAL_ANNOTATION_SEQ
PATH_PHYLO = opt$PATH_PHYLO
SIMU_TYPE = opt$SIMU_TYPE

THEORETICAL_PROBAS = opt$THEORETICAL_PROBAS #path + clono 
current_clono = basename(THEORETICAL_PROBAS)
#THEORETICAL_PROBAS = file.path(THEORETICAL_PROBAS, 

OUTPUT_PATH = file.path(opt$OUTPUT_PATH, current_clono)
if ( !dir.exists( OUTPUT_PATH   )  ){ dir.create( OUTPUT_PATH   , recursive = TRUE, showWarnings = FALSE) }


PHENOS_OF_INTEREST = opt$PHENOS_OF_INTEREST
PHENOS_OF_INTEREST = str_split(PHENOS_OF_INTEREST, pattern = "[,:| ]")[[1]] #creates a vector based on config string provided
PHENOS_OF_INTEREST = PHENOS_OF_INTEREST[PHENOS_OF_INTEREST != ""] #remove emmpty elements if any eg user made a mistake in number of spaces

SAMPLE = basename(PATH_PHYLO)
DATASET = basename( dirname(PATH_PHYLO ) )

PATTERN_TABLE = read.csv(opt$PATTERN_TABLE)
