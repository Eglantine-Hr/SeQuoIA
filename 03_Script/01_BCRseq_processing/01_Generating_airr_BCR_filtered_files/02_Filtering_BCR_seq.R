
## LOADING LIBRARIES 

library(optparse)
library(readr)
library(dplyr)
library(stringr)

## PARAMETERIZE OPTIONS TO LAUNCH THE SCRIPT FOR EACH SAMPLE PASSED AS ARG -

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="bcr file in AIRR format", metavar="character"),
  make_option(c("-c", "--csv"), type="character", default=NULL, 
              help="output file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name", metavar="character"),
  make_option(c("-m", "--multi"), type="character", default=NULL, 
              help="multi chain output file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## MANAGE NULL ARG

if (is.null(opt$file) | is.null(opt$out) | is.null(opt$multi) ){
  print_help(opt_parser)
  stop("3 arguments must be supplied (input file name, output file name, multichain cells output file name).n", call.=FALSE)
}

`%notin%` <- Negate(`%in%`)

## LOAD DATA 

bcr.df <- read_tsv(opt$file)
annotations.df  = read.csv(opt$csv ) 
print("Imported")
print( names(bcr.df ))

write.table(bcr.df, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/bcr.df.txt")
write.table(annotations.df, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_10.03/annotations.df.txt")

## FILTER BCR DATA 
# Filter BCR sequences according to constant region presence and keeping 1 heavy and 1 light chain with the most umi_count
if ("c_call" %in% names(bcr.df) ){ 
	bcr.df <- bcr.df %>% 
  filter( c_call != "NA") %>%  # keep only bcr with productive chain and a constant region #productive == TRUE,
  filter(locus == "IGH" & str_extract(c_call, "^.{3}") == "IGH" | locus == "IGK" & str_extract(c_call, "^.{3}") == "IGK" | locus == "IGL" & str_extract(c_call, "^.{3}") == "IGL") %>% # remove incoherent c_call and locus (e.g : c_call = IGKC, locus = IGH)
  filter(str_extract(c_call, "^.{3}") == "IGH" & str_extract(v_call, "^.{3}") == "IGH" | str_extract(c_call, "^.{3}") == "IGK" & str_extract(v_call, "^.{3}") == "IGK" | str_extract(c_call, "^.{3}") == "IGL" & str_extract(v_call, "^.{3}") == "IGL") %>%  # remove incoherent c_call and v_call (e.g : c_call = IGLC2, v_call = IGKV6D-21) 
  filter(str_extract(c_call, "^.{3}") == "IGH" & str_extract(j_call, "^.{3}") == "IGH" | str_extract(c_call, "^.{3}") == "IGK" & str_extract(j_call, "^.{3}") == "IGK" | str_extract(c_call, "^.{3}") == "IGL" & str_extract(j_call, "^.{3}") == "IGL")  # remove incoherent c_call and j_call (e.g : c_call = IGLC2, v_call = IGKV6D-21)
	}


#print( names(annotations.df ))
#columns_to_select = names(annotations.df) %notin% names(bcr.df)
#print( names(annotations.df ))
annotations.df = annotations.df[, which(  names(annotations.df ) %in% c("sequence_id", "contig_id", "barcode", "cell_id" , "umis") )  ] 
if("contig_id" %in% names(annotations.df) & "sequence_id" %notin%  names(annotations.df) ){ annotations.df = annotations.df %>% dplyr::rename( sequence_id  = contig_id)  }
print( names(annotations.df ))
print(dim(bcr.df ))

bcr.df = left_join(bcr.df ,  annotations.df       ) 
print(dim(bcr.df ))

if(  any(grepl("barcode", names(bcr.df))) ){bcr.df$cell_id = bcr.df$barcode} 


# Filter BCR sequences having more than 1 heavy and 1 light chain conserving chain with the highest umi_count
if ("cell_id" %in% names(bcr.df) ){
	bcr.filtered.df <- bcr.df %>%     
  dplyr::mutate(chains = recode(locus, "IGH" = "Heavy" , "IGK" = "Light",  "IGL" = "Light")) %>% # Rename IGH/IGL/IGK in chains column as heavy or light chain
  arrange(cell_id, chains) %>% # order according to cell_ID (increasing), chain (increasing) and umi_count (decreasing) 
  distinct(cell_id, chains, .keep_all = TRUE)  %>%  # only keep 1 heavy chain and 1 light chain per cell (with the max umi count) #,  - umi_count
  dplyr::select(-chains) # Remove chains column
	} else {bcr.filtered.df <- bcr.df }

# Cells with multiple umis are kept 
#bcr.filtered.df = bcr.filtered.df %>% dplyr::filter(umis > 1 )


# store the df with only cells that present more than 1 chain to tag_them 
multichains.cells.df <- setdiff(bcr.df, bcr.filtered.df)

## SAVE THE OUTPUT FILE
write_tsv(bcr.filtered.df, file = opt$out)
write_tsv(multichains.cells.df, file = opt$multi)
