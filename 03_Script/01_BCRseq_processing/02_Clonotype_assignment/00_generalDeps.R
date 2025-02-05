# ##################################################
# Global declarations and libraries for the analysis
# ##################################################

######## R Libraries

library(pander)
library(readr)
library(digest)
library(magrittr)
library(dplyr)
library(stringr)
library(DescTools)
library(stats)
library(parallel)
library(DT)
library(data.table)
library(htmltools)
library(foreach)
library(doParallel)
library(ggplot2)
#library(Seurat)
library(seqinr)



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
