################################################
# Pseudotime analysis 
################################################


## @knitr GLOBAL_Score_Pseudotime




#  Getting either clusters or cell cycle phase to locate proliferating B cells in the pseudo time representation -----
main_column_of_interest = c(names(SeqName_Features)[grepl("ubset|luster|[Aa]nnot", names(SeqName_Features) )] , 
                            names(SeqName_Features)[grepl("[Pp]hase|[Cc]ycle", names(SeqName_Features) )] )[1]
                            

                            
if ( length(main_column_of_interest[is.na(main_column_of_interest)==F ] ) > 0 ){
	SeqName_Features = SeqName_Features %>% dplyr::mutate(CellInfo = !!main_column_of_interest )
} else{SeqName_Features$CellInfo = "no Pheno Info"}

SeqName_Features = SeqName_Features  %>% dplyr::mutate(CellInfo = ifelse(is.na(CellInfo), "no Metadata", as.character(CellInfo))) %>%
					  dplyr::mutate(ModelingPhase = ifelse( ParentNode == "Germline", "Germline to NCA", "after NCA" ) )




# Ploting 

for (current_score in All_Select_Scores ){
  k = k+1 
  cat("## ", current_score, "\n\n")
  if( nrow(SeqName_Features %>% dplyr::filter( is.na(get(current_score)) ==F & is.finite(get(current_score)) )      )    >2 ){
  	current_plot  = ggplot(data = SeqName_Features, aes_string(x = 'Dist_ToGermline', y = as.character(current_score ) ) ) + #NmutSinceNCA
    		geom_point(aes(fill = CellInfo, color = SeqName, shape = ModelingPhase),  
    			position = position_jitter(w = 0.05, h = 0.05), size = 1.5 , stroke = 0.3) +#geom_smooth(method = "lm") + shape = 21,
    		facet_grid(Chain ~ ., scales ="free_y" ) +   guides(size = "legend", colour = "none") +
    		scale_fill_brewer(palette = "Paired") + 
    		scale_x_continuous(breaks = integer_breaks()) + 
   		stat_cor( method = "spearman", color = "red", geom = "label", label.y.npc="top", label.x.npc = "left") + 
    		ggtitle(paste0(current_score, " along PseudoTime axis") ) + theme_bw() #+  theme(legend.position = "none")
  	print(current_plot)
  }
  
  
  cat("\n\n")
}

