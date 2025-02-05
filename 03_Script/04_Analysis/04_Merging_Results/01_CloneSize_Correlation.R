# This scripts plots correlation between clone size and the mean select score, either immediate or cumulated 
# -----------------------

## @knitr CloneSize_Correlation





# Pre processing -------------------

All_Types_Scores = grep("SelectScore", names(All_Scores_Df), value = T)



#Plotting -------------------


for(current_score in All_Types_Scores){
	cat("## ", current_score, "\n\n")
	
	current_df = All_Scores_Df %>% group_by(Sample, Chain, Clonotype, Clone_Size) %>% 
			dplyr::summarise(ScoreMedian = quantile(as.numeric(get(current_score)), 0.5, na.rm = T) ) 
			
	
			
	current_plot = ggplot(current_df, aes_string(x = "Clone_Size", y = "ScoreMedian", color = "Sample" )  ) +
  			geom_point(size = 0.5) + facet_grid( Chain ~ .) + 
  			stat_cor( method = "spearman", color = "red", geom = "label", label.y.npc="top", label.x.npc = "center") +
  			theme_bw() + ggtitle( paste0( current_score, "vs Clone Size"  )  )  
  			
  	print(current_plot)
  	
  	
  	cat("\n\n\n")
  	
  	
}


