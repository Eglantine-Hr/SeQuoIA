#########################################################################
#######      The aim of this script is to update         	     ######
#######   			10x scRNAseq metadata		     ######
#########################################################################


# Udating metadata with selection scores ---------------------------------------------------------
## @knitr AA_properties


aa_select_scores = read.table(file.path( dirname(ALL_SCORES_TABLE), "mutated_aa_properties.txt"  ) , header = T , row.names = NULL,  fill=TRUE, sep = "," ) 


aa_properties_all = aa_select_scores %>% dplyr::filter(aa.beforeMut != aa.afterMut) %>% dplyr::mutate(Seq_Traceback = paste0(Clono_Chain, " " , SeqName))

max2 = unique(as.numeric(aa_properties_all$SelectScore     ))
min_score = quantile( sort(max2), 0.0001)

max2 = quantile( sort(max2[which(max2 <5)]), 0.99)



cat('##', "Positively selected aa", "{.tabset}", ' \n') 

for(current_region in  c("all", "fwr1" ,"cdr1", "fwr2" ,"cdr2", "fwr3" ,"cdr3" ) ){
	cat('###',current_region,  ' \n')
	
	#cat('###', "Positively selected aa",  ' \n') 
	if(current_region == "all"){current_region <- "cdr1|cdr2|cdr3|fwr1|fwr2|fwr3"}
	aa_properties_all_summary = aa_properties_all %>% dplyr::filter(grepl(current_region, Region)  ) %>% 
				group_by(properties.beforeMut, properties.afterMut) %>% 
                                dplyr::summarise(ScoreMean = mean(SelectScore, na.rm = T)    , #circle fill
                                             ScoreVar = log10(n_distinct(Clono_Chain)),#log10(2-sd(SelectScore, na.rm = T) )   , #circlesize
                                             #n_distinct_clones =  n_distinct(Clono_Chain)   , #text
                                             text_plotly = paste(unique(Seq_Traceback), collapse = "\n") ,
                                             Global = F    #rectheat
                                               ) 
                     

	aa_properties_all_before = aa_properties_all %>% dplyr::filter(grepl(current_region, Region)  ) %>% dplyr::mutate(properties.afterMut = "Total") %>% group_by(properties.beforeMut,properties.afterMut) %>% 
  		dplyr::summarise( ScoreMean = mean(SelectScore, na.rm = T)    , #circle fill
                    		ScoreVar = log10(n_distinct(Clono_Chain)), #log10(2-sd(SelectScore, na.rm = T) )   , #circlesize
                    		#n_distinct_clones =  n_distinct(Clono_Chain)   , #text
                    		text_plotly = "NA", 
                    		Global = T )    #rectheat


	aa_properties_all_after = aa_properties_all %>% dplyr::filter(grepl(current_region, Region)  ) %>% dplyr::mutate(properties.beforeMut = "Total") %>% group_by(properties.beforeMut,properties.afterMut) %>% 
  		dplyr::summarise( ScoreMean = mean(SelectScore, na.rm = T)    , #circle fill
                    		ScoreVar = log10(n_distinct(Clono_Chain)), #log10(2-sd(SelectScore, na.rm = T) )   , #circlesize
                    		#n_distinct_clones =  n_distinct(Clono_Chain)   , #text
                    		text_plotly = "NA",
                    		Global = T )    #rectheat

	aa_properties_all_summary = rbind(aa_properties_all_summary, aa_properties_all_before,aa_properties_all_after )
 
	aa_properties_all_summary[nrow(aa_properties_all_summary)+1,] <-  list("Total", "Total", NA, NA, "NA", T)
 
	
	plot.new()
	
	aaplot = ggplot(aa_properties_all_summary, aes(y = properties.beforeMut, x = properties.afterMut)) +    ## global aes
  	geom_tile(aes(fill = Global)) +         ## to get the rect filled
  	geom_point(aes(colour = ScoreMean, size =ScoreVar, text = text_plotly))  +    ## geom_point for circle illusion
  	scale_size(range = c(3, 15))+             ## to tune the size of circles
  	theme_minimal() + #scale_color_gradient(low = "blue",  high = "orange")  + ## color of the circles
  	scale_color_gradient2( low = "forestgreen",mid = "gold", high = "#f1361b", midpoint = 1, limits = c(min_score, max2 ), na.value = "darkred") + #"#e8ce40"
  	scale_fill_manual(values = alpha( c("transparent", "lemonchiffon3"), 0.5)           ) +
 	 theme(axis.text.x = element_text(angle=45, hjust=1)) +guides(size = "none",  fill="none" ) #+
 	 #geom_text(aes(label = n_distinct_clones )
 	 
 	 #print(aaplot)
 	 #print(htmltools::tagList(list(ggplotly(aaplot , tooltip="text_plotly" )))   )
 	 current_plotly = ggplotly(aaplot , tooltip="text_plotly",  width = 700, height = 500  )
 	 cat(' \n \n')
 	 print(htmltools::tagList(plotly::subplot(current_plotly,  nrows=1)))

 	 

 	
 	# Size selection 
  	aa_sizes_all_after = aa_properties_all %>% dplyr::filter(grepl(current_region, Region)  )  %>% 
  				dplyr::mutate( SizeTransition = paste0( size.beforeMut , "->",  size.afterMut ) , Clono =  str_remove_all( Clono_Chain , "[Hh]eavy_|[Ll]ight_"  ) ) %>%
  				group_by(Chain, SizeTransition) %>%  dplyr::filter(n_distinct(SelectScore) > 1) %>% ungroup() 

	Current_Plot = ggplot(aa_sizes_all_after, aes(x = SizeTransition, y = SelectScore )) + 
  				geom_violin(aes(fill = SizeTransition), width = 0.9, color="transparent", alpha = 0.2, inherit.aes = T) + 
  				geom_jitter(aes(color = Clono),  width = 0.03, height = 0.001, size = 0.55) + 
  				facet_wrap(.~ Chain, scales = "free") +
  				theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust=1) ) +
  				theme(legend.position='none') + ggtitle( "Size Transition Selection" )  #+  add_pvalue(df_p_val, label.size = 1, tip.length=0,  bracket.size = 0.2, step.increase = 0 ) 
  	if(nrow(aa_sizes_all_after) > 5 ){print(Current_Plot)}					 
	
	
	cat(' \n \n')
}
htmltools::tagList(plotly::subplot(current_plotly,  nrows=1))

cat('##', "Negatively selected aa", "{.tabset}", ' \n') 
 aa_inheritence_all_summary  = aa_select_scores %>%
 						 dplyr::filter( aa.beforeMut == aa.afterMut & SelectScore < 0 ) %>% #grepl(current_region, Region)  &
 	 					group_by(Clono_Chain, Chain, SeqName , Region, InheritedFrom ) %>% 
 	 					dplyr::summarise(ScoreMean = mean(SelectScore, na.rm = T), Dist_Germline_to_NCA = mean(Dist_Germline_to_NCA, na.rm=T ) ) 


cat('###', "All",  ' \n')



inheritencePlot = ggplot(aa_inheritence_all_summary, aes(y = ScoreMean, x = InheritedFrom)) +
 	 	geom_violin(alpha = 0.3) + geom_point(aes(color = Clono_Chain) ) + geom_line(aes(group = paste(Clono_Chain, SeqName)  ), color = "gray", alpha = 0.05) + 
 	 	theme_minimal() + theme(legend.position = "none") +
 	 	 facet_grid(Chain ~ Region) + theme(axis.text.x = element_text(angle=45, hjust=1)) + ggtitle("Neg Score vs Inheritance")
 	 	  
print(inheritencePlot)




print(names(aa_select_scores))
print("Time" %in% colnames(aa_select_scores))


if (  "Time" %in% colnames(aa_select_scores)  ){
  aa_inheritence_all_summary  = aa_select_scores %>%
 						 dplyr::filter( aa.beforeMut == aa.afterMut & SelectScore < 0 ) %>% #grepl(current_region, Region)  &
 	 					group_by(Clono_Chain, Chain, SeqName , Region, InheritedFrom, Time )  %>% 
 	 					dplyr::summarise(ScoreMean = mean(SelectScore, na.rm = T), Dist_Germline_to_NCA = mean(Dist_Germline_to_NCA, na.rm=T ) ) %>% 		
 	 					dplyr::filter( is.na(Time)==F )
 if(nrow(aa_inheritence_all_summary) < 2 ){break}	 					
 TimeColumn  =  grep("[Tt]ime", colnames(aa_inheritence_all_summary), value = T)[1] #, value = T
#print(TimeColumn)	 					
 	 					
  print(names(aa_inheritence_all_summary))
  #all_time_points = unique( as.vector(aa_inheritence_all_summary[, names(aa_inheritence_all_summary) %in%  TimeColumn  ] ) )
  all_time_points = aa_inheritence_all_summary[, names(aa_inheritence_all_summary) %in%  TimeColumn  ] %>% dplyr::select(1) %>% distinct()
  all_time_points =  unique( as.vector( unlist( all_time_points[,1]) ) )
  all_time_points = all_time_points[is.na(all_time_points) ==F  & all_time_points != "NA" ]
  print(all_time_points)
  
  for(current_TimePoint in all_time_points ){
  
  
  
    cat(' \n \n')
    
    cat("### ", current_TimePoint, "\n\n")
    
    
    currentDF_plot =  aa_inheritence_all_summary[which(aa_inheritence_all_summary$Dist_Germline_to_NCA >= 4),] %>% dplyr::filter( get(TimeColumn) == current_TimePoint )
    
    inheritencePlot = ggplot(currentDF_plot,
                             aes(y = ScoreMean, x = InheritedFrom)) +
      geom_violin(alpha = 0.3) + geom_point(aes(color = Clono_Chain) ) + geom_line(aes(group = paste(Clono_Chain, SeqName)  ), color = "gray", alpha = 0.05) + 
      theme_minimal() + theme(legend.position = "none") +
      facet_grid(Chain ~ Region) + theme(axis.text.x = element_text(angle=45, hjust=1)) + ggtitle("Neg Score vs Inheritance, putative memory ancestor")
    
    if( nrow(currentDF_plot) > 2 ){ print(inheritencePlot) } 
    
    currentDF_plot =  aa_inheritence_all_summary[which(aa_inheritence_all_summary$Dist_Germline_to_NCA < 1),] %>% dplyr::filter( get(TimeColumn) == current_TimePoint )
    #if( nrow(currentDF_plot) <  ){ next  } 
    inheritencePlot = ggplot(currentDF_plot, aes(y = ScoreMean, x = InheritedFrom)) +
      geom_violin(alpha = 0.3) + geom_point(aes(color = Clono_Chain) ) + geom_line(aes(group = paste(Clono_Chain, SeqName)  ), color = "gray", alpha = 0.05) + 
      theme_minimal() + theme(legend.position = "none") +
      facet_grid(Chain ~ Region) + theme(axis.text.x = element_text(angle=45, hjust=1)) + ggtitle("Neg Score vs Inheritance, putative naive ancestor")
    
    if( nrow(currentDF_plot) > 2 ){ print(inheritencePlot) } 
    #print(inheritencePlot)
    
    
    
    cat(' \n \n')
    
  }
  
  
}




  
