# This scripts provides select scores accross cell groups defined in the config file  
# -----------------------

## @knitr SelectScores_groups


Features_of_Interest = grep("Select|[Ss]eq|Clono|row|barcode|Chain|Node", names(All_Scores_Df), value = T, invert = T)
plotlist = list()




for (current_nomenclature in   Features_of_Interest  ){
	cat("## ", current_nomenclature, "{.tabset}", "\n\n") 
	
	for(current_score in All_Types_Scores ){
	
		
	
	
	
		cat("### ", current_score, "\n\n")
		
		## Preprocessing --------------------------------------------------------------------------------------------------------------------
		current_df = All_Scores_Df %>% 
  				dplyr::mutate(Current_Score = as.numeric(get(current_score)), 
                		Current_Cellular_Feature = as.character(get(current_nomenclature)) ) %>%
                		dplyr::filter( is.finite(Current_Score) & as.character(Current_Cellular_Feature) != "NaN" ) #filter out Inf


		#Features are sorted alphabetically or according to numerical order 
		if( any(grepl(pattern = "[a-z]", current_df$Current_Cellular_Feature))  ){factor_order = sort(unique(current_df$Current_Cellular_Feature))
		}else{	factor_order = as.numeric( gsub("[^0-9]", "", current_df$Current_Cellular_Feature) ) 
      			factor_order = unique(current_df$Current_Cellular_Feature[order(     factor_order     )]  )   }

		current_df = current_df %>% dplyr::mutate(
                		Current_Cellular_Feature = factor(Current_Cellular_Feature, levels = unique(factor_order)  ) , 
                 		text_plotly = paste("\nClonotype: ", Clonotype,
                                     "\nScore: ", current_score, 
                                     #"\nSeqName: ", SeqName,  
                                     "\nNumber of mut. since NCA: ", NmutSinceNCA, 
                                     "\nNumber of mut. since Germline: ", Dist_ToGermline, 
                                     sep="")   )
                                     
             


		## Narrowing time window for instantaneous score  ----------------------------------------------------------------- 
		current_df_intraclonal = current_df  %>%  dplyr::filter(is.finite(Current_Score)==T & is.na(Current_Score) ==F) %>% group_by(Chain, Current_Cellular_Feature) %>% dplyr::filter(n_distinct(Current_Score) >1  )          
       	 if(grepl("Inst",current_score )){current_df = current_df %>% dplyr::filter( as.numeric(MutSinceParent) <2  & NmutSinceNCA > 1  ) #NEW
       	 }else{current_df = current_df %>% dplyr::filter( as.numeric(MutSinceParent) <15 ) } #NEW
       	 

		## Common plot -------------------------------------------------------------------------------------------------------
		current_df$Dist_Germline_NCA =    current_df$Dist_ToGermline - current_df$NmutSinceNCA 
		library(viridis)
		
		Current_Plot = ggplot(current_df, aes_string(x = current_nomenclature, y = current_score )) + 
  				geom_violin(aes(fill = Current_Cellular_Feature), width = 0.9, color="transparent", alpha = 0.2, inherit.aes = T) + 
  				geom_jitter(aes(color = Dist_Germline_NCA),  width = 0.03, height = 0.001, size = 0.55) + #, text = text_plotly
  				facet_wrap(.~ Chain, scales = "free") +
  				theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust=1) ) + #theme(legend.position='none') +
  				 ggtitle( paste0( current_score , " vs ", current_nomenclature )  ) + 
  				scale_color_viridis(option="turbo", limits = c(0, 15)  )  +guides(fill = "none") 
  				
  				
  				
  		## Compute P values ------------------------------------------------------------------------------------------------------------
		current_df = current_df %>%  dplyr::filter(is.finite(Current_Score)==T & is.na(Current_Score) ==F) %>% group_by(Chain, Current_Cellular_Feature) %>% dplyr::filter(n_distinct(Current_Score) >1  )
		if( length(unique(current_df$Current_Cellular_Feature)) < 15 & length(unique(current_df$Current_Cellular_Feature))  > 2 & 
		length(unique(current_df$Current_Score))  > 5 & grepl("[Dd]ist|utSince", current_nomenclature)==F  ){
			#debugging conversion to character factor by rstatix for numerical features (eg. num of muts) 
			
			
			if(any(grepl(pattern="[a-z]",current_df$Current_Cellular_Feature))== F){current_df$Current_Cellular_Feature = as.numeric(current_df$Current_Cellular_Feature)}
                       
			#current_df  =  current_df %>% dplyr::filter( Current_Score <  1.1*as.numeric(quantile( Current_Score, 0.99 ))  )
			df_p_val <- current_df %>% rstatix::group_by(Chain) %>% 
					rstatix::wilcox_test(Current_Score ~ Current_Cellular_Feature, p.adjust.method = "none") %>% #BH
					rstatix::add_xy_position(fun = "max", step.increase = 0.1 ) 
					
			if(any(grepl(pattern="[a-z]",current_df$Current_Cellular_Feature))== F){df_p_val$group1 = as.numeric(df_p_val$group1)
                                                                        			df_p_val$group2 =  as.numeric(df_p_val$group2)}		
					
			#print(head(df_p_val)	)
			#DT::datatable(df_p_val %>% dplyr::arrange(p) )
			
			#print(htmltools::tagList(datatable(df_p_val)))
			if(length(unique(current_df$Current_Cellular_Feature)) < 10 ){Current_Plot = Current_Plot +  stat_pvalue_manual(df_p_val , hide.ns = T,  bracket.nudge.y = -0.5 , tip.length = 0 ) }	
			
		}
		# Correlation plot for distances 
		
		if(grepl("utSinceNCA", current_nomenclature)   ){ 
			Current_Plot  = ggplot(data = current_df %>% group_by(Clonotype, Chain) %>% 
			dplyr::filter(  n_distinct(NmutSinceNCA) > 1 ) %>% group_by(NmutSinceNCA, Chain) %>% dplyr::filter(  n_distinct(SeqName) > 1 ), #dplyr::filter(  n_distinct(NmutSinceNCA) > 1 ) 
					aes_string(x = as.character(current_nomenclature), y = as.character(current_score ) ) ) + #NmutSinceNCA
    					geom_point(aes(fill = Clonotype, color = Clonotype),  
    							position = position_jitter(w = 0.05, h = 0.05), size = 1.5 , stroke = 0.3) +#geom_smooth(method = "lm") + shape = 21,
    					facet_grid(Chain ~ ., scales ="free_y" ) +   guides(size = "legend", colour = "none") +
    					#scale_fill_brewer(palette = "Paired") + 
    					#scale_x_continuous(breaks = integer_breaks()) + 
   					stat_cor( method = "spearman", color = "red", geom = "label", label.y.npc="top", label.x.npc = "left") + 
    					ggtitle(paste0(current_score, " along PseudoTime axis") ) + theme_bw() +  theme(legend.position = "none")
  	
		}


		plot.new()
		## Convert to plotly ---------------------------------------------------------------------------------------
		if(grepl("ehfgrre", current_nomenclature ) ){ #"[Aa]nno|[Ss]ubset|[Cc]ell|[Pp]heno |[Tt]ime"
			#Current_Plot = ggplotly(Current_Plot,   tooltip="text_plotly",  width = 900, height = 300)
			print(htmltools::tagList( ggplotly(Current_Plot,   tooltip="text_plotly",  width = 900, height = 300) )    )} #interactive plots for important features ie pheno 
		#}else{ print(Current_Plot) }
		if( nrow(current_df) > 5 ){ print(Current_Plot) }
		
		if( length(unique(current_df$Current_Cellular_Feature)) < 18 & length(unique(current_df$Current_Cellular_Feature))  > 2 & 
		length(unique(current_df$Current_Score))  > 5 & grepl("[Dd]ist|utSince", current_nomenclature)==F  ){
		#print(knitr::kable( head(df_p_val %>% dplyr::arrange(p) %>% dplyr::select( starts_with("Chain") | starts_with("group") | starts_with("p") )  ,8  ) ) ) }
		print(knitr::kable( head(df_p_val %>% dplyr::arrange(p) %>% dplyr::filter(p.adj <0.1) %>% dplyr::select( starts_with("Chain") | starts_with("group") | starts_with("p") )  ,16  ) ) ) } 
		
		#cat("\n")
		
		
		
		all_values = current_df_intraclonal[, names(current_df_intraclonal) %in%  current_nomenclature  ] %>% dplyr::select(1) %>% distinct()
		all_values =  unique( as.vector( unlist( all_values[,1]) ) )
		#print(all_values)

		
		#print(all_values)
		if(length(all_values)>1 & length(all_values)<15 ) {
			#cat("\n\n")
			pairs_to_test =  combn(all_values, 2)
			#print(pairs_to_test)
			paired_plots_list = list()
			i = 0 
			
			for(k in 1:ncol(pairs_to_test)  ){
				#print(k)
				current_pair = pairs_to_test[,k]
				
	
				paired_df = current_df_intraclonal %>% dplyr::filter( get(current_nomenclature) %in%  current_pair  ) %>%
					group_by(Clonotype, Chain, Current_Cellular_Feature ) %>% dplyr::summarize(Mean_Score = mean(Current_Score, na.rm = T )) %>%
					ungroup() %>% group_by(Clonotype, Chain )  %>% dplyr::filter(n() > 1)
				
				
				#print(dim(paired_df))
				
				
				
				
				if(nrow(paired_df)>7 ){
					Paired_Plot <- ggpaired(paired_df, x = "Current_Cellular_Feature", y = "Mean_Score", id = "Clonotype", 
						#color = "Current_Cellular_Feature", #palette = "jco", 
						line.color = "gray", line.size = 0.4,
						facet.by = "Chain", short.panel.labs = FALSE)
					# Use only p.format as label. Remove method name.
					Paired_Plot = Paired_Plot + stat_compare_means(label = "p.format", paired = TRUE,  method = "wilcox.test", label.y.npc = "bottom") + 
							theme(axis.title.x = element_blank(),  axis.title.y = element_blank(), axis.text.x = element_text(angle=45, hjust=1), text = 								element_text(size=8) ) + expand_limits(y=0.5) 
					
				 	
					i = i+ 1
					#print(i)
					paired_plots_list[[i]] = Paired_Plot
					#print(Paired_Plot)
					if( i%%8 ==0   ){print( multiplot(plotlist = subsetList(paired_plots_list, as.double((i- i%%8 ):i )) , cols = min(ncol(pairs_to_test) , 4) ) ) }
				}
				if(k == ncol(pairs_to_test) & i > 8 ){print( multiplot(plotlist = subsetList(paired_plots_list, as.double((i- i%%8 ):i )) , cols = min(ncol(pairs_to_test) , 4) ) ) }
			} #end of for
			#print( do.call(grid.arrange, c(paired_plots_list, ncol=4) ) )
			#if(ncol(pairs_to_test) < 8 & i > 1 ){print( multiplot(plotlist = paired_plots_list, cols = min(ncol(pairs_to_test) , 4) ) ) }
			#print( multiplot(plotlist = paired_plots_list, cols = 5) )



		
			
			
	
			
		} #end of if
		cat("\n\n")
		

	}

}


