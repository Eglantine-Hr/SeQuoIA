# This scripts provides select scores accross cell groups defined in the config file  
# -----------------------

## @knitr SelectScores_groups


Features_of_Interest = grep("Select|[Ss]eq|Clono|row|barcode|Chain|Node", names(All_Scores_Df), value = T, invert = T)
plotlist = list()




for (current_nomenclature in   Features_of_Interest  ){
	cat("## ", current_nomenclature, "{.tabset}", "\n\n") 
	
	for(current_score in All_Types_Scores ){
	
		
	
	
	
		cat("### ", current_score, "\n")
		
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
       	 if(grepl("Inst",current_score )){current_df = current_df %>% dplyr::filter( as.numeric(MutSinceParent) <2  ) 
       	 }else{current_df = current_df %>% dplyr::filter( as.numeric(MutSinceParent) <15 ) } #NEW
       	 

		## Common plot -------------------------------------------------------------------------------------------------------
		Current_Plot = ggplot(current_df, aes_string(x = current_nomenclature, y = current_score )) + 
  				geom_violin(aes(fill = Current_Cellular_Feature), width = 0.9, color="transparent", alpha = 0.2, inherit.aes = T) + 
  				geom_jitter(aes(color = MutSinceParent, text = text_plotly),  width = 0.03, height = 0.001, size = 0.55) + 
  				facet_wrap(.~ Chain, scales = "free") +
  				theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45, hjust=1) ) +
  				theme(legend.position='none') + ggtitle( paste0( current_score , " vs ", current_nomenclature )  )
  				
  				
  				
  		## Compute P values ------------------------------------------------------------------------------------------------------------
		current_df = current_df %>%  dplyr::filter(is.finite(Current_Score)==T & is.na(Current_Score) ==F) %>% group_by(Chain, Current_Cellular_Feature) %>% dplyr::filter(n_distinct(Current_Score) >1  )
		if( length(unique(current_df$Current_Cellular_Feature)) < 9 & length(unique(current_df$Current_Cellular_Feature))  > 1 ){
			#debugging conversion to character factor by rstatix for numerical features (eg. num of muts) 
			
			#write.table(current_df, "/mnt/DOSI/PMLAB/USERS/Active/Eglantine_HECTOR/bioinfo_projets/test_patterns2/current_df_wilcox.txt")
			#& length(unique( as.numeric(current_df$Current_Score ) )) > 3 
			if(any(grepl(pattern="[a-z]",current_df$Current_Cellular_Feature))== F){current_df$Current_Cellular_Feature = as.numeric(current_df$Current_Cellular_Feature)}
                       
			
			df_p_val <- current_df %>% rstatix::group_by(Chain) %>% 
					rstatix::wilcox_test(Current_Score ~ Current_Cellular_Feature, p.adjust.method = "none") %>% #BH
					rstatix::add_xy_position() 
					
			if(any(grepl(pattern="[a-z]",current_df$Current_Cellular_Feature))== F){df_p_val$group1 = as.numeric(df_p_val$group1)
                                                                        			df_p_val$group2 =  as.numeric(df_p_val$group2)}		
					
			#print(head(df_p_val)	)
			#DT::datatable(df_p_val %>% dplyr::arrange(p) )
			
			#print(htmltools::tagList(datatable(df_p_val)))	
			Current_Plot = Current_Plot +  add_pvalue(df_p_val) 
		}


		plot.new()
		## Convert to plotly ---------------------------------------------------------------------------------------
		if(grepl("ehfgrre", current_nomenclature ) ){ #"[Aa]nno|[Ss]ubset|[Cc]ell|[Pp]heno |[Tt]ime"
			#Current_Plot = ggplotly(Current_Plot,   tooltip="text_plotly",  width = 900, height = 300)
			print(htmltools::tagList( ggplotly(Current_Plot,   tooltip="text_plotly",  width = 900, height = 300) )    ) #interactive plots for important features ie pheno 
		}else{ print(Current_Plot) }
		
		if( length(unique(current_df$Current_Cellular_Feature)) < 9 & length(unique(current_df$Current_Cellular_Feature))  > 1 ){
		print(knitr::kable( head(df_p_val %>% dplyr::arrange(p) %>% dplyr::select(Chain, group1, group2, p.adj ,	p.adj.signif),8  ) ) ) }
		 
		
	
	
		cat("\n\n")

	}

}


