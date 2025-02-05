# Ploting Summary and saving results 
# -----------------------

## @knitr SavingOutputs



 write.table(chain_seq_name, file = file.path(GCTREE_FOLDER, "chain_seq_name.txt") , row.names = T, sep = "\t") #CR step 
 
 # write.table(airr_table, file = file.path(GCTREE_FOLDER, "airr_table.txt") , row.names = T, sep = "\t") #CR step 
 
 if( length(available_chains)>1 ){
 	write.table(cell_seq_IDs, file = file.path(GCTREE_FOLDER, "cell_seq_IDs.txt") , row.names = T, sep = "\t") #pairs step 
 	write.table(lookup_table, file = file.path(GCTREE_FOLDER, "lookup_table.txt") , row.names = T, sep = "\t") #pairs step 
 
 }
 

  
 ### Saving selected trees 

##### If trees are not parsimonious or with several ambiguous intermediary reconstruction (several trees of the same topology), virtual nodes will not be considered

equivalent_topologies = do.call(rbind, all_tree_scores)  %>% dplyr::mutate(GcTree_ll_round= as.character(round(as.numeric(GcTree_ll),6))) %>% 
				group_by(Chain, GcTree_ll_round, Pars_Penalty) %>% dplyr::summarise(n_topo  =  length(GcTree_ll) ) %>% 
				dplyr::mutate( n_topo = ifelse(Pars_Penalty < 0, 100 , n_topo) ) 




Selected_Trees$GcTree_ll_round = as.character(round(as.numeric(Selected_Trees$GcTree_ll),6))
Selected_Trees = left_join( Selected_Trees %>% mutate_all(as.character)  , 
				equivalent_topologies[which(equivalent_topologies$GcTree_ll_round %in% Selected_Trees$GcTree_ll_round ), ] %>% mutate_all(as.character) )
write.table(Selected_Trees, file = file.path(GCTREE_FOLDER, "Selected_Trees.txt") , row.names = T)


### Select fasta with right intermediaries NEW
for(k in 1:nrow(Selected_Trees)){file.copy(file.path(GCTREE_FOLDER, Selected_Trees$Chain[k], paste0("tree_", Selected_Trees$Tree_Index[k], ".fasta") ),
						 file.path( GCTREE_FOLDER, paste0(Selected_Trees$Chain[k], "_tree_", Selected_Trees$Tree_Index[k], ".fasta") ))  }



### Plotting Parisomony forest overview
ggplot(data = Tree_Nk_Table, aes(x =GcTree_ll, y = Mut_Pars, color = Iso_Pars, shape = Productive , alpha = toKeep ) ) +
  geom_point()  + scale_colour_gradientn(colours = c("#12ED3F", "#1252ED") ) +  scale_shape_manual(values = c("TRUE" = 20, "FALSE" = 04) ) + 
  geom_point(data=Selected_Trees, aes(x= as.numeric(GcTree_ll),y=as.numeric(Mut_Pars)),color="red", size=3) +
  theme_bw() +facet_wrap(~Chain, ncol = 2, scales = "free") 
  
  

