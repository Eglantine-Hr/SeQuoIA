#########################################################################
#######      The aim of this script is to vizualise clonal         ######
#######   proportions and estimate total diversity for each sample ######
#########################################################################


# CLONAL COMPOSITION ---------------------------------------------------------
## @knitr Clonal_Diversity




#### Computing clone proportions per chain and sample 
Clonal_Composition_Sample = All_Clones_Df %>% group_by(Sample, Chain ) %>% dplyr::mutate(total_cells = sum(Clone_Size)) %>% #keep track of total number of cells if needed
        ungroup() %>% distinct() %>%
        group_by(Sample, Chain ) %>% dplyr::mutate( total = sum(Clone_Size), proportion = as.numeric(Clone_Size/total )    ) %>%
  ungroup() %>% arrange(desc(proportion)) %>% dplyr::mutate(lab.ypos = cumsum(proportion) - 0.5*proportion) %>%
  dplyr::mutate(Clone_Size = as.numeric(Clone_Size),
                Clone_Cat = case_when(Clone_Size< 5 ~ as.character(Clone_Size),
                                      Clone_Size >4 & Clone_Size<10  ~ "5-10", 
                                      Clone_Size >= 10 & Clone_Size < 20 ~ "10-19" , 
                                      Clone_Size >= 20 & Clone_Size < 50 ~ "20-49" , 
                                      Clone_Size >= 50 & Clone_Size < 100 ~ "50-100" ,
                                      Clone_Size > 100   ~  ">100" ), 
                Clone_Cat = ifelse(cloneID == "NOT_CLONAL", "Not Clonal", Clone_Cat ) )
                                                             
                                                   
Clonal_Composition_Sample$Clone_Cat = factor(Clonal_Composition_Sample$Clone_Cat, levels = unique(Clonal_Composition_Sample$Clone_Cat))

#### Ploting 
my_palette  = c("gray", "#A50026", "#D73027", "#F46D43" ,"#FDAE61", "#FEE090" ,"#FFFFBF", "#E0F3F8", "#ABD9E9" ,"#74ADD1", "#4575B4", "#313695", "#3F007D", "#54278F", "#6A51A3")

ggplot(distinct(Clonal_Composition_Sample), aes(x = 2, y = proportion, fill =Clone_Cat  )) +
  geom_bar(stat = "identity", color = "white", size=0.05 ) + #, color = "black"
  coord_polar(theta = "y", start = 0) +
  #geom_text(aes(y = lab.ypos, label = paste0("n = ", Clone_Size, ", \n", propotion, "%")), color = "white")+  #scale_fill_manual(values = mycols) +
  theme_classic() + #theme(legend.position = "none")   +
  xlim(.5, 2.5) +   facet_grid( Sample ~ Chain)  + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_text(aes(x= 0.5, y=0.5, label = paste0("total=\n", total) ), size=4) +
  scale_fill_manual( values = my_palette )



# CLONAL DIVERSITY ---------------------------------------------------------
## @knitr Chao1_Inference

Clonal_Composition_Sample_metrics =   Clonal_Composition_Sample %>% group_by(Sample, Chain) %>% 
                                      dplyr::mutate(Sobs = n_distinct(Clonotype), 
                                                    detected1 = ifelse(Clonotype == "NOT_CLONAL", Clone_Size, NA ), 
                                                    detected2 = ifelse(Clone_Size ==2, 1, NA ) ) %>%
                                      dplyr::summarise( f1 = sum(detected1, na.rm = T),
                                                        f2 = sum(detected2, na.rm = T),
                                                        Sobs = mean(Sobs, na.rm = T )) %>%
                                      dplyr::mutate(Chao1 = round(Sobs + ( f1*(f1-1) )/(  2*(f2+1)  ))  )


Clonal_Composition_Sample_metrics = data.table::melt(setDT(Clonal_Composition_Sample_metrics %>% dplyr::select(Sample, Chain, Sobs, Chao1)), id.vars = c("Sample", "Chain"), variable.name = "Clonal_Diversity") %>% dplyr::mutate(Clonal_Diversity = ifelse( Clonal_Diversity== "Sobs" , "Observed", "Chao1" ) ) 


ggplot(Clonal_Composition_Sample_metrics, aes(x = paste0(Sample, " " , Chain), y = value, color = Clonal_Diversity  )) + 
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title.x=element_blank() , axis.title.y=element_blank() ) + 
  ggtitle("Total Number of Clones per Sample") + ylim(c(0, max(Clonal_Composition_Sample_metrics$value) ))

