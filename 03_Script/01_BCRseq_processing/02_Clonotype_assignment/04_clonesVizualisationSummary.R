#########################################################################
#######      The aim of this script is to vizualise clustering     ######
#######   cutree and give some datatable with clonotype features   ######
#########################################################################

# CLUSTERING CUTREE VIZUALISATION ---------------------------------------------------------
## @knitr clustering_cutree_vizualisation

# Define clustering cutree to vizualise
cutree_to_viz <- list( "heavy_clones_hc"=heavy_clones_hc,"light_clones_hc"=light_clones_hc)
clustree.cutoff <- c(HC_CUTREE,LC_CUTREE)
# Print the dendrogram plot to
for (i in seq_along(cutree_to_viz)){
  cat("#### ", names(cutree_to_viz[i]), "\n\n")
  plot(cutree_to_viz[[i]], ylab = "Height", xlab = "", label = FALSE) ; abline(h = clustree.cutoff[i], col = 2, lty = 2)
  cat("\n\n")
}

# CLONOTYPE DATATABLE ---------------------------------------------------------
## @knitr clonotype_datatable

# Print datatable
#print(colnames(bcr.df.clono))
DT::datatable( bcr.df.clono[,c("cell_id", "heavy_vj_call_short", "heavy_junction_aa", "light_vj_call_short", "light_junction_aa", "cloneID", "heavy_contig_id", "light_contig_id")], #"contig_id",
               filter="top",
               rownames = FALSE,
               colnames = c("Cell ID", "Heavy V-J gene", "Heavy CDR3", "Light V-J gene", "Light CDR3", "Clone ID", "ContigID.Heavy", "ContigID.Light"),
               extensions = c('Buttons'),
               options = list(autoWidth = TRUE, 
                              dom = 'Blfrtip', 
                              buttons = c('csv', 'excel')
               ))

# CLONOTYPE SUMMARY ---------------------------------------------------------
## @knitr clonotype_summary

# Generate summary per clonotype to identify malignant clonotype and detect eventual errors
clones_summary_df <- bcr.df.clono %>% 
  group_by(cloneID) %>% 
  summarize(Percentage=n()) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100) 

# Print datatable
DT::datatable( clones_summary_df, 
               filter="top",
               rownames = FALSE,
               colnames = c("Clone ID", "Percentage of cells (%)"),
               options = list(dom = 't')
)

print(dim(bcr.df))
print(dim(bcr.df.clono))

# SAVE DATA ---------------------------------------------------------
## @knitr savingData

## SAVE THE OUTPUT FILE
##write_tsv(bcr.df.clono, file = SEURAT_OUTPUT_FILE)
#bcr.df.clono2 =  bcr.df.clono %>%  mutate(raw_clonotype_id =   ifelse( grepl('[0-9]',heavy_raw_clonotype_id) , heavy_raw_clonotype_id, light_raw_clonotype_id )) %>% dplyr::select(cloneID,  raw_clonotype_id )
bcr.df.clono2 =  bcr.df.clono %>% dplyr::select(cloneID) #,  raw_clonotype_id 
#bcr.df <- merge(bcr.df, bcr.df.clono, by ="raw_clonotype_id", all.x =TRUE, all.y = FALSE)
print(head(bcr.df.clono2))
bcr.df = left_join(   bcr.df, bcr.df.clono2  )
print(dim(bcr.df))
bcr.df$Sample = SAMPLE_NAME #NEW

print(AIRR_OUTPUT_FILE)


## Stop cluster parallelization 
stopCluster(cl)
stopImplicitCluster()
