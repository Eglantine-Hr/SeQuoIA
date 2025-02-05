#####################################################################
#######  The aim of this script is to compute distance matrix  ######
#######           between cell based on VJ and CDR3            ######
#####################################################################

# COMPUTE DISTANCE AND DETERMINE CLONALITY ---------------------------------------------------------
## @knitr scoring_matrix_initialization

# Initialize score matrice (compare cells pairwise to identify clones in either light or heavy chain)
# Extract V, J and CDR3 heavy and light chains 



# COMPUTE DISTANCE BETWEEN CELLS BASED ON HEAVY CHAIN ---------------------------------------------------------
## @knitr heavy_based_distance 

# Initialize score matrice (compare cells pairwise to identify clones in either light or heavy chain)
# Extract V, J and CDR3 heavy and light chains 
hclones <- as.character(Heavy_Chain_df$cell_id)
IGHVJ<-Heavy_Chain_df$heavy_vj_call_short
HCDR3<-Heavy_Chain_df$heavy_junction_aa

# Parallelized computation of cell similarity between cells
heavy_clonotype_distance_list <- foreach(i = 1:(length(hclones)), .combine = 'comb', .multicombine=TRUE) %dopar% {
  heavy_clonotype_matrix <- double(length(hclones))
  for(j in 1:length(hclones)) {
    if(MODE=="VJ_CDR3"){
      # For each pairwise comparison, if VJ between 2 cells is different then heavy similarity score is 0 (H_score)
      if (IGHVJ[i] != IGHVJ[j]){
        H_score <- 0
      }
      # For each pairwise comparison, if heavy V/J genes are similar compute CDR3 similarity 
      else {
        # Compute a levenstein similarity score between the 2 cell CDR3 
        H_score = DescTools::StrDist(HCDR3[i], HCDR3[j], method = DISTANCE_METHOD, mismatch = DISTANCE_MISMATCH, gap = DISTANCE_GAP, ignore.case = DISTANCE_CASE_SENSITIVITY)
      }
      heavy_clonotype_matrix[j]<-H_score
      
    }else{
      # Compute a levenstein similarity score between the 2 cell CDR3 
      H_score = DescTools::StrDist(HCDR3[i], HCDR3[j], method = DISTANCE_METHOD, mismatch = DISTANCE_MISMATCH, gap = DISTANCE_GAP, ignore.case = DISTANCE_CASE_SENSITIVITY)
      heavy_clonotype_matrix[j]<-H_score
    }
  }
  list(heavy_clonotype_matrix)
}


# Convert into matrix and add cellID as rownames and colname
heavy_clonotype_matrix <-as.matrix(heavy_clonotype_distance_list[[1]])
rownames(heavy_clonotype_matrix)<-hclones
colnames(heavy_clonotype_matrix)<-hclones

# Compute distance matrix from similarity matrix
heavy_distance_matrix<-parApply(cl,heavy_clonotype_matrix,c(1,2),function(x){1-x})

# Transform matrix into dist object
heavy_distance_matrix <- as.dist(heavy_distance_matrix, diag = TRUE)

# Compute hierarchical clustering between clones from distance obj.
heavy_clones_hc <- hclust(heavy_distance_matrix, method = HC_METHOD)
print(heavy_clones_hc)

# Define cluster through cutree based on height (80% similarity)
heavy_clones_clusters <- cutree(heavy_clones_hc, h = HC_CUTREE)

rm(heavy_clonotype_distance_list, heavy_distance_matrix, heavy_clonotype_matrix) # heavy_clones_hc,
gc()

# COMPUTE DISTANCE BETWEEN CELLS BASED ON LIGHT CHAIN ---------------------------------------------------------
## @knitr light_based_distance 

gc()
IGLVJ<-Light_Chain_df$light_vj_call_short
LCDR3<-Light_Chain_df$light_junction_aa
lclones <- as.character(Light_Chain_df$cell_id)

# Parallelized computation of cell similarity between cells
light_clonotype_distance_list <- foreach(i = 1:(length(lclones)), .combine = 'comb', .multicombine=TRUE) %dopar% {
  light_clonotype_matrix <- double(length(lclones))
  for(j in 1:length(lclones)) {
    if(MODE=="VJ_CDR3"){
      if (IGLVJ[i] != IGLVJ[j]){
        L_score <- 0
      }
      # For each pairwise comparison, if light V/J genes are similar compute CDR3 similarity 
      else {
        # Compute a levenstein similarity score between the 2 cell CDR3 
        L_score = DescTools::StrDist(LCDR3[i], LCDR3[j], method = DISTANCE_METHOD, mismatch = DISTANCE_MISMATCH, gap = DISTANCE_GAP, ignore.case = DISTANCE_CASE_SENSITIVITY)
      }
      light_clonotype_matrix[j]<-L_score
    }else{
      
      # Compute a levenstein similarity score between the 2 cell CDR3 
      L_score = DescTools::StrDist(LCDR3[i], LCDR3[j], method = DISTANCE_METHOD, mismatch = DISTANCE_MISMATCH, gap = DISTANCE_GAP, ignore.case = DISTANCE_CASE_SENSITIVITY)
      light_clonotype_matrix[j]<-L_score
    }
  }
  list(light_clonotype_matrix)
}

# Convert into matrix and add cellID as rownames and colname
light_clonotype_matrix <-as.matrix(light_clonotype_distance_list[[1]])
rownames(light_clonotype_matrix)<-lclones
colnames(light_clonotype_matrix)<-lclones

# Compute distance matrix from similarity matrix
light_distance_matrix<-parApply(cl,light_clonotype_matrix,c(1,2),function(x){1-x})

# Transform matrix into dist object
light_distance_matrix <- as.dist(light_distance_matrix, diag = TRUE)

# Compute hierarchical clustering between clones from distance obj.
light_clones_hc <- hclust(light_distance_matrix, method = LC_METHOD)

# Define cluster through cutree based on height (80% similarity, i.e help to determine group of cells belonging to same clones)
light_clones_clusters <- cutree(light_clones_hc, h = LC_CUTREE)




#Memory optimization 
#rm(  LC_METHOD, light_clonotype_distance_list)
rm(light_clonotype_distance_list, light_distance_matrix,  light_clonotype_matrix) #light_clones_hc,
gc()

gc()
