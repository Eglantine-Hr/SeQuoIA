#####################################################################
#######     The aim of this script is to compute clustering    ######
#######   between cell, and identify clonotype per chain and   ######
#######                 based on both chains                   ######
#####################################################################


# HEAVY CHAIN CLONOTYPE ASSIGNMENT ---------------------------------------------------------
## @knitr heavy_chain_based_assignment 

# Assign clones based on heavy chain
heavy_clones_clusters <- as.data.frame(as.matrix(heavy_clones_clusters))
colnames(heavy_clones_clusters) <- "heavy_clusters"
heavy_clones_clusters$cell_id <- rownames(heavy_clones_clusters)
Heavy_Chain_df <- merge(heavy_clones_clusters,Heavy_Chain_df, by="cell_id")
Heavy_Chain_df$heavy_clusters <- paste0("CLONOTYPE_", Heavy_Chain_df$heavy_clusters)
Heavy_Chain_df <- Heavy_Chain_df %>% 
  group_by(heavy_clusters) %>% 
  mutate(n_clones= n()) %>% 
  ungroup() %>% 
  mutate(heavy_clusters = ifelse(n_clones > 1, heavy_clusters, "NOT_CLONAL"))%>% 
  dplyr::select(-n_clones)

# LIGHT CHAIN CLONOTYPE ASSIGNMENT ---------------------------------------------------------
## @knitr light_chain_based_assignment 

# Assign clones based on light chain
light_clones_clusters <- as.data.frame(as.matrix(light_clones_clusters))
colnames(light_clones_clusters) <- "light_clusters"
light_clones_clusters$cell_id <- rownames(light_clones_clusters)
Light_Chain_df <- merge(light_clones_clusters,Light_Chain_df, by="cell_id")
Light_Chain_df$light_clusters <- paste0("CLONOTYPE_", Light_Chain_df$light_clusters)
Light_Chain_df <- Light_Chain_df %>% 
  group_by(light_clusters) %>% 
  mutate(n_clones= n()) %>% 
  ungroup() %>% 
  mutate(light_clusters = ifelse(n_clones > 1, light_clusters, "NOT_CLONAL")) %>% 
  dplyr::select(-n_clones)


# COMBINING LIGHT AND HEAVY CLONOTYPE ASSIGNMENT --------------------------
## @knitr both_chains_based_clonotype_assignment 

# Merge by cell_id both heavy and light chains

if(exists("Light_Chain_df") ){ merged.bcr.df <- merge(Heavy_Chain_df, Light_Chain_df, by ="cell_id", all.x =TRUE, all.y = TRUE) # combiner les 2 tables
}else{merged.bcr.df = Heavy_Chain_df}
print(dim(merged.bcr.df))

# Combine clonotype
merged.bcr.df$combined.clono <- ifelse((merged.bcr.df$heavy_clusters=="NOT_CLONAL"|merged.bcr.df$light_clusters=="NOT_CLONAL"),"NOT_CLONAL", 
                                       ifelse((is.na(merged.bcr.df$heavy_clusters)|is.na(merged.bcr.df$light_clusters)), NA, paste0("HEAVY_",merged.bcr.df$heavy_clusters,"_LIGHT_",merged.bcr.df$light_clusters)))

print(dim(merged.bcr.df))
# Extract top clonotype
top_clono <- unlist(dimnames( sort(table(merged.bcr.df$combined.clono), decreasing = TRUE)))[!unlist(dimnames( sort(table(merged.bcr.df$combined.clono), decreasing = TRUE))) %in% "NOT_CLONAL"]
top_clono_df <- unique(merged.bcr.df[match(top_clono, merged.bcr.df$combined.clono),c("light_clusters","heavy_clusters", "combined.clono")])


# Re-attach NA clone if one chain has the same clonotype
for(i in 1:nrow(merged.bcr.df)){
  if(!is.na(merged.bcr.df$combined.clono[i])){
    merged.bcr.df$combined.clono[i] <- merged.bcr.df$combined.clono[i]
  } else if(is.na(merged.bcr.df$combined.clono[i]) & is.na(merged.bcr.df$heavy_clusters[i]) & length(top_clono_df$combined.clono[merged.bcr.df$light_clusters[i] == top_clono_df$light_clusters])>0){
    merged.bcr.df$combined.clono[i] <- top_clono_df$combined.clono[merged.bcr.df$light_clusters[i] == top_clono_df$light_clusters][1] 
  } else if(is.na(merged.bcr.df$combined.clono[i]) & is.na(merged.bcr.df$heavy_clusters[i]) & length(top_clono_df$combined.clono[merged.bcr.df$light_clusters[i] == top_clono_df$light_clusters])==0){
    merged.bcr.df$combined.clono[i] <- paste0("LIGHT_",merged.bcr.df$light_clusters[i])
  } else if(is.na(merged.bcr.df$combined.clono[i]) & is.na(merged.bcr.df$light_clusters[i]) & length(top_clono_df$combined.clono[merged.bcr.df$heavy_clusters[i] == top_clono_df$heavy_clusters])>0){
    merged.bcr.df$combined.clono[i] <- top_clono_df$combined.clono[merged.bcr.df$heavy_clusters[i] == top_clono_df$heavy_clusters][1]
  } else if(is.na(merged.bcr.df$combined.clono[i]) & is.na(merged.bcr.df$light_clusters[i]) & length(top_clono_df$combined.clono[merged.bcr.df$heavy_clusters[i] == top_clono_df$heavy_clusters])==0){
    merged.bcr.df$combined.clono[i] <- paste0("HEAVY_",merged.bcr.df$heavy_clusters[i])
  }
}
print(dim(merged.bcr.df))
# Rename Clone ID & add info if cell present both light and heavy chain
merged.bcr.df <- merged.bcr.df %>% 
  arrange(combined.clono) %>% 
  group_by(combined.clono) %>% 
  mutate(n.combined.clono = n()) %>% 
  mutate(combined.clono = ifelse(n.combined.clono > 1, combined.clono, "NOT_CLONAL"))  %>% # rename as not clonal heavy or light remaning "clonotype" with less than 1 occurance
  as.data.table()
merged.bcr.df$n.combined.clono <- NULL
merged.bcr.df[ , cloneID := rleid(combined.clono)]

print(dim(merged.bcr.df))

print(exists("current_donor"))
print(current_donor)
#if (exists("current_donor")){ merged.bcr.df =merged.bcr.df %>% dplyr::mutate(cloneID = ifelse(combined.clono=="NOT_CLONAL","NOT_CLONAL", paste0(current_donor, "CLONOTYPE_",cloneID))) 
#}else{merged.bcr.df =merged.bcr.df %>% dplyr::mutate(cloneID = ifelse(combined.clono=="NOT_CLONAL","NOT_CLONAL", paste0(SAMPLE_NAME, "_CLONOTYPE_",cloneID)))}
merged.bcr.df =merged.bcr.df %>% dplyr::mutate(cloneID = ifelse(combined.clono=="NOT_CLONAL","NOT_CLONAL", paste0( "CLONOTYPE_",cloneID)))

bcr.df.clono <- merged.bcr.df %>% 
  #mutate(cloneID = ifelse(combined.clono=="NOT_CLONAL","NOT_CLONAL", paste0(SAMPLE_NAME, "_CLONOTYPE_",cloneID))) %>% 
  mutate(chains = case_when(
    is.na(heavy_vj_call_short) != TRUE & is.na(light_vj_call_short) == TRUE ~ "heavy_chain",
    is.na(heavy_vj_call_short) == TRUE & is.na(light_vj_call_short) != TRUE ~ "light_chain",
    is.na(heavy_vj_call_short) != TRUE & is.na(light_vj_call_short) != TRUE ~ "both")) %>% 
  dplyr::select(-c(heavy_clusters, light_clusters, combined.clono))


print("DIM bcr.df.clono")
print(dim(bcr.df.clono))


# Add clone ID to AIRR filtered dataframe
#bcr.df <- merge(bcr.df, bcr.df.seurat[,c("cell_id", "cloneID")], by = "cell_id")

# Remove -1 from cell_id in the df that will be used as metadata in seurat object
#bcr.df.seurat$cell_id <- gsub("-1", "", bcr.df.seurat$cell_id)



#Memory optimization 
rm(heavy_clones_clusters, light_clones_clusters, heavy_clusters, light_clusters, merged.bcr.df )
gc()
