eval_clusters_distMatrix <- function(distMatrix,df_ASV_Clusters){
  # colnames(df_ASV_Clusters) must be ASV and Clusters
  
  #distMatrix = normalizeDistMatrix(hammingDistMatrix)
  #df_ASV_Clusters = hamming_clusters_alpha0.25$dFrameAsv_hclust %>% select(ASV,nClusters22) %>% rename('Clusters'=2)
  
  df_ASV_Clusters = df_ASV_Clusters %>% rename('Clusters'=2)
  numberOfClusters = max(df_ASV_Clusters$Clusters)
  vet_number_obs_per_clusters <- df_ASV_Clusters %>% group_by(Clusters) %>% summarise(n()) %>% pull()
  idx_list <- list()
  distMatrix_subseted <- list()
  medoid_idx = rep(NA,length = numberOfClusters)
  medoid_idx_within = rep(NA,length = numberOfClusters)
  
  for(i in 1:numberOfClusters){
    # getting the index of each observation for each clusters
    idx_clust_i = which(df_ASV_Clusters$Clusters==i)
    
    # getting the index of each observation for each clusters --- storing it
    idx_list[[i]] <- idx_clust_i
    
    # subsetting the distance matrix
    distMatrix_subseted[[i]] = distMatrix[idx_clust_i,idx_clust_i]
    
    #getting the medoid id for each cluster 
    #first we store the subsetted idx for the medoid
    if(length(idx_clust_i)<2){
      medoid_idx_within[i] = 1
    }else{
      medoid_idx_within[i] = cluster::pam(x = as.dist(distMatrix_subseted[[i]]),k = 1)$id.med
    }
    
    #than we store the idx in the not subsetted version.
    medoid_idx[i] = idx_clust_i[medoid_idx_within[i]]
   # cat('\n Done for cluster ', i,'---------------------- \n')
  }
  
  ## Sum within pairwise dist
  sum_within_pairwise_dist = lapply(distMatrix_subseted,function(x) {sum(x)/2}) %>% unlist()
  
  ## Metrics within -------------------------------------------------------------------------
  sum_within_dist2medoid <- rep(NA,length = numberOfClusters)
  avg_within_dist2medoid <- rep(NA,length = numberOfClusters)
  max_within_dist2medoid <- rep(NA,length = numberOfClusters)
  max_within_dist <- rep(NA,length = numberOfClusters)
  min_within_dist <- rep(NA,length = numberOfClusters)
  for(i in 1:numberOfClusters){
    if(is.null(dim(distMatrix_subseted[[i]]))){
      #value if the cluster is composed by only one observation
      sum_within_dist2medoid[i] = 0
      avg_within_dist2medoid[i] = 0
      max_within_dist2medoid[i] = 0
      max_within_dist[i]=0
      min_within_dist[i]=0
    }else{
      sum_within_dist2medoid[i] = sum(distMatrix_subseted[[i]][medoid_idx_within[i],])
      avg_within_dist2medoid[i] = sum_within_dist2medoid[i]/(nrow(distMatrix_subseted[[i]])-1)
      max_within_dist2medoid[i] = max(distMatrix_subseted[[i]][medoid_idx_within[i],])
      max_within_dist[i]=max(distMatrix_subseted[[i]])
      min_within_dist[i]=min(distMatrix_subseted[[i]][upper.tri(distMatrix_subseted[[i]])])
    }
  }
  
  ## Metrics between -------------------------------------------------------------------------
  combinations_between_clusters = choose(numberOfClusters,2)
  sum_dist_between_medoid <- rep(NA,length = numberOfClusters)
  avg_dist_between_medoid <- rep(NA,length = numberOfClusters)
  max_dist_between_medoid <- rep(NA,length = numberOfClusters)
  min_dist_between_medoid <- rep(NA,length = numberOfClusters)
  max_dist_between_dist <- rep(NA,length = numberOfClusters)
  min_dist_between_dist <- rep(NA,length = numberOfClusters)
  
  ### Calculating  ---------------------------------------------------------------------------
  dist_matrix_medoids = distMatrix[medoid_idx,medoid_idx]
  
  ## sum of distances between medoids
  sum_dist_between_medoid = sum(dist_matrix_medoids)/2
  ## avg of distances between medoids
  avg_dist_between_medoid = sum_dist_between_medoid/combinations_between_clusters
  
  ## max and min between medoids
  max_dist_between_medoid = max(dist_matrix_medoids)
  min_dist_between_medoid = min(dist_matrix_medoids[upper.tri(dist_matrix_medoids)])
  
  ## max between clusters and min between clusters
  min_dist_between_clusters = matrix(0,nrow = numberOfClusters,ncol=numberOfClusters)
  for(i in 1:numberOfClusters){
    for(j in 1:numberOfClusters){
      if(i>j){
        dist_aux = distMatrix[idx_list[[i]],idx_list[[j]]]
        min_dist_between_clusters[i,j] <- min(dist_aux)
      }
    }
  }
  
  min_dist_between_clusters=min_dist_between_clusters[lower.tri(min_dist_between_clusters)]
  avg_min_dist_between_clusters = mean(min_dist_between_clusters)
  sum_min_dist_between_clusters = sum(min_dist_between_clusters)
  
  ## arranging output
  out <- list()
  ## avg_within2medoid / avg_between_medoid
  weighted_mean_avg_within_to_medoid = 
    sum(avg_within_dist2medoid * (vet_number_obs_per_clusters/
                                    sum(vet_number_obs_per_clusters)))
  
  ## to medoid
  out$vet_number_obs_per_clusters=vet_number_obs_per_clusters
  out$vet_avg_within_to_medoid = avg_within_dist2medoid
  out$mean_avg_within_to_medoid = mean(avg_within_dist2medoid)
  out$weighted_mean_avg_within_to_medoid = weighted_mean_avg_within_to_medoid
  out$avg_dist_between_medoid = avg_dist_between_medoid
  out$ratio_avg_within_between_medoids = weighted_mean_avg_within_to_medoid/avg_dist_between_medoid
  
  ## within
  out$max_max_dist_within = max(max_within_dist)
  out$sum_weighted_max_dist_within = sum(max_within_dist*(vet_number_obs_per_clusters/sum(vet_number_obs_per_clusters)))
  out$min_max_dist_within = min(max_within_dist)
  out$avg_max_dist_within = mean(max_within_dist)
  
  ## min between clusters
  out$min_dist_between_clusters     = min_dist_between_clusters
  out$avg_min_dist_between_clusters = avg_min_dist_between_clusters
  out$sum_min_dist_between_clusters = sum_min_dist_between_clusters
  
  ## Previous measuremens
  
  out$within_between_avg = (sum(avg_within_dist2medoid)/numberOfClusters) / avg_dist_between_medoid
  out$numer_within_between_avg = (sum(avg_within_dist2medoid)/numberOfClusters)
  out$denom_within_between_avg = avg_dist_between_medoid
  
  return(out)
}
