eval_adist_clustering <- function(compositionDF){
  ## @compositionDF is a dataframe containing the following collumns:
  ### 1 - ASVnames 
  ### 2 - Clusters
  ### 3 - p The composition.
  
  #Step 1 - Calculate the aitchison distance
  aitMatrix <- robCompositions::aDist(x = compositionDF[,-c(1,2)])
  
  #Step 2 - Find the medoid
  medoid = cluster::pam(x = as.dist(aitMatrix),k = 1)
  medoid_id <- as.numeric(medoid$id.med)
  
  dist_to_OverAllMedoid <- sum(as.matrix(aitMatrix)[medoid_id,])
  allPairWiseDist <- sum(aitMatrix)/2
  allPairWiseDist_squared <- sum(as.matrix(aitMatrix)^2)/2

  #Step 3 
  ## now splitting per clustering -----------------------
  nclusters <- length(unique(compositionDF$Cluster))
  
  dist_within_clust_i <- numeric(length = nclusters)
  medoid_clust_i <- numeric(length = nclusters)
  nobs_clust_i <- numeric(length = nclusters)
  pairwise_dist_clust_i <- numeric(length = nclusters)
  pairwise_dist_clust_i_squared <- numeric(length = nclusters)
  max_dist_within <- numeric(length = nclusters)
  max_dist2medoid_within <- numeric(length = nclusters)
  min_dist_clust_i_other <- numeric(length = nclusters)
  
  for(i in 1:nclusters){
    idx_clust_i <- which(compositionDF$Cluster==i)
    
    if(length(idx_clust_i)>=2){
      
      nobs_clust_i[i] <- length(idx_clust_i)
      aux_Adist <- as.matrix(aitMatrix)[idx_clust_i,idx_clust_i]
      pairwise_dist_clust_i[i] <- sum(aux_Adist)/2
      pairwise_dist_clust_i_squared[i] <- sum(as.matrix(aux_Adist)^2)/2
      medoid_clust_i[i] = cluster::pam(x = as.dist(aux_Adist),k = 1)$id.med
      dist_within_clust_i[i] <- sum(aux_Adist[medoid_clust_i[i],])
      max_dist_within[i] <- max(aux_Adist)
      max_dist2medoid_within[i] <- max(aux_Adist[medoid_clust_i[i],])
      min_dist_clust_i_other[i] <- min(as.matrix(aitMatrix)[idx_clust_i,-idx_clust_i])
      ## keeping the medoid
      medoid_clust_i[i] <- idx_clust_i[medoid_clust_i[i]]
      
    }else{
      nobs_clust_i[i] <- 1
      aux_Adist <- 0
      pairwise_dist_clust_i[i] <- 0
      pairwise_dist_clust_i_squared[i] <- 0
      medoid_clust_i[i] = idx_clust_i[1]
      dist_within_clust_i[i] <- 0
      max_dist_within[i] <- 0
      max_dist2medoid_within[i] <- 0
      min_dist_clust_i_other[i] <- min(as.matrix(aitMatrix)[idx_clust_i,-idx_clust_i])
    }
  }
  
  #This is the numerator, the average distance between observations and its cluster medoid.
  numerator = sum(dist_within_clust_i/nobs_clust_i)/nclusters
  
  #Now, we need to work with the denominator.
  
  denominator = (sum(as.matrix(aitMatrix)[medoid_clust_i,medoid_clust_i])/2) / choose(nclusters,2)
  
  ######## permanova :
  SST = allPairWiseDist_squared / nrow(compositionDF)
  SSW = sum(pairwise_dist_clust_i_squared/nobs_clust_i)
  SSA = SST - SSW
  Fstat_num = (SSA/(nclusters-1))
  Fstat_denon = (SSW/(nrow(compositionDF)-nclusters))
  Fstat = Fstat_num/Fstat_denon
  ######## maxdist - points 
  max_max_dist_within        = max(max_dist_within)
  avg_max_dist_within        = mean(max_dist_within)
  #min_min_dist_clust_i_other = min(min_dist_clust_i_other)
  avg_min_dist_clust_i_other = mean(min_dist_clust_i_other)
  max_min_dist_clust_i_other = max(min_dist_clust_i_other)
  
  
  output<- list()
  # outputs with no medoid reference
  output$sum_pairwise = allPairWiseDist
  output$pairwise_dist_within = pairwise_dist_clust_i
  output$sum_pairwise_dist_within = sum(pairwise_dist_clust_i)
  
  # outputs with medoid reference
  output$dist_to_overall = dist_to_OverAllMedoid
  output$mean_avg_within = numerator
  output$avg_within = dist_within_clust_i/nobs_clust_i
  output$avg_between_medoids = denominator
  output$obs_per_cluster = nobs_clust_i
  output$ratio_avg_within_between = numerator/denominator
  
  output$max_dist_within = max_dist_within
  output$max_dist2medoid_within = max_dist2medoid_within
  output$min_dist_clust_i_other = min_dist_clust_i_other
 
  output$max_max_dist_within        = max_max_dist_within       
  output$avg_max_dist_within        = avg_max_dist_within       
  #output$min_min_dist_clust_i_other = min_min_dist_clust_i_other
  output$avg_min_dist_clust_i_other = avg_min_dist_clust_i_other
  output$max_min_dist_clust_i_other = max_min_dist_clust_i_other
  
  output$SST = SST
  output$SSW = SSW
  output$SSA = SSA
  output$Fstat = Fstat
  output$Fstat_num = Fstat_num
  output$Fstat_denon = Fstat_denon
  return(output)
}

#list_from_eval_adist_clustering = list_DisSum_0.1

summarise_adist_clustering <- function(list_results_eval){
  list_results_eval %>% plyr::ldply(function(el){
    data.frame(
      nclust = length(el$avg_within),
      #max_max_min_min_ratio = el$max_max_min_min_ratio,
      Fstat = el$Fstat,
      SSA = el$SSA,
      #SST = el$SST,
      #SSW = el$SSW,
      SSW = el$SSW,
      mean_avg_within = el$mean_avg_within,
      avg_between_medoids = el$avg_between_medoids,
      ratio_avg_within_between = el$ratio_avg_within_between,
      max_max_dist_within = el$max_max_dist_within,
      max_min_dist_clust_i_other = el$max_min_dist_clust_i_other,
      avg_max_dist_within = el$avg_max_dist_within,
      avg_min_dist_clust_i_other = el$avg_min_dist_clust_i_other,
      ratio_avg_max_min = el$avg_max_dist_within/el$avg_min_dist_clust_i_other
      #min_min_dist_clust_i_other = el$min_min_dist_clust_i_other
    )
  }) %>%
    return()
}


