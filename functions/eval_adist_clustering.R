######## -----------------------------------
## Evaluating fit for compositioins --------
######## -----------------------------------

# #Generating a mock-data-set only to verify
# set.seed(1234)
# testingData <- rbind(
#   dirmult::rdirichlet(n=,5,alpha=c(0.5,1,1,5)),
#   dirmult::rdirichlet(n=,5,alpha=c(5,1,1,0.5)),
#   dirmult::rdirichlet(n=,5,alpha=c(0.5,5,5,0.5))
#   )
# 
# 
# testingData <- data.frame(testingData)
# colnames(testingData) <- paste0('Samp',1:4)
# 
# testingData <- data.frame(
#   name=paste0('name',1:15),
#   Cluster = c(rep(1,5),rep(2,5),rep(3,5)),
#   testingData)
# 
# compositionDF <- testingData

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
  
  #Step 3 
  ## now splitting per clustering -----------------------
  nclusters <- length(unique(compositionDF$Cluster))
  
  dist_within_clust_i <- numeric(length = nclusters)
  medoid_clust_i <- numeric(length = nclusters)
  nobs_clust_i <- numeric(length = nclusters)
  pairwise_dist_clust_i <- numeric(length = nclusters)
  
  for(i in 1:nclusters){
    idx_clust_i <- which(compositionDF$Cluster==i)
    
    if(length(idx_clust_i)>=2){
      
      nobs_clust_i[i] <- length(idx_clust_i)
      aux_Adist <- as.matrix(aitMatrix)[idx_clust_i,idx_clust_i]
      pairwise_dist_clust_i[i] <- sum(aux_Adist)/2
      medoid_clust_i[i] = cluster::pam(x = as.dist(aux_Adist),k = 1)$id.med
      dist_within_clust_i[i] <- sum(aux_Adist[medoid_clust_i[i],])
      medoid_clust_i[i] <- idx_clust_i[medoid_clust_i[i]]
      
    }else{
      nobs_clust_i[i] <- 1
      aux_Adist <- 0
      pairwise_dist_clust_i[i] <- 0
      medoid_clust_i[i] = idx_clust_i[1]
      dist_within_clust_i[i] <- 0
      medoid_clust_i[i] <- idx_clust_i[1]
    }
  }
  
  #This is the numerator, the average distance between observations and its cluster medoid.
  numerator = sum(dist_within_clust_i/nobs_clust_i)/nclusters
  
  #Now, we need to work with the denominator.
  
  denominator = (sum(as.matrix(aitMatrix)[medoid_clust_i,medoid_clust_i])/2) / choose(nclusters,2)
  
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
  return(output)
}

#list_from_eval_adist_clustering = list_DisSum_0.1

#nested_eval_list = list_DisSum_0.1
summarise_AVG_Dist <- function(nested_eval_list,name_alpha){
  
  list2df <- function(eval_list){
    return(df_summary = plyr::ldply(eval_list, function(el){
      data.frame(
        ratio_avg_within_between = el$ratio_avg_within_between, 
        mean_avg_within = el$mean_avg_within,
        avg_between_medoids = el$avg_between_medoids
      )},
      .id = "test") %>% 
        mutate(nclusters=as.integer( 1:n()+1)))
  }
  vetnames = names(nested_eval_list)
  
  df_list <- list()
  for( i in 1:length(nested_eval_list)){
    df_list[[i]] <- list2df(nested_eval_list[[i]]) %>% 
      mutate(alpha_val=name_alpha,scenario=vetnames[i])
  }
  
  df_output <- data.table::rbindlist(df_list)
  
  return(df_output)
  
}


