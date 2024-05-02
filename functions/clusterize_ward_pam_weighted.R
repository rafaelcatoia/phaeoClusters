clusterize_ward_pam_weighted <- function(distMatrix,maxnclust=30,dFrameAsv,weights_vec){
  ##### Hierarquical
  
  #function to get the clusters from wcKMedoids
  generate_output_sequence <- function(input_vector) {
    unique_numbers <- unique(input_vector)
    output_sequence <- rep(NA, length(input_vector))
    
    for (i in 1:length(unique_numbers)) {
      output_sequence[input_vector == unique_numbers[i]] <- i
    }
    
    return(output_sequence)
  }
  
  vet_nclust <- 1:maxnclust
  col_names <- c('ASV',ifelse(vet_nclust<10,
                              paste0('nClusters0',vet_nclust),
                              paste0('nClusters',vet_nclust)))
  
  dFrameAsv_hclust <- data.frame(
    ASV = dFrameAsv %>% unlist(),
    Clust1 = 1
  ) 
  
  ##### Pam
  dFrameAsv_pam <- data.frame(
    ASV = dFrameAsv %>% unlist(),
    Clust1 = 1
  )
  
  hclust_obj <- hclust(d = as.dist(distMatrix) ,method = 'ward.D',members = weights_vec)
  
  for(i in 2:maxnclust){
  obj = WeightedCluster::wcKMedoids(diss = as.dist(distMatrix),k = i,weights = weights_vec,method = 'KMedoids',cluster.only = T)
  
    dFrameAsv_hclust = cbind.data.frame(dFrameAsv_hclust,cutree(hclust_obj,k = i))
    dFrameAsv_pam = cbind.data.frame(dFrameAsv_pam,generate_output_sequence(obj))
  }
  
  rownames(dFrameAsv_hclust) <- NULL
  colnames(dFrameAsv_hclust) = col_names
  
  rownames(dFrameAsv_pam) <- NULL
  colnames(dFrameAsv_pam) <- col_names
  
  output <- list()
  
  output$dFrameAsv_hclust = dFrameAsv_hclust
  output$dFrameAsv_pam = dFrameAsv_pam
  
  return(output)
  
}
