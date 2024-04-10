clusterize_ward_pam <- function(distMatrix,maxnclust=30,dFrameAsv){
  ##### Hierarquical
  vet_nclust <- 1:maxnclust
  col_names <- c('ASV',ifelse(vet_nclust<10,
                              paste0('nClusters0',vet_nclust),
                              paste0('nClusters',vet_nclust)))
  
  dFrameAsv_hclust <- data.frame(
    ASV = dFrameAsv %>% unlist()
  ) 
  
  ##### Pam
  dFrameAsv_pam <- data.frame(
    ASV = dFrameAsv %>% unlist()
  )
  
  hclust_obj <- hclust(d = as.dist(distMatrix) ,method = 'ward.D')
  
  for(i in 1:maxnclust){
    dFrameAsv_hclust = cbind.data.frame(dFrameAsv_hclust,cutree(hclust_obj,k = i))
    dFrameAsv_pam = cbind.data.frame(dFrameAsv_pam,cluster::pam(x = as.dist(distMatrix),k = i, cluster.only=TRUE))
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
