############# Create and Evaluate Clusters::

create_eval <- function(dfComposition,inducedMatrixDist=NULL,alpha_=0){
  dfComposition = dfComposition %>% as.data.frame()
  #distMatrix_ = vegan::vegdist(dfComposition %>% select(-name),method = 'aitchison')
  distMatrix_ = robCompositions::aDist(dfComposition %>% select(-name))
  normAit = norm(as.matrix(distMatrix_),type="2")
  distMatrix_ = distMatrix_/normAit
  
  if(!is.null(inducedMatrixDist)){
    distMatrix_ = alpha_*as.matrix(inducedMatrixDist) + (1-alpha_)*as.matrix(distMatrix_)
  }
  
  obj_clusters <- clusterize_ward_pam(
    distMatrix = as.dist(distMatrix_),
    maxnclust = 50,dFrameAsv = dfComposition %>% select(name))
  
  nclusters <- ncol(obj_clusters$dFrameAsv_hclust)-1
  list_results_hclust <- list()
  list_results_medoid <- list()
  for( i in 1:nclusters){
    list_results_hclust[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      obj_clusters$dFrameAsv_hclust[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dfComposition))
    
    list_results_medoid[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      obj_clusters$dFrameAsv_pam[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dfComposition))
    
    cat('iteration------:',i,'\n')
  }
  return(list(list_results_hclust=list_results_hclust,list_results_medoid=list_results_medoid))
}


