#list_with_clusters_alpha0 = readRDS(file = paste0(savingdir,'/','list_with_clusters_alpha0'))
#list_with_clusters_alpha0.10 = readRDS(file = paste0(savingdir,'/','list_with_clusters_alpha0.10'))
#list_with_clusters_alpha0.25 = readRDS(file = paste0(savingdir,'/','list_with_clusters_alpha0.25'))
#list_with_clusters_alpha0.50 = readRDS(file = paste0(savingdir,'/','list_with_clusters_alpha0.50'))
#
#list_with_clusters = list_with_clusters_alpha0
#ait_distMatrix <- inducedDist$normalizedAitDist

permanova_custom <- function(
    list_with_clusters,
    ait_distMatrix,cpus=7,
    B=100){
  
  # hardcoded but... why not?
  number_clusters <- max(list_with_clusters_alpha0$dFrameAsv_hclust[,-1])
  list_permanova_hclust <-list()
  list_permanova_pam <-list()
  for (i in 2:number_clusters) {
    
    dfGroups <- data.frame(groups_ = as.factor(list_with_clusters$dFrameAsv_hclust[,i+1]))
    list_permanova_hclust[[i]] <- vegan::adonis2(
      formula = as.dist(ait_distMatrix)~groups_,data = dfGroups,
      permutations = B,parallel = cpus)
    list_permanova_hclust[[i]]$nclust <- i
    
  
    dfGroups <- data.frame(groups_=as.factor(list_with_clusters$dFrameAsv_pam[,i+1]))
    list_permanova_pam[[i]] <- vegan::adonis2(
      formula = as.dist(ait_distMatrix)~groups_,data = dfGroups,
      permutations = B,parallel = cpus)
    list_permanova_pam[[i]]$nclust <- i
    cat('iteration--------------- ',i,' of ', number_clusters,'\n')
  }
  return(list(hclust=list_permanova_hclust,
              pam=list_permanova_pam))
}
