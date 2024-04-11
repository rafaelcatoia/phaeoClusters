getwd()
setwd('/Users/rafaelcatoia/Desktop/repos/clustersPhaeo/phaeoClusters/saved_files/')
list_with_clusters_alpha0 <- readRDS('list_with_clusters_alpha0')
list_with_clusters_alpha0$dFrameAsv_hclust


randomize_points <- function(){
  
}
# 
# fit_cluster <- function(){
#   
# }
# 
# df_clusters <- list_with_clusters_alpha0$dFrameAsv_hclust[,c(1,4)]
# df_clusters$cluster_id <- df_clusters$nClusters03
# 
# dist_given_cluster <- function(df_clusters){
#   out <- as.matrix(dist(df_clusters$cluster_id))
#   out <- ifelse(out == 0, 0, 1)
#   rownames(out) <- df_clusters$ASV
#   colnames(out) <- df_clusters$ASV
#   return(out)
# }
# 
# dist_clusters <- array(NA, c(nrow(df), ncol(df), B))
# for(i in 1:B){
#   df_ <- randomize_points(df)
#   df_clusters <- fit_cluster(df_, k)
#   dist_clusters[,,i] <- dist_given_cluster(df_clusters)
# }
# stability_matrix <- apply(dist_clusters, 3, mean)







