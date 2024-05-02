########------------------------------------------------------
# distMatrix = hammDistMatrix

distMatrix_dimReduc_coords <- function(distMatrix,perplexityTsne=250,neighUmap=50){
  
  out <- list()
  ## MDS ----------------------------------------------------------------
  mds_obj <- cmdscale(d = as.dist(distMatrix),
                      k = 3,eig = T,add=T)
  
  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
  
  output_df = data.frame(
    MDS1= mds_obj$points[,1],
    MDS2= mds_obj$points[,2],
    MDS3= mds_obj$points[,3]
  )
  
  
  
  # nMDS ----------------------------------------------------------------
  # For next steps
  non_metricMDS = vegan::metaMDS(comm = as.dist(distMatrix),k = 4)
  
  output_df = output_df %>% 
    mutate(
      nMDS1=non_metricMDS$points[,1],
      nMDS2=non_metricMDS$points[,2],
      nMDS3=non_metricMDS$points[,3])
  
  ## tsne ----------------------------------------------------------------
  
  tsneCoords = Rtsne::Rtsne(X = as.dist(distMatrix),is_distance = TRUE, dims = 2,
                     perplexity=perplexityTsne, max_iter = 2000,learning=100)
  
  output_df$tsne1 = tsneCoords$Y[,1]
  output_df$tsne2 = tsneCoords$Y[,2]
  
  ## umap ----------------------------------------------------------------
  library(umap)
  custom.config <- umap.defaults
  custom.config$n_neighbors=neighUmap
  custom.config$min_dist=0.05
  custom.config$spread = 0.5
  custom.config$random_state=10
  custom.config$input = 'dist'
  
  set.seed(1234)
  umap_obj <- umap(d = as.matrix(distMatrix),config = custom.config)
  
  output_df = output_df %>%
    mutate(umap1=umap_obj$layout[,1],umap2=umap_obj$layout[,2])
  
  return(output_df)
}
