########------------------------------------------------------
# distMatrix = hammDistMatrix

distMatrix_dimReduction <- function(distMatrix,perplexityTsne=250,neighUmap=50){
  
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
  
  #out$MDS2d = ggplot(
  #  data = output_df,
  #  mapping = aes(x=MDS1,y=MDS2))+
  #  geom_point(alpha=0.7,size=3) +
  #  theme_minimal(base_size = 12) +
  #  xlab(paste('MDS1 -',pct_explained[1],'%'))+
  #  ylab(paste('MDS2 -',pct_explained[2],'%'))+
  #  theme(legend.position = 'bottom')+
  #  ggtitle('MDS')
  
  ## tsne ----------------------------------------------------------------
  
  tsneCoords = Rtsne::Rtsne(X = as.dist(distMatrix),is_distance = TRUE, dims = 2,
                     perplexity=perplexityTsne, max_iter = 2000,learning=100)
  output_df$tsne1 = tsneCoords$Y[,1]
  output_df$tsne2 = tsneCoords$Y[,2]
  
  
  #out$TSNE_2d <- ggplot(
  #  data = output_df,
  #  mapping = aes(x=tsne1,y=tsne2))+
  #  geom_point(alpha=0.5,size=3) +
  #  theme_minimal(base_size = 12)+
  #  theme(legend.position = 'bottom')+
  #  ggtitle('tsne')
  
  
  ## umap ----------------------------------------------------------------
  library(umap)
  custom.config <- umap.defaults
  custom.config$n_neighbors=neighUmap
  custom.config$min_dist=0.05
  custom.config$spread = 0.5
  custom.config$random_state=10
  custom.config$input = 'dist'
  
  set.seed(1234)
  umap_obj <- umap(d = distMatrix,config = custom.config)
  
  output_df = output_df %>%
    mutate(umap1=umap_obj$layout[,1],umap2=umap_obj$layout[,2])
  
  #out$umap_2d <- ggplot(
  #  data = output_df,
  #  mapping = aes(x=umap1,y=umap2))+
  #  geom_point(alpha=0.5,size=3) +
  #  theme_minimal(base_size = 12)+
  #  theme(legend.position = 'bottom')+
  #  ggtitle(paste0('umap ; min_dist=0.05 ; spread = 0.5; n_neighbors =',neighUmap))
  #
  
  out$df_dimrec_axis<- output_df
  
  return(out)
}
