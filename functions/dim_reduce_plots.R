#dfComposition = dframe_allASVs$ASV_composition
#install.packages('lsa')
dim_reduce_plots <- function(dfComposition,colId='name',perplexityTsne=250,neighUmap=50){
  ## objects that we are going to use
  colnames_vars = dfComposition %>% select(-any_of(colId)) %>% colnames()
  aitDistHere = robCompositions::aDist(dfComposition %>% select(-any_of(colId)))
  clrTransf = robCompositions::cenLR(dfComposition %>% select(-any_of(colId)))
  brayDist = vegan::vegdist(dfComposition %>% select(-any_of(colId)),method = 'bray')
  # cosineSimm <- lsa::cosine(x = dfComposition %>% select(-any_of(colId)) %>% as.matrix() %>% t())
  # cosineSimmCLR <- lsa::cosine(x = clrTransf$x.clr %>% as.matrix() %>% t())
  # cosineDist = 1 - cosineSimm
  # cosineDistCLR = 1 - cosineSimmCLR
  #setting the output
  out<-list()
  
  
  ###### FOR Aitchison ----------------------------------------------------
  ##### MDS----------------------------------------------------------------
  mds_obj <- cmdscale(d = as.dist(aitDistHere),
                      k = 3,eig = T,add=T)
  
  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
  
  dfComposition <- dfComposition %>% mutate(
    MDS1 = mds_obj$points[,1],
    MDS2 = mds_obj$points[,2],
    MDS3 = mds_obj$points[,3]
  )
  
  out$MDS2d_ASVs_Ait = ggplot(
    data = dfComposition,
    mapping = aes(x=MDS1,y=MDS2))+
    geom_point(alpha=0.7,size=3) +
    theme_minimal(base_size = 12) +
    xlab(paste('MDS1 -',pct_explained[1],'%'))+
    ylab(paste('MDS2 -',pct_explained[2],'%'))+
    theme(legend.position = 'bottom')+
    ggtitle('AitDist ASV compositions')
  
  ############# ------------------------------------------------------------
  library(Rtsne)
  set.seed(1234)
  tsneCoords = Rtsne(X = as.dist(aitDistHere),is_distance = TRUE, dims = 2,
                     perplexity=perplexityTsne, max_iter = 2000,learning=100)
  dfComposition$tsne1 = tsneCoords$Y[,1]
  dfComposition$tsne2 = tsneCoords$Y[,2]
  
  
  out$TSNE_ASVs_Ait_2d <- ggplot(
    data = dfComposition,
    mapping = aes(x=tsne1,y=tsne2))+
    geom_point(alpha=0.5,size=3) +
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('AitDist ASV compositions - tsne')
  
  #######---------------------------------------------------------------------
  ## umap
  library(umap)
  custom.config <- umap.defaults
  custom.config$n_neighbors=neighUmap
  custom.config$min_dist=0.05
  custom.config$spread = 0.5
  custom.config$random_state=10
  
  set.seed(1234)
  umap_obj <- umap(d = clrTransf$x.clr,config = custom.config)
  
  dfComposition = dfComposition %>%
    mutate(umap1=umap_obj$layout[,1],umap2=umap_obj$layout[,2])
  
  out$umap_ASVs_Ait_2d <- ggplot(
    data = dfComposition,
    mapping = aes(x=umap1,y=umap2))+
    geom_point(alpha=0.5,size=3) +
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('AitDist ASV compositions - umap - n_neighbors=50 - min_dist=0.05 - spread = 0.5')
  
  
  ###### FOR Bray ----------------------------------------------------
  ##### MDS----------------------------------------------------------------
  mds_obj <- cmdscale(d = as.dist(brayDist),
                      k = 3,eig = T,add=T)
  
  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
  
  dfComposition <- dfComposition %>% mutate(
    MDS1 = mds_obj$points[,1],
    MDS2 = mds_obj$points[,2],
    MDS3 = mds_obj$points[,3]
  )
  
  out$MDS2d_ASVs_Bray = ggplot(
    data = dfComposition,
    mapping = aes(x=MDS1,y=MDS2))+
    geom_point(alpha=0.7,size=3) +
    theme_minimal(base_size = 12) +
    xlab(paste('MDS1 -',pct_explained[1],'%'))+
    ylab(paste('MDS2 -',pct_explained[2],'%'))+
    theme(legend.position = 'bottom')+
    ggtitle('Bray ASV compositions')
  
  ############# ------------------------------------------------------------
  library(Rtsne)
  set.seed(1234)
  tsneCoords = Rtsne(X = as.dist(brayDist),is_distance = TRUE, dims = 2,
                     perplexity=perplexityTsne, max_iter = 2000,learning=100)
  dfComposition$tsne1 = tsneCoords$Y[,1]
  dfComposition$tsne2 = tsneCoords$Y[,2]
  
  
  out$TSNE_ASVs_Bray_2d <- ggplot(
    data = dfComposition,
    mapping = aes(x=tsne1,y=tsne2))+
    geom_point(alpha=0.5,size=3) +
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('Bray ASV compositions - tsne')
  
  
  ###  ###### FOR Cosine ----------------------------------------------------
  ###  ##### MDS----------------------------------------------------------------
  ###  mds_obj <- cmdscale(d = as.dist(cosineDistCLR),
  ###                      k = 3,eig = T,add=T)
  ###  
  ###  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
  ###  
  ###  dfComposition <- dfComposition %>% mutate(
  ###    MDS1 = mds_obj$points[,1],
  ###    MDS2 = mds_obj$points[,2],
  ###    MDS3 = mds_obj$points[,3]
  ###  )
  ###  
  ###  out$MDS2d_ASVs_Cosine = ggplot(
  ###    data = dfComposition,
  ###    mapping = aes(x=MDS1,y=MDS2))+
  ###    geom_point(alpha=0.7,size=3) +
  ###    theme_minimal(base_size = 12) +
  ###    xlab(paste('MDS1 -',pct_explained[1],'%'))+
  ###    ylab(paste('MDS2 -',pct_explained[2],'%'))+
  ###    theme(legend.position = 'bottom')+
  ###    ggtitle('Bray ASV compositions')
  ###  
  ###  ############# ------------------------------------------------------------
  ###  library(Rtsne)
  ###  set.seed(1234)
  ###  tsneCoords = Rtsne(X = as.dist(cosineDistCLR),is_distance = TRUE, dims = 2,
  ###                     perplexity=perplexityTsne, max_iter = 2000,learning=100)
  ###  dfComposition$tsne1 = tsneCoords$Y[,1]
  ###  dfComposition$tsne2 = tsneCoords$Y[,2]
  ###  
  ###  
  ###  out$TSNE_ASVs_Cosine_2d <- ggplot(
  ###    data = dfComposition,
  ###    mapping = aes(x=tsne1,y=tsne2))+
  ###    geom_point(alpha=0.5,size=3) +
  ###    theme_minimal(base_size = 12)+
  ###    theme(legend.position = 'bottom')+
  ###    ggtitle('cosineDistCLR ASV compositions - tsne')
  ###  
  ###  
  ###  ###### FOR Cosine ----------------------------------------------------
  ###  ##### MDS----------------------------------------------------------------
  ###  mds_obj <- cmdscale(d = as.dist(cosineDistCLR),
  ###                      k = 3,eig = T,add=T)
  ###  
  ###  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
  ###  
  ###  dfComposition <- dfComposition %>% mutate(
  ###    MDS1 = mds_obj$points[,1],
  ###    MDS2 = mds_obj$points[,2],
  ###    MDS3 = mds_obj$points[,3]
  ###  )
  ###  
  ###  out$MDS2d_ASVs_Cosine = ggplot(
  ###    data = dfComposition,
  ###    mapping = aes(x=MDS1,y=MDS2))+
  ###    geom_point(alpha=0.7,size=3) +
  ###    theme_minimal(base_size = 12) +
  ###    xlab(paste('MDS1 -',pct_explained[1],'%'))+
  ###    ylab(paste('MDS2 -',pct_explained[2],'%'))+
  ###    theme(legend.position = 'bottom')+
  ###    ggtitle('Bray ASV compositions')
  ###  
  ###  ############# ------------------------------------------------------------
  ###  library(Rtsne)
  ###  set.seed(1234)
  ###  tsneCoords = Rtsne(X = as.dist(cosineDistCLR),is_distance = TRUE, dims = 2,
  ###                     perplexity=perplexityTsne, max_iter = 2000,learning=100)
  ###  dfComposition$tsne1 = tsneCoords$Y[,1]
  ###  dfComposition$tsne2 = tsneCoords$Y[,2]
  ###  
  ###  
  ###  out$TSNE_ASVs_Cosine_2d <- ggplot(
  ###    data = dfComposition,
  ###    mapping = aes(x=tsne1,y=tsne2))+
  ###    geom_point(alpha=0.5,size=3) +
  ###    theme_minimal(base_size = 12)+
  ###    theme(legend.position = 'bottom')+
  ###    ggtitle('cosineDist ASV compositions - tsne')
  
return(out)
  
}
