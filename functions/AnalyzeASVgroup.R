
#generates all the plots and analysis.

AnalyzeASVgroup <- function(
    dfLonger, #GrumpDataset a column of ASV cluster
    clusterVar='Cluster', 
    sampleVar='SampleKey',
    abioticVars = c('Salinity','Oxygen','Temperature','Silicate','NO2','NO3','NOx','NH3')
) # Name of the key for each sample
{
  
  ## First, lets compare the ASV groupping and the Species
  ## Changing the Species to Missing
  dfLonger = dfLonger %>%  mutate(Species=ifelse(Species=='Phaeocystis_sp.',NA,Species))
  
  ## Redefining the columns that we need
  dfLonger$Cluster = dfLonger %>% select(one_of(clusterVar)) %>% pull()
  dfLonger$Cluster = as.numeric(dfLonger$Cluster)
  dfLonger = dfLonger %>% mutate(
    Cluster = factor(
      ifelse(Cluster < 10,paste0('ASVG_0',Cluster),paste0('ASVG_',Cluster))
    )
  )
  
  dfLonger$ID_Sample = dfLonger %>% select(one_of(sampleVar)) %>% pull()
  SampleKeysNames <- unique(dfLonger$ID_Sample)
  
  
  ### Setting Longhurst Colors
  longhurst_colours <- c("Polar - Boreal Polar Province (POLR)"="#fe0080",
                         "Polar - N. Pacific Epicontinental Province"="#ffcbd7",
                         "Westerlies - Pacific Subarctic Gyres Province (East)"="#6e000b",
                         "Coastal - NW Atlantic Shelves Province"="#f75a00",
                         "Westerlies - Gulf Stream Province"="#734900",
                         "Westerlies - N. Atlantic Subtropical Gyral Province (East) (STGE)"="#ffc274",
                         "Westerlies - N. Pacific Polar Front Province"="#ffe1a4",
                         "Westerlies - N. Atlantic Subtropical Gyral Province (West) (STGW)"="#b49200",
                         "Trades - N. Pacific Tropical Gyre Province"="#cebd00",
                         "Trades - N. Atlantic Tropical Gyral Province (TRPG)"="#688900",
                         "Coastal - E. India Coastal Province"="#8ffd32",
                         "Trades - Western Tropical Atlantic Province"="#cfffa6",
                         "Trades - N. Pacific Equatorial Countercurrent Province"="#0d1c00",
                         "Trades - Pacific Equatorial Divergence Province"="#005d07",
                         "Trades - South Atlantic Gyral Province (SATG)"="#7fffae",
                         "Westerlies - S. Pacific Subtropical Gyre Province"="#00af64",
                         "Trades - Indian Monsoon Gyres Province"="#007859",
                         "Trades - Archipelagic Deep Basins Province"="#c5fcff",
                         "Coastal - East Australian Coastal Province"="#004a4f",
                         "Coastal - E. Africa Coastal Province"="#01b0be", 
                         "Coastal - Benguela Current Coastal Province"="#0058b3",
                         "Trades - Indian S. Subtropical Gyre Province"="#b8aeff",
                         "Westerlies - S. Subtropical Convergence Province"="#21003f", 
                         "Coastal - New Zealand Coastal Province"="#eac7ff",
                         "Coastal - SW Atlantic Shelves Province"="#7b00a5",
                         "Westerlies - Subantarctic Province"="#d95dff",
                         "Polar - Antarctic Province"="#630037",
                         "Polar - Austral Polar Province"="#ff8fb9")
  
  longhurst_order <- rev(c("Polar - Boreal Polar Province (POLR)",
                           "Polar - N. Pacific Epicontinental Province",
                           "Westerlies - Pacific Subarctic Gyres Province (East)",
                           "Coastal - NW Atlantic Shelves Province",
                           "Westerlies - Gulf Stream Province",
                           "Westerlies - N. Atlantic Subtropical Gyral Province (East) (STGE)",
                           "Westerlies - N. Pacific Polar Front Province",
                           "Westerlies - N. Atlantic Subtropical Gyral Province (West) (STGW)",
                           "Trades - N. Pacific Tropical Gyre Province",
                           "Trades - N. Atlantic Tropical Gyral Province (TRPG)",
                           "Coastal - E. India Coastal Province",
                           "Trades - Western Tropical Atlantic Province",
                           "Trades - N. Pacific Equatorial Countercurrent Province",
                           "Trades - Pacific Equatorial Divergence Province",
                           "Trades - South Atlantic Gyral Province (SATG)",
                           "Westerlies - S. Pacific Subtropical Gyre Province",
                           "Trades - Indian Monsoon Gyres Province",
                           "Trades - Archipelagic Deep Basins Province",
                           "Coastal - East Australian Coastal Province",
                           "Coastal - E. Africa Coastal Province",
                           "Coastal - Benguela Current Coastal Province",
                           "Trades - Indian S. Subtropical Gyre Province",
                           "Westerlies - S. Subtropical Convergence Province",
                           "Coastal - New Zealand Coastal Province",
                           "Coastal - SW Atlantic Shelves Province",
                           "Westerlies - Subantarctic Province",
                           "Polar - Antarctic Province",
                           "Polar - Austral Polar Province"))
  
  dfLonger$Longhurst_Long <- factor(
    dfLonger$Longhurst_Long,longhurst_order)
  
  ### Now gettin the colors for the samples
  
  colorASVgs<-viridis::turbo(n = length(unique(dfLonger$Cluster)))
  
  ## output
  out<-list()
  
  ## Checkings:
  cat(paste0('In this dataset we have: ',length(unique(dfLonger$ASV_name)),' ASVs, and ',
             length(unique(dfLonger$ID_Sample)),' samples \n'))
  
  ### Creating the freq df
  df_freq = dfLonger %>% select(ASV_name) %>% distinct() %>% 
    left_join(dfLonger %>%select(ASV_name,Species,Cluster) %>% distinct()) %>% 
    group_by(Species,Cluster) %>% 
    summarise(Freq=n()) %>% 
    ungroup() %>% group_by(Cluster) %>%  
    mutate(Pct=Freq/sum(Freq),Cluster=as.factor(Cluster)) %>% 
    arrange(Cluster,Species)
  
  out$barplot1 <- df_freq %>% ggplot(aes(x=Cluster,y=Freq,fill=Species))+
    geom_col()+#coord_flip()+
    theme_minimal(base_size = 14)+
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  out$barplot2 <- df_freq %>% ggplot(aes(x=Cluster,y=Pct,fill=Species))+
    geom_col()+#coord_flip()+
    theme_minimal(base_size = 14) +
    theme(legend.position = 'bottom',
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  out$Freq <- df_freq %>% pivot_wider(
    id_cols = Species,
    names_from = Cluster,
    values_from = Freq,
    values_fill = 0) %>% 
    knitr::kable() %>% kableExtra::kable_styling()
  
  out$RelativeFreq <- df_freq %>% pivot_wider(
    id_cols = Species,
    names_from = Cluster,
    values_from = Pct,values_fill = 0)  %>% 
    knitr::kable() %>% kableExtra::kable_styling()
  
  ## First step is to aggregate by the ASV groupping
  ## Here we have the stacked version of the RA for all ASV Groupping that we created.
  dfBiotic = dfLonger %>%
    select(ID_Sample,ASV_name,Raw.Sequence.Counts) %>% 
    pivot_wider(id_cols = c('ASV_name'),
                values_from ='Raw.Sequence.Counts',
                names_from='ID_Sample',values_fill = 0) %>%
    select(ASV_name,sort(names(.))) %>% 
    left_join(dfLonger %>% select(ASV_name,Cluster) %>% distinct()) %>% 
    group_by(Cluster) %>% 
    summarise(across(where(is.numeric),sum)) %>% 
    pivot_longer(cols = all_of(SampleKeysNames),
                 names_to = 'ID_Sample',
                 values_to = 'RawSeqCounts') %>%
    group_by(ID_Sample) %>% 
    transmute(ID_Sample=ID_Sample,
              Cluster=Cluster,
              RA=RawSeqCounts/sum(RawSeqCounts),
              RawSeqCounts=RawSeqCounts,
              TotalRawSeqCounts = sum(RawSeqCounts)) %>% 
    arrange(ID_Sample,Cluster) %>% ungroup
  
  ## Lets now get the Most Abundant Group for each sample
  
  dfBioticMostAbundant <- dfBiotic %>% 
    group_by(ID_Sample) %>%
    slice_max(RA,with_ties = FALSE) %>% ungroup
  
  out$MostAbundantFreq = dfBioticMostAbundant %>% 
    group_by(Cluster) %>% 
    summarise(Freq=n(),Pct=n()/nrow(dfBioticMostAbundant))
  
  ###########################################################
  # Now, lets introduce some abiotic factors
  ###########################################################
  
  dfAbioMostAbundant <- dfBioticMostAbundant %>% 
    left_join(dfLonger %>% select(ID_Sample,Latitude,Longitude,Depth,
                                  Cruise,Salinity,Oxygen,Temperature,
                                  Silicate,NO2,NO3,NOx,NH3) %>% distinct()) %>% 
    mutate(LogDepth=log(Depth+1))
  
  #######################
  ## Map vizualizations
  #######################
  out$mapView = mapview::mapview(dfAbioMostAbundant,
                                 xcol = "Longitude",
                                 ycol = "Latitude", crs = 4269,
                                 zcol = 'Cluster',
                                 col.regions = colorASVgs,
                                 legend = TRUE,
                                 grid = F)
  
  
  out$MissingAbiotic <- colSums(is.na(dfAbioMostAbundant))/nrow(dfAbioMostAbundant)
  
  
  ######################################################
  # Matrix plot for each abiotic var
  ######################################################
  out$MatrixAbioticsUnivar = dfAbioMostAbundant %>% 
    filter(Oxygen>-10) %>% 
    pivot_longer(cols = all_of(abioticVars),
                 names_to = 'AbioticVar') %>% 
    filter(!is.na(value)) %>% 
    ggplot(aes(y=RA,x=value,color=Cluster))+
    geom_point(alpha=0.50)+
    facet_grid(AbioticVar~Cluster,scales = 'free_y')+
    theme_minimal()+
    theme(legend.position = 'bottom')+
    coord_flip() +
    scale_fill_manual(values = colorASVgs)+
    scale_color_manual(values = colorASVgs)
  
  
  gpairs_lower <- function(g){
    g$plots <- g$plots[-(1:g$nrow)]
    g$yAxisLabels <- g$yAxisLabels[-1]
    g$nrow <- g$nrow -1
    
    g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
    g$xAxisLabels <- g$xAxisLabels[-g$ncol]
    g$ncol <- g$ncol - 1
    
    g
  }
  
  out$MatrixAbioticsBiVar <- gpairs_lower(GGally::ggpairs(
    data = dfAbioMostAbundant %>% filter(Oxygen>-10) %>%
      select(Cluster,RA,all_of(abioticVars)),
    
    aes(color=Cluster,alpha=0.10),
    upper = NULL,diag=NULL,columns = 3:10,
    lower=list(continuous="points", combo="facethist", discrete="facetbar", na="points")
  )+
    theme_minimal()+scale_color_manual(values = colorASVgs))
  
  ######### Lat long ---------------------------
  
  out$LatDepth = dfAbioMostAbundant %>% 
    select(Latitude,LogDepth,Longitude,Cluster,RA,all_of(abioticVars)) %>% 
    ggplot(aes(y=LogDepth,Latitude,size=RA,color=Cluster))+
    geom_point(alpha=0.5) + 
    scale_y_reverse()+
    theme_minimal()+
    theme(legend.position = 'bottom')+
    scale_color_manual(values = colorASVgs)+
    scale_fill_manual(values = colorASVgs)
  
  out$LatLong = dfAbioMostAbundant %>% ungroup %>%
    select(Latitude,LogDepth,Longitude,Cluster,RA,all_of(abioticVars)) %>% 
    ggplot(aes(y=Latitude,Longitude,size=RA,color=Cluster))+
    geom_point(alpha=0.5) + 
    theme_minimal()+
    theme(legend.position = 'bottom')+
    scale_color_manual(values = colorASVgs)+
    scale_fill_manual(values = colorASVgs)
  
  out$DepthLong=dfAbioMostAbundant %>% 
    select(Latitude,LogDepth,Longitude,Cluster,RA,all_of(abioticVars)) %>% 
    ggplot(aes(y=LogDepth,Longitude,size=RA,color=Cluster))+
    geom_point(alpha=0.5) + 
    scale_y_reverse()+
    theme_minimal()+
    theme(legend.position = 'bottom')+
    scale_color_manual(values = colorASVgs)+
    scale_fill_manual(values = colorASVgs)
  
  ############ ------- LongHurst 
  out$FreqTable_LongHurst_ASVG = dfAbioMostAbundant %>% 
    select(ID_Sample,Latitude,LogDepth,Longitude,Cluster,RA) %>% 
    left_join(dfLonger %>% select(ID_Sample,Longhurst_Long) %>% distinct()) %>%
    group_by(Longhurst_Long) %>% mutate(FreqLongHust=n()) %>% ungroup %>% 
    group_by(Longhurst_Long,Cluster) %>%
    mutate(PctASVG_LH = n()/FreqLongHust) %>% ungroup %>% 
    select(Longhurst_Long,Cluster, FreqLongHust,PctASVG_LH) %>% distinct() %>% 
    arrange(Longhurst_Long,Cluster)
  
  
  out$BarPlotLongHurst_ASVG <- out$FreqTable_LongHurst_ASVG %>% 
    ggplot(aes(y=Longhurst_Long,fill=Cluster,x=PctASVG_LH))+
    geom_bar(stat='identity')+
    theme_minimal()+
    theme(legend.position = 'bottom')+
    scale_color_manual(values = colorASVgs)+
    scale_fill_manual(values = colorASVgs)
  
  out$FreqTable_ASVG_LH = dfAbioMostAbundant %>% 
    select(ID_Sample,Latitude,LogDepth,Longitude,Cluster,RA) %>% 
    left_join(dfLonger %>% select(ID_Sample,Longhurst_Long) %>% distinct()) %>%
    group_by(Cluster) %>% mutate(Freq=n()) %>% ungroup %>% 
    group_by(Cluster,Longhurst_Long) %>%
    mutate(PctASVG_LH = n()/Freq) %>% ungroup %>% 
    select(Cluster,Longhurst_Long, Freq,PctASVG_LH) %>% distinct() %>% 
    arrange(Cluster,Longhurst_Long)
  
  out$BarPlot_ASVG_LH <- out$FreqTable_ASVG_LH %>% 
    ggplot(aes(y=Cluster,fill=Longhurst_Long,x=PctASVG_LH))+
    geom_bar(stat='identity')+
    theme_minimal()+
    theme(legend.position = 'bottom')+
    scale_fill_manual(values = longhurst_colours)
  
  
  
  ######### Relative abundance per Longhurst 
  out$BarPlot_ASVG_LH_RA = dfLonger %>%
    group_by(Longhurst_Long) %>% 
    mutate(TotalSeq=sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    group_by(Longhurst_Long,Cluster) %>% 
    mutate(SeqPerASVG=sum(Raw.Sequence.Counts)) %>% 
    select(Longhurst_Long,Cluster,TotalSeq,SeqPerASVG) %>% distinct() %>% 
    mutate(PctASVG = SeqPerASVG/TotalSeq) %>% 
    ggplot(aes(y=Longhurst_Long,fill=Cluster,x=PctASVG))+
    geom_bar(stat='identity')+
    scale_fill_manual(values = colorASVgs)+
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('RA of each ASVG inside each Longhurst')
  
  out$BarPlot_ASVG_LH_RA2 = dfLonger %>%
    group_by(Cluster) %>% 
    mutate(TotalSeq=sum(Raw.Sequence.Counts)) %>% ungroup %>% 
    group_by(Cluster,Longhurst_Long) %>% 
    mutate(SeqPerASVG=sum(Raw.Sequence.Counts)) %>% 
    select(Longhurst_Long,Cluster,TotalSeq,SeqPerASVG) %>% distinct() %>% 
    mutate(PctLH = SeqPerASVG/TotalSeq) %>% 
    ggplot(aes(y=Cluster,fill=Longhurst_Long,x=PctLH))+
    geom_bar(stat='identity')+
    scale_fill_manual(values = longhurst_colours)+
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('RA of each Longhurst inside each ASVG')
  
  ###################################
  ################################### 
  #### ASV Compositions :
  df_asv_composition = tidy_grump(dfLonger)$ASV_composition
  df_asv_composition = df_asv_composition %>% rename('ASV_name'=1) %>%
    left_join(dfLonger %>% select(ASV_name,Cluster) %>% distinct()) %>% data.frame()
  
  selectcols  = colnames(df_asv_composition)
  selectcols  = selectcols[which(selectcols%!in%c('ASV_name','Cluster'))]
  aitDistHere = robCompositions::aDist(df_asv_composition[,selectcols])
  
  mds_obj <- cmdscale(d = aitDistHere,
                      k = 3,eig = T,add=T)
  
  pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
  
  df_asv_composition <- df_asv_composition %>% mutate(
    MDS1 = mds_obj$points[,1],
    MDS2 = mds_obj$points[,2],
    MDS3 = mds_obj$points[,3]
  )
  
  out$MDS2d_ASVs_Ait = ggplot(
    data = df_asv_composition,
    mapping = aes(x=MDS1,y=MDS2,
                  color=Cluster))+
    geom_point(alpha=0.7,size=3) +
    theme_minimal(base_size = 12) +
    xlab(paste('MDS1 -',pct_explained[1],'%'))+
    ylab(paste('MDS2 -',pct_explained[2],'%'))+
    theme(legend.position = 'bottom')+
    ggtitle('AitDist ASV compositions')
    
    
  ###########################################
  #tsneplots
  library(Rtsne)
  set.seed(1234)
  tsneCoords = Rtsne(X = as.dist(aitDistHere),is_distance = TRUE, dims = 2,
                     perplexity=150, max_iter = 2000,learning=100)
  df_asv_composition$tsne1 = tsneCoords$Y[,1]
  df_asv_composition$tsne2 = tsneCoords$Y[,2]
  
  out$TSNE_ASVs_Ait_2d <- ggplot(
    data = df_asv_composition,
    mapping = aes(x=tsne1,y=tsne2,
                  color=Cluster))+
    geom_point(alpha=0.5,size=3) +
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('AitDist ASV compositions - tsne')
  
  set.seed(1234)
  tsneCoords = Rtsne(X = as.dist(aitDistHere),is_distance = TRUE, dims = 3,
                     perplexity=100, max_iter = 2000,learning=100)
  df_asv_composition$tsne1 = tsneCoords$Y[,1]
  df_asv_composition$tsne2 = tsneCoords$Y[,2]
  df_asv_composition$tsne3 = tsneCoords$Y[,3]
  df_asv_composition = df_asv_composition %>% mutate(Cluster_num=as.numeric(Cluster))
  
  axx <- list(
    title = "tsne1"
  )
  
  axy <- list(
    title = "tsne2"
  )
  
  axz <- list(
    title = "tsne3"
  )
  
  out$TSNE_ASVs_Ait_3d <- plotly::plot_ly(
    x=df_asv_composition$tsne1,
    y=df_asv_composition$tsne2,
    z=df_asv_composition$tsne3,
    color = as.numeric(df_asv_composition$Cluster_num),
    type="scatter3d",
    mode="markers") %>%
    #plotly::layout(
    #  scene = list(aspectratio = list(x = 4, y = 4, z = 3))) %>% 
    #list(xaxis=list(title='MDS1'),yaxis=list(title='MDS2'),zaxis=list(title='MDS3'))) %>% 
    plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
    plotly::hide_colorbar() 
  ###################-------------------------------------------------------------------------
  ##### Umap --------------------------
  ###################-------------------------------------------------------------------------
  library(umap)
  custom.config <- umap.defaults
  custom.config$n_neighbors=25
  custom.config$min_dist=0.05
  custom.config$spread = 0.5
  custom.config$random_state=10
  
  clr_transf <- robCompositions::cenLR(df_asv_composition[,selectcols])
  #test_duplicated = duplicated(clr_transf) there are no duplicates!
  #dim(clr_transf$x.clr)
  clr_transf = bind_cols(
    ASV_name = df_asv_composition$ASV_name,
    Cluster = df_asv_composition$Cluster,
    clr_transf$x.clr)
  
  
  set.seed(1234)
  umap_obj <- umap(d = clr_transf %>% select(-ASV_name,-Cluster),config = custom.config)
  #umap_obj <- umap(d = as.matrix(aitDistHere),input="dist")
  #plot(umap_obj$layout)
  df_asv_composition = df_asv_composition %>%
    mutate(umap1=umap_obj$layout[,1],umap2=umap_obj$layout[,2])

  
  out$umap_ASVs_Ait_2d <- ggplot(
    data = df_asv_composition,
    mapping = aes(x=umap1,y=umap2,
                  color=Cluster))+
    geom_point(alpha=0.5,size=3) +
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('AitDist ASV compositions - umap - n_neighbors=50 - min_dist=0.1 -spread = 0.5')
  
  
  custom.config$n_components=3
  set.seed(1234)
  umap_obj <- umap(d = clr_transf %>% select(-ASV_name,-Cluster),config = custom.config)
  #umap_obj <- umap(d = as.matrix(aitDistHere),input="dist")
  #plot(umap_obj$layout)
  df_asv_composition = df_asv_composition %>%
    mutate(umap1=umap_obj$layout[,1],umap2=umap_obj$layout[,2],umap3=umap_obj$layout[,3])
  
  
  axx <- list(
    title = "umap1"
  )
  
  axy <- list(
    title = "umap2"
  )
  
  axz <- list(
    title = "umap3"
  )
  
  out$umap_ASVs_Ait_3d <- plotly::plot_ly(
    x=df_asv_composition$umap1,
    y=df_asv_composition$umap2,
    z=df_asv_composition$umap3,
    color = as.numeric(df_asv_composition$Cluster_num),
    type="scatter3d",
    mode="markers") %>%
    #plotly::layout(
    #  scene = list(aspectratio = list(x = 4, y = 4, z = 3))) %>% 
    #list(xaxis=list(title='MDS1'),yaxis=list(title='MDS2'),zaxis=list(title='MDS3'))) %>% 
    plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
    plotly::hide_colorbar() 
  
  
  #######################################-------------------------------------------------------------------------------------------
  ### Now using the Sample compositions ---composition
  
  
  dfBiotic_Wider = dfBiotic %>% 
    select(ID_Sample,Cluster,RA) %>%
    pivot_wider(names_from = Cluster,values_from = RA)
  
  dfBiotic_Wider_MDS <- dfBiotic_Wider %>% 
    mutate(across(where(is.numeric),function(x) x+0.00000001)) %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) 
  
  
  aitDistHere <- robCompositions::aDist(x = dfBiotic_Wider_MDS[,-1])
  
  set.seed(1234)
  tsneCoords = Rtsne(X = as.dist(aitDistHere),is_distance = TRUE, dims = 2,
                     perplexity=100, max_iter = 2000,learning=100,)
  dfBiotic_Wider_MDS$tsne1 = tsneCoords$Y[,1]
  dfBiotic_Wider_MDS$tsne2 = tsneCoords$Y[,2]
  dfBiotic_Wider_MDS = dfBiotic_Wider_MDS %>%  left_join(
    dfAbioMostAbundant %>% 
      select(ID_Sample,Cluster,RA))
  
  out$TSNE_Sample_Ait_2d <- ggplot(
    data = dfBiotic_Wider_MDS,
    mapping = aes(x=tsne1,y=tsne2,
                  color=Cluster,size=RA))+
    geom_point(alpha=0.35) +
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('AitDist Sample compositions - tsne')
  
  set.seed(1234)
  tsneCoords = Rtsne(X = as.dist(aitDistHere),is_distance = TRUE, dims = 3,
                     perplexity=100, max_iter = 2000,learning=100)
  dfBiotic_Wider_MDS$tsne1 = tsneCoords$Y[,1]
  dfBiotic_Wider_MDS$tsne2 = tsneCoords$Y[,2]
  dfBiotic_Wider_MDS$tsne3 = tsneCoords$Y[,3]
  dfBiotic_Wider_MDS = dfBiotic_Wider_MDS %>% mutate(Cluster_num=as.numeric(Cluster))
  
  axx <- list(
    title = "tsne1"
  )
  
  axy <- list(
    title = "tsne2"
  )
  
  axz <- list(
    title = "tsne3"
  )
  
  out$TSNE_Sample_Ait_3d <- plotly::plot_ly(
    x=dfBiotic_Wider_MDS$tsne1,
    y=dfBiotic_Wider_MDS$tsne2,
    z=dfBiotic_Wider_MDS$tsne3,
    color = as.numeric(dfBiotic_Wider_MDS$Cluster_num),
    type="scatter3d",
    mode="markers") %>%
    #plotly::layout(
    #  scene = list(aspectratio = list(x = 4, y = 4, z = 3))) %>% 
    #list(xaxis=list(title='MDS1'),yaxis=list(title='MDS2'),zaxis=list(title='MDS3'))) %>% 
    plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
    plotly::hide_colorbar() 
  
  ###################-------------------------------------------------------------------------
  ##### Umap --------------------------
  ###################-------------------------------------------------------------------------

  clr_transf <- robCompositions::cenLR(x = dfBiotic_Wider_MDS %>% select(starts_with('ASVG_')) %>% data.frame())
  #test_duplicated = duplicated(clr_transf) there are no duplicates!
  #dim(clr_transf$x.clr)
  clr_transf = bind_cols(
    ID_Sample = dfBiotic_Wider_MDS$ID_Sample,
    clr_transf$x.clr)
  
  custom.config <- umap.defaults
  custom.config$n_neighbors=25
  custom.config$min_dist=0.05
  custom.config$spread = 0.5
  custom.config$random_state=10
  
  set.seed(1234)
  umap_obj <- umap(d = clr_transf %>% select(-ID_Sample),config = custom.config)
  #umap_obj <- umap(d = as.matrix(aitDistHere),input="dist")
  #plot(umap_obj$layout)
  dfBiotic_Wider_MDS = dfBiotic_Wider_MDS %>%
    mutate(umap1=umap_obj$layout[,1],umap2=umap_obj$layout[,2])
  
  
  out$umap_Sample_Ait_2d <- ggplot(
    data = dfBiotic_Wider_MDS,
    mapping = aes(x=umap1,y=umap2,size=RA,color=Cluster))+
    geom_point(alpha=0.35) +
    theme_minimal(base_size = 12)+
    theme(legend.position = 'bottom')+
    ggtitle('AitDist Sample compositions - umap - n_neighbors=15 - min_dist=0.25-spread = 0.5')
  
  
  custom.config$n_components=3
  set.seed(1234)
  umap_obj <- umap(d = clr_transf %>% select(-ID_Sample),config = custom.config)
  #umap_obj <- umap(d = as.matrix(aitDistHere),input="dist")
  #plot(umap_obj$layout)
  dfBiotic_Wider_MDS = dfBiotic_Wider_MDS %>%
    mutate(umap1=umap_obj$layout[,1],umap2=umap_obj$layout[,2],umap3=umap_obj$layout[,3])
  
  
  axx <- list(
    title = "umap1"
  )
  
  axy <- list(
    title = "umap2"
  )
  
  axz <- list(
    title = "umap3"
  )
  
  out$umap_Sample_Ait_3d <- plotly::plot_ly(
    x=dfBiotic_Wider_MDS$umap1,
    y=dfBiotic_Wider_MDS$umap2,
    z=dfBiotic_Wider_MDS$umap3,
    color = as.numeric(dfBiotic_Wider_MDS$Cluster_num),
    type="scatter3d",
    mode="markers") %>%
    #plotly::layout(
    #  scene = list(aspectratio = list(x = 4, y = 4, z = 3))) %>% 
    #list(xaxis=list(title='MDS1'),yaxis=list(title='MDS2'),zaxis=list(title='MDS3'))) %>% 
    plotly::layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>% 
    plotly::hide_colorbar() 
  
  ############################################################_----------------------------------------------
  ############################################################-----------------------------------------------
  ##### Ordination with aitichison
  dfBiotic_Wider = dfBiotic %>% 
    select(ID_Sample,Cluster,RA) %>%
    pivot_wider(names_from = Cluster,values_from = RA)
  
  dfBiotic_Wider_MDS <- dfBiotic_Wider %>% 
    mutate(across(where(is.numeric),function(x) x+0.00000001)) %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) 
  
   mds_obj <- cmdscale(d = robCompositions::aDist(dfBiotic_Wider_MDS[,-1]),
                       k = 3,eig = T,add=T)
   pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
   
   dfBiotic_Wider_MDS <- dfBiotic_Wider_MDS %>% mutate(
     MDS1 = mds_obj$points[,1],
     MDS2 = mds_obj$points[,2],
     MDS3 = mds_obj$points[,3]
   )
   
   dfBiotic_Wider_MDS = dfBiotic_Wider_MDS %>% 
     left_join(dfLonger %>% select(ID_Sample,Longhurst_Long) %>% distinct())
   
   out$MDS2d_Ait = ggplot(
     data = dfBiotic_Wider_MDS,
     mapping = aes(x=MDS1,y=MDS2,
                   color=Longhurst_Long))+
     geom_point(alpha=0.7,size=3) +
     theme_minimal(base_size = 12) +
     xlab(paste('MDS1 -',pct_explained[1],'%'))+
     ylab(paste('MDS2 -',pct_explained[2],'%'))+
     scale_fill_manual(values = longhurst_colours)+
     theme(legend.position = 'bottom')
   
   
   out$MDS2d_Ait_MostAbundant = 
     dfBiotic_Wider_MDS %>% 
     left_join(dfAbioMostAbundant %>% 
                 select(ID_Sample,Cluster,RA)) %>% 
     ggplot(
       mapping = aes(x=MDS1,y=MDS2,
                     color=Cluster,size=RA))+
     geom_point(alpha=0.7) +
     theme_minimal() +
     xlab(paste('MDS1 -',pct_explained[1],'%'))+
     ylab(paste('MDS2 -',pct_explained[2],'%'))+
     scale_fill_manual(values = colorASVgs)+
     theme(legend.position = 'bottom')
   
   ##### Ordination with bray
   dfBiotic_Wider = dfBiotic %>% 
     select(ID_Sample,Cluster,RA) %>%
     pivot_wider(names_from = Cluster,values_from = RA)
   
   mds_obj <- cmdscale(d = vegan::vegdist(dfBiotic_Wider[,-1],method = 'bray'),
                       k = 3,eig = T,add=T)
   pct_explained = round(100 * mds_obj$eig/sum(mds_obj$eig),1)
   
   dfBiotic_Wider <- dfBiotic_Wider %>% mutate(
     MDS1 = mds_obj$points[,1],
     MDS2 = mds_obj$points[,2],
     MDS3 = mds_obj$points[,3]
   )
   
   dfBiotic_Wider = dfBiotic_Wider %>% 
     left_join(dfLonger %>% select(ID_Sample,Longhurst_Long) %>% distinct())
   
   
   out$MDS2d_Bray = ggplot(
     data = dfBiotic_Wider,
     mapping = aes(x=MDS1,y=MDS2,
                   color=Longhurst_Long,
                   label=Longhurst_Long))+
     geom_point(alpha=0.7,size=4) +
     theme_minimal() +
     xlab(paste('MDS1 -',pct_explained[1],'%'))+
     ylab(paste('MDS2 -',pct_explained[2],'%'))+
     scale_fill_manual(values = longhurst_colours)+
     theme(legend.position = 'bottom')
   
   out$MDS2d_Bray_MostAbundantdf =
     dfBiotic_Wider %>% 
     left_join(dfAbioMostAbundant %>% 
                 select(ID_Sample,Cluster,RA)) %>% 
     ggplot(
       mapping = aes(x=MDS1,y=MDS2,
                     color=Cluster,size=RA))+
     geom_point(alpha=0.7) +
     theme_minimal() +
     xlab(paste('MDS1 -',pct_explained[1],'%'))+
     ylab(paste('MDS2 -',pct_explained[2],'%'))+
     scale_color_manual(values = colorASVgs)+
     theme(legend.position = 'bottom')
  
  #################################################
  ### Plots with the stacked version.
  #################################################
  
  dfBiotic_longer = dfBiotic %>% 
    left_join(dfLonger %>% select(ID_Sample,Latitude,Longitude,Depth,
                                  Cruise,Salinity,Oxygen,Temperature,
                                  Silicate,NO2,NO3,NOx,NH3) %>% distinct()) %>% 
    mutate(LogDepth=log(Depth+1))
  
  out$LatDepth_stacked =  dfBiotic_longer %>% filter(RA>0) %>% 
    select(Latitude,LogDepth,Longitude,Cluster,RA,all_of(abioticVars)) %>% 
    ggplot(aes(y=LogDepth,Latitude,size=RA,color=Cluster))+
    geom_point(alpha=0.25) + 
    scale_y_reverse()+
    theme_minimal()+
    facet_wrap(~Cluster,ncol=4)+
    theme(legend.position = 'bottom')
  
  out$LatLong_stacked = dfBiotic_longer %>% filter(RA>0) %>% 
    select(Latitude,LogDepth,Longitude,Cluster,RA,all_of(abioticVars)) %>% 
    ggplot(aes(y=Longitude,Latitude,size=RA,color=Cluster))+
    geom_point(alpha=0.25) + 
    theme_minimal()+
    facet_wrap(~Cluster,ncol=4)+
    theme(legend.position = 'bottom')
  
  out$DepthLong_stacked =dfBiotic_longer %>% filter(RA>0) %>% 
    select(Latitude,LogDepth,Longitude,Cluster,RA,all_of(abioticVars)) %>% 
    ggplot(aes(y=LogDepth,Longitude,size=RA,color=Cluster))+
    geom_point(alpha=0.25) + 
    scale_y_reverse()+
    theme_minimal()+
    facet_wrap(~Cluster,ncol=4)+
    theme(legend.position = 'bottom')
  
  ######################################################
  # Need To add the new ~conditinal~ probability here yet.
  ######################################################
  
  return(out)
}
