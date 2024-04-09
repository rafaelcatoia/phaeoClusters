# create the induced splitting matrix
library('dplyr')

inducedDist <- function(dFrame,c1 = 1, c2 = 2, c3=3, c4=4, normtype = '2',compMatrix){
  
  #c1 is within species - but no unnassinged
  #c2 is between species - but no unnassinged
  #c3 is between species  and unnassinged
  #c4 is between ASV 
  
  output <- list()
  
  ## Function we'll need here
  distMatrixToDF <- function(dist_matrix,vetNames=''){
    # Ensure input is a proper matrix
    if (!is.matrix(dist_matrix)) {
      stop("Input must be a matrix.")
    }
    
    # Ensure it's square and symmetric (characteristics of a distance matrix)
    if (nrow(dist_matrix) != ncol(dist_matrix)) {
      stop("Input matrix must be square.")
    }
    if (!all(dist_matrix == t(dist_matrix))) {
      stop("Input matrix must be symmetric.")
    }
    
    # Get lower triangular indices (excluding the diagonal)
    ind <- which(lower.tri(dist_matrix), arr.ind = TRUE)
    
    # Create the data frame
    df <- data.frame(
      item1 = ind[, 1],
      item2 = ind[, 2],
      distance = as.vector(dist_matrix[ind])
    )
    
    if( any(vetNames!='') ){
      df = df %>% transmute(
        ASV1 = as.vector(vetNames[ind[, 2]]),
        ASV2 = as.vector(vetNames[ind[, 1]]),
        distance = distance
      )
    }
    return(df)
  }
  
  ###################################--
  ASVorder <- unique(compMatrix[,1])
  ASVorder <- ASVorder[order(ASVorder)]
  ##stacked version of distance matrix
  AitMatrix = robCompositions::aDist(compMatrix[,-1])
  aitStacked <- distMatrixToDF(dist_matrix = as.matrix(AitMatrix),vetNames = compMatrix[,1])
  
  myDist = expand.grid(ASV1=unique(dFrame$ASV_name),ASV2=unique(dFrame$ASV_name)) %>% 
    left_join(dFrame %>% select(ASV_name,Species) %>% rename('ASV1'=1,'Species_ASV1'=2) %>% distinct()) %>% 
    left_join(dFrame %>% select(ASV_name,Species) %>% rename('ASV2'=1,'Species_ASV2'=2) %>% distinct()) %>% 
    mutate(Species_ASV1=ifelse(is.na(Species_ASV1),'noSpecies',Species_ASV1),
           Species_ASV2=ifelse(is.na(Species_ASV2),'noSpecies',Species_ASV2)) %>% 
    
    mutate(
      mdist = ifelse(
        Species_ASV1==Species_ASV2,c1,c2)
    ) %>% 
    mutate(
      mdist=ifelse(Species_ASV1 == 'noSpecies' | Species_ASV2 == 'noSpecies',c3,mdist)
    ) %>% 
    mutate(
      mdist=ifelse(Species_ASV1 == 'noSpecies' & Species_ASV2 == 'noSpecies',c4,mdist)
    ) %>% 
    mutate(
      mdist=ifelse(ASV1 == ASV2,0,mdist)
    )
  
  
  output$myDist_checking = myDist %>% select(Species_ASV1,Species_ASV2,mdist) %>% 
    arrange(Species_ASV1,Species_ASV2) %>% 
    distinct() %>% knitr::kable() %>% kableExtra::kable_styling()
  
  ## inserting AitDist
  myDist = myDist %>%
    left_join(aitStacked) %>% 
    left_join(aitStacked %>% rename('ASV2'=1,'ASV1'=2,'dist2'=3)) %>% 
    mutate(distAit = pmax(distance,dist2,na.rm = T)) %>% 
    select(ASV1,ASV2,Species_ASV1,Species_ASV2,mdist,distAit) %>% 
    mutate(distAit=ifelse(is.na(distAit),0,distAit))
  ##
  
  ## Transforming in matrix
  myDist_matrix = myDist %>% arrange(ASV1,ASV2) %>% 
    pivot_wider(
      names_from  = ASV1,
      id_cols=ASV2,
      values_from = mdist,
      values_fill = 0
    ) %>% arrange(ASV2) %>% select(ASV2,any_of(ASVorder))
  
  normMyDist <- norm(myDist_matrix[,-1],type = normtype)
  
  aitDist_matrix = myDist %>% arrange(ASV1,ASV2) %>% 
    pivot_wider(
      names_from  = ASV1,
      id_cols=ASV2,
      values_from = distAit
    ) %>% arrange(ASV2) %>% select(ASV2,any_of(ASVorder))
  normAitDist <- norm(aitDist_matrix[,-1],type = normtype)
  
  output$normalizedMyDist = myDist_matrix[-1]/normMyDist
  output$normalizedAitDist = aitDist_matrix[-1]/normAitDist
  return(output)
}