################# 
# small funcitons
#################

library(dplyr)

# Not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# Dframe<- data.table::fread('/Users/rafaelcatoia/Desktop/repos/phaeocystis/grump.phaeocystis_asv_long.csv') %>% 
#  filter(Cruise!="MOSAiC" )

# Create composition and filter two ASV with no abundance
tidy_grump <- function(Dframe,removeMOSAiC=T,vet_ASVs2remove=NULL,binnary=F){
  #ASVs_toRemove = c('GRUMP203261','GRUMP203863') #ASVs that had raw sequence counts equal to zero
  #Dframe = Dframe %>% filter(ASV_name%!in%ASVs_toRemove) # filtering those asvs
  Dframe = Dframe %>% filter(Raw.Sequence.Counts>0)
  if(removeMOSAiC){
    Dframe = Dframe %>% filter(Cruise!="MOSAiC" )
  }
  
  if(!is.null(vet_ASVs2remove)){
    vet_ASVs2remove = unlist(vet_ASVs2remove)
    Dframe = Dframe %>% filter(ASV_name%!in%vet_ASVs2remove)
  }
  
  idAsvs = unique(Dframe$ASV_name) # listing all ASVs
  idSamples = unique(Dframe$SampleID) # listing all the sample names
  
  ## just creating new sample names to be easiear to plot.
  vetNumbers <- 1:length(idSamples)
  SampleKeys = ifelse(vetNumbers<10,paste('000',vetNumbers,sep=''),
                      ifelse(vetNumbers<100,paste('00',vetNumbers,sep=''),
                             ifelse(vetNumbers<1000,paste('0',vetNumbers,sep=''),paste(vetNumbers))))
  
  SampleKeysNames = unique(paste(factor(Dframe$SampleID,
                                        labels = paste('S',SampleKeys,sep=''))))
  SampleKeysNames = SampleKeysNames[order(SampleKeysNames)]
  Dframe$SampleKey = paste(factor(Dframe$SampleID,labels = paste('S',SampleKeys,sep='')))
  
  Dframe = Dframe %>%  mutate(Species=ifelse(Species=='Phaeocystis_sp.',NA,Species))
  
  ASVs_Species = Dframe %>% select(ASV_name,Species) %>% distinct() %>% arrange(Species)
  
  #### to be able to make the
  value_toAdd = min(Dframe$Raw.Sequence.Counts)/1000
  
  output =  list()
    
  ## Id with all ASVs
  id_ASVs = idAsvs
  
  ## Id with all SamplesId and sample names
  id_Samples = Dframe %>%  select(SampleKey,SampleID) %>% distinct()
  
  ## ASV composition (composition on samples, transposing and than composition on asvs)
  Sample_Composition = Dframe %>% #filter(is.na(Species)) %>% 
    select(SampleKey,ASV_name,Raw.Sequence.Counts) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts) %>% 
    pivot_wider(id_cols = c('SampleKey'),
                values_from ='Raw.Sequence.Counts',
                names_from='ASV_name',values_fill = 0) %>%
    arrange(SampleKey) %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    tibble %>% select("SampleKey", sort(colnames(.)))
  
  
  ## ASV composition (composition on samples, transposing and than composition on asvs)
  Sample_Composition_Filled = Dframe %>% #filter(is.na(Species)) %>% 
    select(SampleKey,ASV_name,Raw.Sequence.Counts) %>% 
    mutate(Raw.Sequence.Counts=Raw.Sequence.Counts+value_toAdd) %>% 
    pivot_wider(id_cols = c('SampleKey'),
                values_from ='Raw.Sequence.Counts',
                names_from='ASV_name',values_fill = value_toAdd) %>%
    arrange(SampleKey) %>% 
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
    tibble %>% select("SampleKey", sort(colnames(.)))
  
  ASV_composition = Sample_Composition_Filled %>% pivot_longer(cols = any_of(idAsvs)) %>% ## Stacking to transpose 
    pivot_wider(id_cols='name',
                values_from = value,
                names_from='SampleKey') %>% ## Wider again, ASVs in the row
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% # each row is an ASVs and the columns are the sample. The rows now add to 1.
    arrange(name)
  
  
  dframe = Dframe
  
  if(binnary){
  binnary_ASV_df = Dframe %>% #filter(is.na(Species)) %>% 
    select(ASV_name,SampleKey) %>% distinct() %>% mutate(binValue = 1) %>%
    arrange(ASV_name) %>%
    pivot_wider(id_cols = 'ASV_name',
                values_from ='binValue',
                names_from='SampleKey',values_fill = 0) %>% 
    arrange(ASV_name) %>% 
    tibble %>% select("ASV_name", sort(colnames(.)))
  
  output$binnary_ASV_df = binnary_ASV_df
  }
    
  rawCounts_ASV_df = Dframe %>%
    select(ASV_name,SampleKey,Raw.Sequence.Counts) %>% distinct() %>% 
    arrange(ASV_name) %>%
    pivot_wider(id_cols = 'ASV_name',
                values_from ='Raw.Sequence.Counts',
                names_from='SampleKey',values_fill = 0) %>% 
    arrange(ASV_name) %>% 
    tibble %>% select("ASV_name", sort(colnames(.)))
  
  
  ## arranging output 
  output$id_ASVs_Species = ASVs_Species
  output$id_Samples = id_Samples
  output$dframe = Dframe
  output$Sample_Composition = Sample_Composition
  output$ASV_composition = ASV_composition
  output$rawCounts_ASV_df = rawCounts_ASV_df
  output$Sample_Composition_Filled = Sample_Composition_Filled
  return(output)
  
}

###### Hamming Dist using dot product
hammingDist <- function(X) {
  D <- (1 - X) %*% t(X)
  D + t(D)
}

###### normalizing dist funciton
normalizeDistMatrix <- function(distMatrix,typeNorm='2'){
  distMatrix <- as.matrix(distMatrix)
  normD = norm(distMatrix,type = typeNorm)
  return(distMatrix/normD)
}


## Justins hamming dist -- 
hammdist <- function(p1, p2){
  sum(p1!=p2)
}
hamm_distmat <- function(dfmat){
  distmat = matrix(NA, ncol = nrow(dfmat), nrow = nrow(dfmat))
  for(irow1 in 1:nrow(dfmat)){
    for(irow2 in 1:nrow(dfmat)){
      if(irow1 < irow2){
        distmat[irow1, irow2] = hammdist(dfmat[irow1,], dfmat[irow2,])
      }
    }
  }
  return(distmat)
}
