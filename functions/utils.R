################# 
# small funcitons
#################

library(dplyr)

# Not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# Dframe<- data.table::fread('/Users/rafaelcatoia/Desktop/repos/phaeocystis/grump.phaeocystis_asv_long.csv') %>% 
#  filter(Cruise!="MOSAiC" )

# Create composition and filter two ASV with no abundance
tidy_grump <- function(Dframe,removeMOSAiC=T,vet_ASVs2remove=NULL){
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
  
  ASVs_with_Species = Dframe %>% filter(!is.na(Species)) %>% select(ASV_name,Species) %>% distinct()
  
  #### to be able to make the
  value_toAdd = min(Dframe$Raw.Sequence.Counts)/1000
  
  output = list(
    
    ## Id with all ASVs
    id_ASVs = idAsvs,
    
    ## Id with all SamplesId and sample names
    id_Samples = Dframe %>%  select(SampleKey,SampleID) %>% distinct(),
    
    ## ASV composition (composition on samples, transposing and than composition on asvs)
    ASV_composition = Dframe %>% #filter(is.na(Species)) %>% 
      select(SampleKey,ASV_name,Raw.Sequence.Counts) %>% 
      pivot_wider(id_cols = c('SampleKey'),
                  values_from ='Raw.Sequence.Counts',
                  names_from='ASV_name',values_fill = value_toAdd) %>%
      arrange(SampleKey) %>% 
      mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% ### Here we have the compositions on the samples
      #mutate(across(where(is.numeric),function(x) {x+0.00000000001}))%>%  ### Adding a small value in all, in order to have the CLR transformation / aitichison distance
      #mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% ### composition over the samples gain
      pivot_longer(cols = any_of(idAsvs)) %>% ## Stacking to transpose 
      pivot_wider(id_cols='name',
                  values_from = value,
                  names_from='SampleKey') %>% ## Wider again, ASVs in the row
      mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% # each row is an ASVs and the columns are the sample. The rows now add to 1.
      arrange(name) %>% data.frame(),
    
    
    code2comp = "select(SampleKey,ASV_name,Raw.Sequence.Counts) %>% 
      pivot_wider(id_cols = c('SampleKey'),
                  values_from ='Raw.Sequence.Counts',
                  names_from='ASV_name',values_fill = 0) %>%
      arrange(SampleKey) %>% 
      mutate(across(where(is.numeric))/
               rowSums(across(where(is.numeric)))) %>% ### Here we have the compositions on the samples
      mutate(across(where(is.numeric),function(x) {x+0.00000000001}))%>%  ### Adding a small value in all, in order to have the CLR transformation / aitichison distance
      mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% ### composition over the samples gain
      pivot_longer(cols = any_of(idAsvs)) %>% ## Stacking to transpose 
      pivot_wider(id_cols='name',
                  values_from = value,
                  names_from='SampleKey') %>% ## Wider again, ASVs in the row
      mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% # each row is an ASVs and the columns are the sample. The rows now add to 1.
      arrange(name) %>% data.frame()",
    
    dframe = Dframe
  )
  
  return(output)
  
}

