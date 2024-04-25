### packages ---
library(dplyr) ; library(tidyr) ; library(ggplot2) 

### loading functions
root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")
savingdir_oceanSlices = root$find_file("saved_files_oceanSlices")
datapath = root$find_file(paste0(datadir,'/','grump.phaeocystis_asv_long.csv'))
files_vec <- list.files(funsdir)

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

### loading data
dframe = data.table::fread(input = datapath) %>% filter(Cruise!="MOSAiC",Raw.Sequence.Counts>0)
# dframe_allASVs = tidy_grump(Dframe = dframe)

### ------------------------------------
## asvs to remove
asvs_to_remove_oneSample = readRDS(file = paste0(savingdir_oceanSlices,'/','asvs_to_remove_oneSample'))
asvs_to_remove_twoSample = readRDS(file = paste0(savingdir_oceanSlices,'/','asvs_to_remove_twoSample'))

## If we want to remove ASVs that are rare
## include filter(ASV_name %!in% asvs_to_remove_XXXSample) where XXX is one or two 

###### ------------------------------------------------------------- data handling
# First quantiles of depth
distinct_depth <- dframe %>% select(SampleID,Depth) %>% distinct()
qtls_0.1 = quantile(distinct_depth$Depth, seq(0,1,0.1))
qtlsNames = 1:(length(qtls_0.1)-1)
qtlsNames = ifelse(qtlsNames<10,paste0('DGR0',qtlsNames),paste0('DGR',qtlsNames))

LatSlices = seq(-75,85,10)
qtlsNames_lat = 1:(length(LatSlices)-1)
qtlsNames_lat = ifelse(qtlsNames_lat<10,paste0('LatSlice0',qtlsNames_lat),paste0('LatSlice',qtlsNames_lat))


### here we have the stacked grump by unitMeasurement (ocean sites/ocean slices)
### unitMeasurement, ASV and Raw Sequence Counts

## This is using LH and Depth ------------------------------------------------------------------------------------

stacked_LH_CatDepth = dframe  %>%  ## insert the filter here ------
  mutate(CatDepth=cut(Depth,breaks = qtls_0.1,include.lowest = T,right = T,labels = qtlsNames)) %>% 
  mutate(UnitMeasurement = paste0(Longhurst_Short,'_',CatDepth)) %>% 
  select(UnitMeasurement,ASV_name,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(UnitMeasurement,ASV_name) %>% 
  summarise(Raw.Sequence.Counts = sum(Raw.Sequence.Counts)) %>% ungroup() %>% ### Now getting the composition/RA by adding all RawSequence
  group_by(UnitMeasurement) %>% 
  mutate(AllRawCounts = sum(Raw.Sequence.Counts)) %>% ungroup() %>% 
  select(UnitMeasurement,ASV_name,Raw.Sequence.Counts) %>%  ungroup()

## finding the value to add to be able to find CLR and Ait Dist
min(stacked_LH_CatDepth$Raw.Sequence.Counts)
valueToAdd =min(stacked_LH_CatDepth$Raw.Sequence.Counts)*1e-06

wide_LH_CatDepth = stacked_LH_CatDepth %>% 
  select(UnitMeasurement,Raw.Sequence.Counts,ASV_name) %>% 
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts+valueToAdd) %>% 
  pivot_wider(id_cols = 'UnitMeasurement',
              values_from ='Raw.Sequence.Counts',
              names_from='ASV_name',values_fill = valueToAdd) %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% #Compositions of unitMeasurements were created
  pivot_longer(cols = -1) %>% #Stacking to transpose 
  pivot_wider(id_cols='name',values_from = value,
              names_from='UnitMeasurement') %>% ## Wider again, ASVs in the row
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% # each row is an ASVs and the columns are combinations LH and depth slices. The rows now 
  arrange(name) %>% data.frame()

## This is using and Depth ------------------------------------------------------------------------------------

dfLonger_LatDepth_Group = dframe  %>% ### add filter of ASVs here
  mutate(CatDepth=cut(Depth,breaks = qtls_0.1,include.lowest = T,right = T,labels = qtlsNames)) %>% 
  mutate(LatSlices = cut(Latitude,breaks = LatSlices,include.lowest = T,right = T,labels = qtlsNames_lat)) %>% 
  mutate(UnitMeasurement = paste0(LatSlices,'_',CatDepth)) %>% 
  select(UnitMeasurement,ASV_name,Raw.Sequence.Counts) %>% distinct() %>% 
  group_by(UnitMeasurement,ASV_name) %>% 
  summarise(Raw.Sequence.Counts = sum(Raw.Sequence.Counts)) %>% ### Now getting the composition/RA by adding all RawSequence
  ungroup() %>% 
  group_by(UnitMeasurement) %>% 
  mutate(AllRawCounts = sum(Raw.Sequence.Counts)) %>% ungroup()


valueToAdd =min(dfLonger_LatDepth_Group$Raw.Sequence.Counts)*1e-06


dfWiderComposition = dfLonger_LatDepth_Group %>% 
  select(UnitMeasurement,Raw.Sequence.Counts,ASV_name) %>% 
  mutate(Raw.Sequence.Counts=Raw.Sequence.Counts+valueToAdd) %>% 
  pivot_wider(id_cols = 'UnitMeasurement',
              values_from ='Raw.Sequence.Counts',
              names_from='ASV_name',values_fill = valueToAdd) %>% 
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% ## now transposing
  pivot_longer(cols = -1) %>% #Stacking to transpose 
  pivot_wider(id_cols='name',values_from = value,
              names_from='UnitMeasurement') %>% ## Wider again, ASVs in the row
  mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% # each row is an ASVs and the columns are the sample. The rows now 
  arrange(name) %>% data.frame()

