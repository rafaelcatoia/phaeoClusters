##### ------------------------------------------------------------
## Clustering ASV over Ocean Slices
##### ------------------------------------------------------------

# packages
library(dplyr) ; library(tidyr) ; library(ggplot2) 
#library(doParallel) ; library(foreach) ; library(doSNOW)

# files
root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")
savingdirOceanSlices = root$find_file("saved_files_oceanSlices")

datapath = root$find_file(paste0(datadir,'/','grump.phaeocystis_asv_long.csv'))
files_vec <- list.files(funsdir)
currentwd <- getwd()


for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

# data 
dframe = data.table::fread(input = datapath) %>% filter(Cruise!="MOSAiC",Raw.Sequence.Counts>0)
dframe_allASVs = tidy_grump(Dframe = dframe)

## Ocean slices composition
Lat_Depth_Composition_df=readRDS(file = paste0(savingdir,'/','Lat_Depth_Composition_df'))
LH_Depth_Composition_df=readRDS(file = paste0(savingdir,'/','LH_Depth_Composition_df'))

##----------------------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------------------
## Creating clusters and evaluating metrics::---------------------------------------------------------------------------------------------------------
runningLatDepth = F
###### Lat_Depth -------------------------------------------------------------------------------------------------------------------------------------

inducedDist_ = inducedDist(
  dFrame = dframe_allASVs$dframe,
  c1 = 1,c2=1000,c3=10,c4=10,
  compMatrix = dframe_allASVs$ASV_composition)

if(runningLatDepth){
  ########### For Alpha = 0 

  obj_eval <- create_eval(dfComposition = Lat_Depth_Composition_df)
  saveRDS(obj_eval$list_results_hclust,file = paste0(savingdirOceanSlices,'/','list_LatDepth_eval_alpha0_hclust'))
  saveRDS(obj_eval$list_results_medoid,file = paste0(savingdirOceanSlices,'/','list_LatDepth_eval_alpha0_pam'))
  
  ########### For Alpha = 0.1
  obj_eval <- create_eval(dfComposition = Lat_Depth_Composition_df,inducedMatrixDist = inducedDist_$normalizedMyDist,alpha_ = 0.1)
  saveRDS(obj_eval$list_results_hclust,file = paste0(savingdirOceanSlices,'/','list_LatDepth_eval_alpha0.1_hclust'))
  saveRDS(obj_eval$list_results_medoid,file = paste0(savingdirOceanSlices,'/','list_LatDepth_eval_alpha0.1_pam'))
  
  ########### For Alpha = 0.25
  obj_eval <- create_eval(dfComposition = Lat_Depth_Composition_df,inducedMatrixDist = inducedDist_$normalizedMyDist,alpha_ = 0.25)
  saveRDS(obj_eval$list_results_hclust,file = paste0(savingdirOceanSlices,'/','list_LatDepth_eval_alpha0.25_hclust'))
  saveRDS(obj_eval$list_results_medoid,file = paste0(savingdirOceanSlices,'/','list_LatDepth_eval_alpha0.25_pam'))
}

##----------------------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------------------
## Creating clusters and evaluating metrics::---------------------------------------------------------------------------------------------------------
runningLHDepth = F
###### LH_Depth --------------------------------------------------------------------------------------------------------------------------------------
if(runningLHDepth){
  
  ########### For Alpha = 0
  obj_eval <- create_eval(dfComposition = LH_Depth_Composition_df)
  saveRDS(obj_eval$list_results_hclust,file = paste0(savingdirOceanSlices,'/','list_LH_Depth_eval_alpha0_hclust'))
  saveRDS(obj_eval$list_results_medoid,file = paste0(savingdirOceanSlices,'/','list_LH_Depth_eval_alpha0_pam'))
  
  ########### For Alpha = 0.1
  obj_eval <- create_eval(dfComposition = LH_Depth_Composition_df,inducedMatrixDist = inducedDist_$normalizedMyDist,alpha_ = 0.1)
  saveRDS(obj_eval$list_results_hclust,file = paste0(savingdirOceanSlices,'/','list_LH_Depth_eval_alpha0.1_hclust'))
  saveRDS(obj_eval$list_results_medoid,file = paste0(savingdirOceanSlices,'/','list_LH_Depth_eval_alpha0.1_pam'))
  
  ########### For Alpha = 0.25
  obj_eval <- create_eval(dfComposition = LH_Depth_Composition_df,inducedMatrixDist = inducedDist_$normalizedMyDist,alpha_ = 0.25)
  saveRDS(obj_eval$list_results_hclust,file = paste0(savingdirOceanSlices,'/','list_LH_Depth_eval_alpha0.25_hclust'))
  saveRDS(obj_eval$list_results_medoid,file = paste0(savingdirOceanSlices,'/','list_LH_Depth_eval_alpha0.25_pam'))
}



#### --------------------------------------------------------------------------------------------------------------------------------------------- ####
#### --------------------------------------------------------------------------------------------------------------------------------------------- ####
#### --------------------------------------------------------------------------------------------------------------------------------------------- ####
#### --------------------------------------------------------------------------------------------------------------------------------------------- ####
#### Plots --------------------------------------------------------------------------------------------------------------------------------------- ####

list_LatDepth_eval_alpha0_hclust    <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LatDepth_eval_alpha0_hclust"   ))
list_LatDepth_eval_alpha0.1_hclust  <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LatDepth_eval_alpha0.1_hclust" ))
list_LatDepth_eval_alpha0.25_hclust <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LatDepth_eval_alpha0.25_hclust"))
list_LatDepth_eval_alpha0_pam       <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LatDepth_eval_alpha0_pam"      ))
list_LatDepth_eval_alpha0.1_pam     <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LatDepth_eval_alpha0.1_pam"    ))
list_LatDepth_eval_alpha0.25_pam    <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LatDepth_eval_alpha0.25_pam"   ))
list_LH_Depth_eval_alpha0_hclust    <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LH_Depth_eval_alpha0_hclust"   ))
list_LH_Depth_eval_alpha0.1_hclus   <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LH_Depth_eval_alpha0.1_hclust" ))
list_LH_Depth_eval_alpha0.25_hclust <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LH_Depth_eval_alpha0.25_hclust"))
list_LH_Depth_eval_alpha0_pam       <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LH_Depth_eval_alpha0_pam"      ))
list_LH_Depth_eval_alpha0.1_pam     <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LH_Depth_eval_alpha0.1_pam"    ))
list_LH_Depth_eval_alpha0.25_pam    <- readRDS(file = paste0(savingdirOceanSlices,'/',"list_LH_Depth_eval_alpha0.25_pam"))


eval_clusters_full = bind_rows(
  summarise_adist_clustering(list_LatDepth_eval_alpha0_hclust   ) %>% mutate(sliceType='LatDepth',method='hclust',alpha=0),
  summarise_adist_clustering(list_LatDepth_eval_alpha0.1_hclust ) %>% mutate(sliceType='LatDepth',method='hclust',alpha=0.10),
  summarise_adist_clustering(list_LatDepth_eval_alpha0.25_hclust) %>% mutate(sliceType='LatDepth',method='hclust',alpha=0.25),
  summarise_adist_clustering(list_LatDepth_eval_alpha0_pam      ) %>% mutate(sliceType='LatDepth',method='pam',alpha=0),
  summarise_adist_clustering(list_LatDepth_eval_alpha0.1_pam    ) %>% mutate(sliceType='LatDepth',method='pam',alpha=0.10),
  summarise_adist_clustering(list_LatDepth_eval_alpha0.25_pam   ) %>% mutate(sliceType='LatDepth',method='pam',alpha=0.25),
  
  summarise_adist_clustering(list_LH_Depth_eval_alpha0_hclust   ) %>% mutate(sliceType='LH_Depth',method='hclust',alpha=0),
  summarise_adist_clustering(list_LH_Depth_eval_alpha0.1_hclus  ) %>% mutate(sliceType='LH_Depth',method='hclust',alpha=0.10),
  summarise_adist_clustering(list_LH_Depth_eval_alpha0.25_hclust) %>% mutate(sliceType='LH_Depth',method='hclust',alpha=0.25),
  summarise_adist_clustering(list_LH_Depth_eval_alpha0_pam      ) %>% mutate(sliceType='LH_Depth',method='pam',alpha=0),
  summarise_adist_clustering(list_LH_Depth_eval_alpha0.1_pam    ) %>% mutate(sliceType='LH_Depth',method='pam',alpha=0.10),
  summarise_adist_clustering(list_LH_Depth_eval_alpha0.25_pam   ) %>% mutate(sliceType='LH_Depth',method='pam',alpha=0.25)
)



####


eval_clusters_full %>%
  select(nclust,alpha,method,sliceType,ratio_avg_within_between,avg_between_medoids,mean_avg_within) %>% 
  filter(nclust>2) %>% mutate(alpha=factor(alpha,ordered = T,labels = c(0,0.10,0.25))) %>% 
  pivot_longer(cols = -c(1,2,3,4),names_to = 'metric') %>%
  ggplot(aes(x = nclust, y = value,col=method,linetype=alpha)) +
  geom_line() +
  theme_minimal(base_size = 12) +
  facet_grid(metric~sliceType, scales = "free_y")+
  theme(legend.position = 'bottom')

eval_clusters_full %>%
  select(nclust,alpha,method,sliceType,max_max_dist_within,avg_max_dist_within,avg_min_dist_clust_i_other,ratio_avg_max_min) %>% 
  filter(nclust>3) %>% mutate(alpha=factor(alpha,ordered = T,labels = c(0,0.10,0.25))) %>% 
  pivot_longer(cols = -c(1,2,3,4),names_to = 'metric') %>%
  ggplot(aes(x = nclust, y = value,col=method,linetype=alpha)) +
  geom_line() +
  theme_minimal(base_size = 12) +
  facet_grid(metric~sliceType, scales = "free_y")+
  theme(legend.position = 'bottom')


#### --------------------------------------------------------------------------------------------------------------------------------------------- ####
#### --------------------------------------------------------------------------------------------------------------------------------------------- ####
#### --------------------------------------------------------------------------------------------------------------------------------------------- ####
#### --------------------------------------------------------------------------------------------------------------------------------------------- ####
#### Plots --------------------------------------------------------------------------------------------------------------------------------------- ####
aitDist = vegan::vegdist(x = LH_Depth_Composition_df %>% select(-name),method = 'aitchison')
hclust_obj <- hclust(d = as.dist(aitDist),method = 'ward.D')
plot(hclust_obj,labels = F)


aitDist = vegan::vegdist(x = Lat_Depth_Composition_df %>% select(-name),method = 'aitchison')
hclust_obj <- hclust(d = as.dist(aitDist),method = 'ward.D')
plot(hclust_obj,labels = F)


testing <- robCompositions::clustCoDa_qmode(x = Lat_Depth_Composition_df %>% select(-name))
# Warning message:
#   In robustbase::covMcd(z) : The covariance matrix has become singular during
# the iterations of the MCD algorithm.
# There are 1415 observations (in the entire dataset of 1482 obs.) lying on the hyperplane with equation a_1*(x_i1 - m_1) + ... + a_p*(x_ip - m_p) = 0
# with (m_1, ..., m_p) the mean of these observations and coefficients a_i from the vector 
# a <- c(0.4453373, 0.0172303, -0.606222, 0.1394027, -7.71e-05,
#        -7.78e-05, -0.1340071, 0.1715656, 0.3741709, 0.0032482, 0.0032747, 0.0033017, 0.0033291, 0.0033569, 0.0033853, 0.0034141, 0.0034434, 0.0034732,
#        0.0035035, 0.0035344, 0.0035658, 0.0035978, 0.0036303, 0.0036635, 0.0036973, 0.0037317, 0.0037667, 0.0038024, 0.0038388, 0.0038759, 0.0039137,
#        0.0039523, 0.0039916, 0.0040317, 0.0040726, 0.0041144, 0.004157, 0.0042006, 0.004245, 0.0042904, 0.0043368, 0.0043842, 0.0044327, 0.0044822,
#        0.0045328, 0.0045846, 0.0046376, 0.0046919, 0.0047474, 0.0048043, 0.0048625, 0.0049222, 0.0049833, 0.005046, 0.0051103, 0.0051762, 0.0052439,
#        0.0053133, 0.0053847, 0.0054579, 0.0055332,  [... truncated]



testing <- robCompositions::clustCoDa(x = Lat_Depth_Composition_df %>% select(-name),k=5)

max(Lat_Depth_Composition_df[,-1])
       