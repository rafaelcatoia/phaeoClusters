######## --- last script 

library(dplyr) ; library(tidyr) ; library(ggplot2) 
library(doParallel) ; library(foreach) ; library(doSNOW)

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

datapath = root$find_file(paste0(datadir,'/','grump.phaeocystis_asv_long.csv'))
files_vec <- list.files(funsdir)
currentwd <- getwd()

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

getwd()
funcitons_parallel = c("clusterize_ward_pam",
                       "eval_adist_clustering",
                       "inducedDist",
                       "summarise_adist_clustering",
                       "tidy_grump")

dframe = data.table::fread(input = datapath)
dframe = tidy_grump(Dframe = dframe)

### creating induced matrix
inducedDist_ = inducedDist(
  dFrame = dframe$dframe,
  c1 = 1,c2=1000,c3=10,c4=10,
  compMatrix = dframe$ASV_composition)

#### --------------------------
running=T

if(running){
  ## pure biotic factors = 
  list_with_clusters_alpha0 <- clusterize_ward_pam(
    distMatrix = inducedDist_$normalizedAitDist,
    maxnclust = 50,
    dFrameAsv = dframe$ASV_composition$name)
  
  alpha_=0.1
  
  distMatrix_ = alpha_*inducedDist_$normalizedMyDist+(1-alpha_)*inducedDist_$normalizedAitDist
  
  list_with_clusters_alpha0.10 <- clusterize_ward_pam(
    distMatrix = distMatrix_,
    maxnclust = 50,
    dFrameAsv = dframe$ASV_composition$name)
  
  alpha_=0.25
  
  distMatrix_ = alpha_*inducedDist_$normalizedMyDist+(1-alpha_)*inducedDist_$normalizedAitDist
  
  list_with_clusters_alpha0.25 <- clusterize_ward_pam(
    distMatrix = distMatrix_,
    maxnclust = 50,
    dFrameAsv = dframe$ASV_composition$name)
  
  alpha_=0.5
  
  distMatrix_ = alpha_*inducedDist_$normalizedMyDist+(1-alpha_)*inducedDist_$normalizedAitDist
  
  list_with_clusters_alpha0.50 <- clusterize_ward_pam(
    distMatrix = distMatrix_,
    maxnclust = 50,
    dFrameAsv = dframe$ASV_composition$name)
  
  saveRDS(object = list_with_clusters_alpha0,file = paste0(savingdir,'/','list_with_clusters_alpha0'))
  saveRDS(object = list_with_clusters_alpha0.10,file = paste0(savingdir,'/','list_with_clusters_alpha0.10'))
  saveRDS(object = list_with_clusters_alpha0.25,file = paste0(savingdir,'/','list_with_clusters_alpha0.25'))
  saveRDS(object = list_with_clusters_alpha0.50,file = paste0(savingdir,'/','list_with_clusters_alpha0.50'))
}


####### --- Evaluating metrics --------------------------------------- Hclust
running = T

if(running){
  nclusters <- ncol(list_with_clusters_alpha0$dFrameAsv_hclust)-1
  list_results <- list()
  
  for( i in 1:nclusters){
    list_results[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      list_with_clusters_alpha0$dFrameAsv_hclust[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dframe$ASV_composition))
    cat('iteration------:',i,'\n')
  }
  saveRDS(list_results,file = paste0(savingdir,'/','list_results_eval_alpha0_hclust'))
}

if(running){
  nclusters <- ncol(list_with_clusters_alpha0.10$dFrameAsv_hclust)-1
  list_results <- list()
  
  for( i in 1:nclusters){
    list_results[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      list_with_clusters_alpha0.10$dFrameAsv_hclust[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dframe$ASV_composition))
    cat('iteration------:',i,'\n')
  }
  saveRDS(list_results,file = paste0(savingdir,'/','list_results_eval_alpha0.10_hclust'))
}

if(running){
  nclusters <- ncol(list_with_clusters_alpha0.25$dFrameAsv_hclust)-1
  list_results <- list()
  
  for( i in 1:nclusters){
    list_results[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      list_with_clusters_alpha0.25$dFrameAsv_hclust[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dframe$ASV_composition))
    cat('iteration------:',i,'\n')
  }
  saveRDS(list_results,file = paste0(savingdir,'/','list_results_eval_alpha0.25_hclust'))
}

if(running){
  nclusters <- ncol(list_with_clusters_alpha0.50$dFrameAsv_hclust)-1
  list_results <- list()
  
  for( i in 1:nclusters){
    list_results[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      list_with_clusters_alpha0.50$dFrameAsv_hclust[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dframe$ASV_composition))
    cat('iteration------:',i,'\n')
  }
  saveRDS(list_results,file = paste0(savingdir,'/','list_results_eval_alpha0.50_hclust'))
}

####### --- Evaluating metrics --------------------------------------- pam

running =T

if(running){
  nclusters <- ncol(list_with_clusters_alpha0$dFrameAsv_pam)-1
  list_results <- list()
  
  for( i in 1:nclusters){
    list_results[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      list_with_clusters_alpha0$dFrameAsv_pam[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dframe$ASV_composition))
    cat('iteration------:',i,'\n')
  }
  saveRDS(list_results,file = paste0(savingdir,'/','list_results_eval_alpha0_pam'))
}

if(running){
  nclusters <- ncol(list_with_clusters_alpha0.10$dFrameAsv_pam)-1
  list_results <- list()
  
  for( i in 1:nclusters){
    list_results[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      list_with_clusters_alpha0.10$dFrameAsv_pam[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dframe$ASV_composition))
    cat('iteration------:',i,'\n')
  }
  saveRDS(list_results,file = paste0(savingdir,'/','list_results_eval_alpha0.10_pam'))
}

if(running){
  nclusters <- ncol(list_with_clusters_alpha0.25$dFrameAsv_pam)-1
  list_results <- list()
  
  for( i in 1:nclusters){
    list_results[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      list_with_clusters_alpha0.25$dFrameAsv_pam[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dframe$ASV_composition))
    cat('iteration------:',i,'\n')
  }
  saveRDS(list_results,file = paste0(savingdir,'/','list_results_eval_alpha0.25_pam'))
}

if(running){
  nclusters <- ncol(list_with_clusters_alpha0.50$dFrameAsv_pam)-1
  list_results <- list()
  
  for( i in 1:nclusters){
    list_results[[i]] <-  eval_adist_clustering(compositionDF = left_join(
      list_with_clusters_alpha0.50$dFrameAsv_pam[,c(1,i+1)] %>% rename('name'=1,'Cluster'=2),
      dframe$ASV_composition))
    cat('iteration------:',i,'\n')
  }
  saveRDS(list_results,file = paste0(savingdir,'/','list_results_eval_alpha0.50_pam'))
}

############# Finding K -----------------------------------------------------------------------------------------
running = T
if (running){
  set.seed(1234)
  alpha0_bootstrap = finding_optimal_k_parallel(
    dframe_asv = dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
    B = 500,split_pct = 0.5,maxnclust_ = 50,vec_functions_fromEnv = funcitons_parallel)
  saveRDS(alpha0_bootstrap,file = paste0(savingdir,'/','alpha0_bootstrap'))
}

running = T
if (running){
  set.seed(1234)
  alpha0_bootstrap = finding_optimal_k_parallel(
    dframe_asv = dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
    B = 500,split_pct = 0.5,maxnclust_ = 50,vec_functions_fromEnv = funcitons_parallel,
    inducedDistMatrix = inducedDist_$normalizedMyDist,alpha_ = 0.10)
  saveRDS(alpha0_bootstrap,file = paste0(savingdir,'/','alpha0.10_bootstrap'))
}


running = T
if (running){
  set.seed(1234)
  alpha0_bootstrap = finding_optimal_k_parallel(
    dframe_asv = dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
    B = 500,split_pct = 0.5,maxnclust_ = 50,vec_functions_fromEnv = funcitons_parallel,
    inducedDistMatrix = inducedDist_$normalizedMyDist,alpha_ = 0.25)
  saveRDS(alpha0_bootstrap,file = paste0(savingdir,'/','alpha0.25_bootstrap'))
}

running = T
if (running){
  set.seed(1234)
  alpha0_bootstrap = finding_optimal_k_parallel(
    dframe_asv = dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
    B = 500,split_pct = 0.5,maxnclust_ = 50,vec_functions_fromEnv = funcitons_parallel,
    inducedDistMatrix = inducedDist_$normalizedMyDist,alpha_ = 0.50)
  saveRDS(alpha0_bootstrap,file = paste0(savingdir,'/','alpha0.50_bootstrap'))
}