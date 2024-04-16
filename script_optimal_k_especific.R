######## --- last script 
# setwd('phaeoClusters/')
library(dplyr) ; library(tidyr) ; library(ggplot2) 
library(doParallel) ; library(foreach) ; library(doSNOW)

root <- rprojroot::has_file(".git/index")
datadir = root$find_file("data")
funsdir = root$find_file("functions")
savingdir = root$find_file("saved_files")

datapath = root$find_file(paste0(datadir,'/','grump.phaeocystis_asv_long.csv'))
files_vec <- list.files(funsdir)

for( i in 1:length(files_vec)){
  source(root$find_file(paste0(funsdir,'/',files_vec[i])))
}

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


#dframe_asv <- dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts)


finding_optimal_k_parallel_alt <- function(
    dframe_asv,
    B=10,split_pct = 0.75,
    maxnclust_ = 25,
    vec_functions_fromEnv,alpha_=NULL,inducedDistMatrix=NULL){
  
  #B is the number of splits
  #split_pct is the percentage of samples we will use to fit, therefore (1-0.75) will be used to evaluate that.
  
  ## calculating the sample compositions
  df_comp_all <- dframe_asv %>% 
    pivot_wider(id_cols = c('SampleKey'),
                values_from ='Raw.Sequence.Counts',
                names_from='ASV_name',values_fill = 0) %>%
    arrange(SampleKey) %>% 
    mutate(across(where(is.numeric))/
             rowSums(across(where(is.numeric)))) %>% ### Here we have the compositions on the samples
    mutate(across(where(is.numeric),function(x) {x+0.00000000001}))%>%  ### Adding a small value in all, in order to have the CLR transformation / aitichison distance
    mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% data.frame()
  
  ## selecting the sample to train
  nsamples <- nrow(df_comp_all)
  vet_idx <- 1:nsamples
  
  vet_asv_names <- dframe_asv %>% select(ASV_name) %>% distinct() %>% pull()
  vet_asv_names <- vet_asv_names[order(vet_asv_names)]
  
  ## now spplitting::::::::::::::::::::::::: 
  cl <- makeCluster(detectCores()-1)
  registerDoSNOW(cl)
  pb=txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  t0 <- Sys.time()
  
  results <- foreach(i = 1:B,
                     .packages=c("dplyr", "cluster", "dplyr", "robCompositions", "tidyr","knitr","kableExtra"),
                     .options.snow = opts,
                     .export = vec_functions_fromEnv) %dopar% {
                       
                       sample_train_idx <- sample(vet_idx,size = round(nsamples*split_pct))
                       sample_train_idx <- sample_train_idx[order(sample_train_idx)]
                       
                       comp_fit <- df_comp_all[sample_train_idx,] %>%
                         pivot_longer(cols = -'SampleKey') %>%
                         pivot_wider(names_from='SampleKey',values_from = 'value') %>% 
                         mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
                         arrange(name)
                       #comp_fit[1:5,1:5]
                       
                       comp_eval <- df_comp_all[-sample_train_idx,] %>%
                         pivot_longer(cols = -'SampleKey') %>%
                         pivot_wider(names_from='SampleKey',values_from = 'value') %>% 
                         mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
                         arrange(name)
                       #comp_eval[1:5,1:5]
                       
                       
                       ##### clustering
                       aitDist_fit <- robCompositions::aDist(x = comp_fit[,-1])
                       norm_fit <- aitDist_fit/norm(as.matrix(aitDist_fit),type = '2')
                       
                       if(!is.null(inducedDistMatrix)){
                         norm_fit1 <- (1-alpha_[1]) * as.matrix(norm_fit) + alpha_[1]*as.matrix(inducedDistMatrix)
                         norm_fit2 <- (1-alpha_[2]) * as.matrix(norm_fit) + alpha_[2]*as.matrix(inducedDistMatrix)
                         norm_fit3 <- (1-alpha_[3]) * as.matrix(norm_fit) + alpha_[3]*as.matrix(inducedDistMatrix)
                       }
                       
                       obj_clusterize <- clusterize_ward_pam(as.dist(norm_fit),dFrameAsv =vet_asv_names, maxnclust = maxnclust_)
                       obj_clusterize1 <- clusterize_ward_pam(as.dist(norm_fit1),dFrameAsv =vet_asv_names, maxnclust = maxnclust_)
                       obj_clusterize2 <- clusterize_ward_pam(as.dist(norm_fit2),dFrameAsv =vet_asv_names, maxnclust = maxnclust_)
                       obj_clusterize3 <- clusterize_ward_pam(as.dist(norm_fit3),dFrameAsv =vet_asv_names, maxnclust = maxnclust_)
                       
                       #### evaluating cluster
                       list_results_hclust <- list()
                       list_results_pam <- list()
                       
                       list_results_hclust1 <- list()
                       list_results_pam1 <- list()
                       
                       list_results_hclust2 <- list()
                       list_results_pam2 <- list()
                       
                       list_results_hclust3 <- list()
                       list_results_pam3 <- list()
                       
                       for( j in 1:maxnclust_){
                         ############ ---------------------------------------------------------------
                         list_results_hclust[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize$dFrameAsv_hclust[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         list_results_pam[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize$dFrameAsv_pam[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         ############ -----------------------------------------------------------------
                         
                         list_results_hclust1[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize1$dFrameAsv_hclust[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         list_results_pam1[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize1$dFrameAsv_pam[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         list_results_hclust2[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize2$dFrameAsv_hclust[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         list_results_pam2[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize2$dFrameAsv_pam[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         list_results_hclust3[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize3$dFrameAsv_hclust[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         list_results_pam3[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize3$dFrameAsv_pam[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         cat('iteration------:',j,'\n')
                       }
                       
                       list_df_results <- bind_rows(
                         summarise_adist_clustering(list_results_hclust) %>% mutate(method='hclustD',alpha0=0,iteration = i),
                         summarise_adist_clustering(list_results_pam) %>% mutate(method='pam',alpha0=0,iteration = i),
                         summarise_adist_clustering(list_results_hclust1) %>% mutate(method='hclustD',alpha0=0.1,iteration = i),
                         summarise_adist_clustering(list_results_pam1) %>% mutate(method='pam',alpha0=0.1,iteration = i),
                         summarise_adist_clustering(list_results_hclust2) %>% mutate(method='hclustD',alpha0=0.25,iteration = i),
                         summarise_adist_clustering(list_results_pam2) %>% mutate(method='pam',alpha0=0.25,iteration = i),
                         summarise_adist_clustering(list_results_hclust3) %>% mutate(method='hclustD',alpha0=0.50,iteration = i),
                         summarise_adist_clustering(list_results_pam3) %>% mutate(method='pam',alpha0=0.50,iteration = i)
                       )
                       rm(list_results_hclust)
                       rm(list_results_pam)
                       rm(list_results_hclust1)
                       rm(list_results_pam1)
                       rm(list_results_hclust2)
                       rm(list_results_pam2)
                       rm(list_results_hclust3)
                       rm(list_results_pam3)
                       list_df_results
                     }
  stopCluster(cl)
  rm(cl)
  return(results)
}

############# Finding K -----------------------------------------------------------------------------------------

set.seed(1234)
cv_clustering = finding_optimal_k_parallel_alt(
  dframe_asv = dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
  B = 100,split_pct = 0.50,maxnclust_ = 50,vec_functions_fromEnv = funcitons_parallel,
  alpha_ = c(0.1,0.25,0.5),inducedDistMatrix = inducedDist_$normalizedMyDist)
saveRDS(cv_clustering,file = paste0(savingdir,'/','cv_clustering'))


