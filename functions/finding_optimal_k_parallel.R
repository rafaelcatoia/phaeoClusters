#dframe_asv <- dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts)


finding_optimal_k_parallel <- function(
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
                       sample_train_idx[order(sample_train_idx)]
                       
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
                         norm_fit <- (1-alpha_) * as.matrix(aitDist_fit) + alpha_*as.matrix(inducedDistMatrix)
                       }
                       
                       obj_clusterize <- clusterize_ward_pam(norm_fit,dFrameAsv =vet_asv_names, maxnclust = maxnclust_)
                       
                       
                       
                       
                       #### evaluating cluster
                       list_results_hclust <- list()
                       list_results_pam <- list()
                       for( j in 1:maxnclust_){
                         list_results_hclust[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize$dFrameAsv_hclust[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         list_results_pam[[j]] <-  eval_adist_clustering(compositionDF = left_join(
                           obj_clusterize$dFrameAsv_pam[,c(1,j+1)] %>% rename('name'=1,'Cluster'=2),
                           comp_eval))
                         
                         
                         cat('iteration------:',j,'\n')
                       }
                       
                       list_df_results <- bind_rows(
                         summarise_adist_clustering(list_results_hclust) %>% mutate(method='hclustD',iteration = i),
                         summarise_adist_clustering(list_results_pam) %>% mutate(method='pam',iteration = i)
                       )
                       rm(list_results_hclust)
                       rm(list_results_pam)
                       list_df_results
                     }
  stopCluster(cl)
  rm(cl)
  return(results)
}
