asv_stable_calc <- function(
    dframe_asv,
    B=10,
    split_pct = 0.75,
    nclust=10,
    vec_functions_fromEnv,
    splittingInducerMatrix=NULL){
  
  #### Function that will probably be used only here.
  dist_given_cluster <- function(df_ASv_ClusterMembership){
    out <- as.matrix(dist(df_ASv_ClusterMembership$cluster_id))
    out <- ifelse(out == 0, 1, 0)
    rownames(out) <- df_ASv_ClusterMembership$ASV
    colnames(out) <- df_ASv_ClusterMembership$ASV
    return(out)
  }
  
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
  
  ## preparing for parallelization
  cl <- makeCluster(detectCores()-1)
  registerDoSNOW(cl)
  pb=txtProgressBar(max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  t0 <- Sys.time()
  results <- foreach(
    i = 1:B,
    .packages=c("dplyr", "cluster", "dplyr", "robCompositions", "tidyr","knitr","kableExtra"),
    .options.snow = opts,
    .export = vec_functions_fromEnv) %dopar% {
      
      sample_train_idx <- sample(vet_idx,size = round(nsamples*split_pct))
      #sample_train_idx <- sample_train_idx[order(sample_train_idx)]
      
      comp_fit <- df_comp_all[sample_train_idx,] %>%
        pivot_longer(cols = -'SampleKey') %>%
        pivot_wider(names_from='SampleKey',values_from = 'value') %>% 
        mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
        arrange(name)
      #comp_fit[1:5,1:5]
      
      ##### clustering
      aitDist_fit <- robCompositions::aDist(x = comp_fit[,-1])
      norm_fit <- aitDist_fit/norm(as.matrix(aitDist_fit),type = '2')
      
      ##### creating the mixture of distance matrices 
      if(!is.null(splittingInducerMatrix)){
        
      }
      ##### creating the clusters 
      
      hclust_obj <- cbind.data.frame(ASV = vet_asv_names,
                                     cluster_id=cutree(hclust(d = as.dist(aitDist_fit) ,method = 'ward.D'),k=nclust))
      
      pam_obj <- cbind.data.frame(ASV = vet_asv_names,
                                  cluster_id = cluster::pam(x = as.dist(aitDist_fit),k = nclust, cluster.only=TRUE))
      
      hclust_matrix <- dist_given_cluster(hclust_obj)
      
      pam_matrix <- dist_given_cluster(pam_obj)
      
      list_results = list(hclust = hclust_matrix, pam=pam_matrix)
      list_results
    }
  stopCluster(cl)
  rm(cl)
  return(results)
}


####### Testing 

# testing <- asv_stable_calc(
#   dframe_asv=dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
#   B=250,
#   split_pct = 0.5,
#   nclust=10,
#   vec_functions_fromEnv=funcitons_parallel,
#   splittingInducerMatrix=NULL
# )



