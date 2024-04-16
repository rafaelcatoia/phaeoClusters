#dframe_asv <- dframe$dframe %>% select(SampleKey,Raw.Sequence.Counts,ASV_name)

asv_stable_calc_notparallel <- function(
    dframe_asv,
    B=10,
    split_pct = 0.75,
    nclust=10,
    splittingInducerMatrix=NULL,alpha_=0.1){
  
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
  
  output_hclust <- matrix(0,nrow=length(vet_asv_names),ncol=length(vet_asv_names))
  output_pam <- matrix(0,nrow=length(vet_asv_names),ncol=length(vet_asv_names))
  for ( i in 1:B){

    cat('Iteration-------------- ', i , 'of ' , B, "-----", paste0(Sys.time()),'\n')
    
    
    sample_train_idx <- sample(vet_idx,size = round(nsamples*split_pct))
    sample_train_idx <- sample_train_idx[order(sample_train_idx)]
    
    comp_fit <- df_comp_all[sample_train_idx,] %>%
      pivot_longer(cols = -'SampleKey') %>%
      pivot_wider(names_from='SampleKey',values_from = 'value') %>% 
      mutate(across(where(is.numeric))/rowSums(across(where(is.numeric)))) %>% 
      arrange(name)
    
    aitDist_fit <- robCompositions::aDist(x = comp_fit[,-1])
    norm_fit <- aitDist_fit/norm(as.matrix(aitDist_fit),type = '2')
    
    ##### creating the mixture of distance matrices 
    if(!is.null(splittingInducerMatrix)){
      norm_fit <- (1-alpha_) * as.matrix(norm_fit) + alpha_*as.matrix(splittingInducerMatrix)
    }
    
    hclust_obj <- cbind.data.frame(ASV = vet_asv_names,
                                   cluster_id=cutree(hclust(d = as.dist(norm_fit) ,method = 'ward.D'),k=nclust))
    
    pam_obj <- cbind.data.frame(ASV = vet_asv_names,
                                cluster_id = cluster::pam(x = as.dist(norm_fit),k = nclust, cluster.only=TRUE))
    
    output_hclust =  output_hclust + dist_given_cluster(hclust_obj)
    # output_hclust[1:5,1:5]
    output_pam = output_pam + dist_given_cluster(pam_obj)
    # output_pam[1:5,1:5]
  }
  return(
    list(hclust_ = output_hclust/B,pam_=output_pam/B)
  )
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



