summarise_eval_clusters_distMatrix <- function(list_results_eval){
  list_hclust = list_results_eval$hclust
  list_pam = list_results_eval$pam
 
  df_out <-bind_rows(
  list_hclust %>% plyr::ldply(function(el){
    data.frame(
      nclust = length(el$vet_number_obs_per_clusters),
      weighted_mean_avg_within_to_medoid = el$weighted_mean_avg_within_to_medoid,
      avg_dist_between_medoid = el$avg_dist_between_medoid,
      ratio_avg_within_between_medoids = el$ratio_avg_within_between_medoids,
      max_max_dist_within = el$max_max_dist_within,
      avg_max_dist_within = el$avg_max_dist_within,
      sum_weighted_max_dist_within = el$sum_weighted_max_dist_within,
      min_max_dist_within = el$min_max_dist_within,
      avg_min_dist_between_clusters = el$avg_min_dist_between_clusters,
      sum_min_dist_between_clusters = el$sum_min_dist_between_clusters,
      within_between_avg = el$within_between_avg,
      numer_within_between_avg = el$numer_within_between_avg,
      denom_within_between_avg = el$denom_within_between_avg,
      method='hclust'
    )
  }) %>% return() , 
  
  list_pam %>% plyr::ldply(function(el){
    data.frame(
      nclust = length(el$vet_number_obs_per_clusters),
      weighted_mean_avg_within_to_medoid = el$weighted_mean_avg_within_to_medoid,
      avg_dist_between_medoid = el$avg_dist_between_medoid,
      ratio_avg_within_between_medoids = el$ratio_avg_within_between_medoids,
      max_max_dist_within = el$max_max_dist_within,
      avg_max_dist_within = el$avg_max_dist_within,
      sum_weighted_max_dist_within = el$sum_weighted_max_dist_within,
      min_max_dist_within = el$min_max_dist_within,
      avg_min_dist_between_clusters = el$avg_min_dist_between_clusters,
      sum_min_dist_between_clusters = el$sum_min_dist_between_clusters,
      within_between_avg = el$within_between_avg,
      numer_within_between_avg = el$numer_within_between_avg,
      denom_within_between_avg = el$denom_within_between_avg,
      method='pam'
    )
  }) %>%
    return()
  )
  return(df_out)
}
