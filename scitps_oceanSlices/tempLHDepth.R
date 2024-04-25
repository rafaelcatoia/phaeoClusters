library(dplyr) ; library(tidyr) ; library(ggplot2) 
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
saveRDS(obj_eval$list_results_hclust,file = paste0(savingdirOceanSlices,'/','list_LH_Depth_eval_alpha0.1_hclust'))
saveRDS(obj_eval$list_results_medoid,file = paste0(savingdirOceanSlices,'/','list_LH_Depth_eval_alpha0.1_pam'))