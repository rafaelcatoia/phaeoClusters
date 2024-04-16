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

inducedDist_ = inducedDist(
  dFrame = dframe$dframe,
  c1 = 1,c2=1000,c3=10,c4=10,
  compMatrix = dframe$ASV_composition)

#################################################################################
running = T
# if(running){
#   stable_obj_alpha0 <- asv_stable_calc_notparallel(
#     dframe_asv=dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
#     B=500,
#     split_pct = 0.5,
#     nclust=10,
#     splittingInducerMatrix=NULL
#   )
#   saveRDS(object = stable_obj_alpha0,
#           file = paste0(savingdir,'/','stable_obj_alpha0'))
# }

######## --- stable with alpha = 0.1
#if(running){
#  stable_obj_alpha0.10 <- asv_stable_calc_notparallel(
#    dframe_asv=dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
#    B=500,
#    split_pct = 0.5,
#    nclust=10,
#    splittingInducerMatrix=inducedDist_$normalizedMyDist,alpha_ = 0.1)
#  saveRDS(object = stable_obj_alpha0.10,
#          file = paste0(savingdir,'/','stable_obj_alpha0.10'))
#}
# 
# ######## --- stable with alpha = 0.25
#if(running){
#  stable_obj_alpha0.25 <- asv_stable_calc_notparallel(
#    dframe_asv=dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
#    B=500,
#    split_pct = 0.5,
#    nclust=10,
#    splittingInducerMatrix=inducedDist_$normalizedMyDist,alpha_ = 0.25)
#  saveRDS(object = stable_obj_alpha0.25,
#          file = paste0(savingdir,'/','stable_obj_alpha0.25'))
#}
# 
 ######## --- stable with alpha = 0.50
 if(running){
   stable_obj_alpha0.50 <- asv_stable_calc_notparallel(
     dframe_asv=dframe$dframe %>% select(SampleKey,ASV_name,Raw.Sequence.Counts),
     B=500,
     split_pct = 0.8,
     nclust=22,
     splittingInducerMatrix=inducedDist_$normalizedMyDist,alpha_ = 0.50)
   saveRDS(object = stable_obj_alpha0.50,
           file = paste0(savingdir,'/','stable_obj_alpha0.50_n22'))
 }# 