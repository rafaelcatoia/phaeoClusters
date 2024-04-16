reorder_to_matrix_plot<-function(stable_obj,nsubclust=22){
  
  pallete10 <- colorRampPalette(c('white','purple4'))(4)
  
  hclust_stable_df_order <- data.frame(
  ASV = rownames(stable_obj$hclust_),
  subclust = cutree(hclust(d = as.dist(stable_obj_alpha0_n22$hclust_)),k=nsubclust)
  )
  hclust_stable_df_order = hclust_stable_df_order %>% arrange(subclust)
  
  dm_hclust = ggdendroplot::hmReady(df = stable_obj$hclust_) %>% 
    mutate(simm=cut(value,ordered_result = T,breaks = seq(0,1,0.25),include.lowest = T,right = T))
  
  dm_pam = ggdendroplot::hmReady(df = stable_obj$pam_) %>% 
    mutate(simm=cut(value,ordered_result = T,breaks = seq(0,1,0.25),include.lowest = T,right = T))
  
  dm_hclust$rowid=factor(dm_hclust$rowid,levels = hclust_stable_df_order$ASV,ordered = T)
  dm_hclust$variable=factor(dm_hclust$variable,levels = hclust_stable_df_order$ASV,ordered = T)
  
  dm_pam$rowid=factor(dm_pam$rowid,levels = hclust_stable_df_order$ASV,ordered = T)
  dm_pam$variable=factor(dm_pam$variable,levels = hclust_stable_df_order$ASV,ordered = T)
  
  p_hclust <- dm_hclust %>% ggplot(aes(x=rowid, y=variable, fill=simm)) + 
    geom_tile() + 
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())+
    scale_fill_manual(values = pallete10)
  
  p_pam <- dm_pam %>% ggplot(aes(x=rowid, y=variable, fill=simm)) + 
    geom_tile() + 
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())+
    scale_fill_manual(values = pallete10)
  
  return(list(p_pam=p_pam,p_hclust=p_hclust))
}
