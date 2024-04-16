plot_knowns_matrix <- function(dm_aux,df_ordered_ASVs){
  
  pallete10 <- colorRampPalette(c('white','purple4'))(4)
  
  df_plot <- dm_aux %>% filter(rowid%in%df_ordered_ASVs$ASV_name,
                 variable%in%df_ordered_ASVs$ASV_name)
  df_plot$rowid <- factor(df_plot$rowid, levels = df_ordered_ASVs$ASV_name, ordered = T)
  df_plot$variable <- factor(df_plot$variable, levels = df_ordered_ASVs$ASV_name, ordered = T)
  
  df_plot = df_plot %>%
    left_join(df_ordered_ASVs %>% rename('rowid'=1,'species1'=2)) %>% 
    left_join(df_ordered_ASVs %>% rename('variable'=1,'species2'=2))
  
  vet_limits <- df_ordered_ASVs %>% group_by(Species) %>% 
    summarise(Freq=n()) %>% mutate(Cumm_Freq = cumsum(Freq)) %>% 
    select(Cumm_Freq) %>% pull()
  
  vet_limits= vet_limits+0.5
  
  plot_out <- df_plot %>% 
    ggplot(aes(x=rowid, y=variable, fill=simm)) + 
    geom_tile() + 
    theme_minimal() +
    scale_fill_manual(values = pallete10)+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())+
    geom_vline(xintercept = vet_limits)+
    geom_hline(yintercept = vet_limits)+
    xlab('ASVs')+
    ylab('ASVs')
  
  return(plot_out)
}
