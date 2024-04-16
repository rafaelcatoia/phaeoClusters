plot_knowns_vs_unknowns <- function(dm_aux,df_ordered_ASVs,specific_species=F){
  
  pallete10 <- colorRampPalette(c('white','purple4'))(4)
  
  df_plot <- dm_aux %>% filter(variable%in%df_ordered_ASVs$ASV_name)
  df_plot$variable <- factor(df_plot$variable, levels = df_ordered_ASVs$ASV_name, ordered = T)
  
  rowid_order <- df_plot %>%
    transmute(rowid=as.character(rowid)) %>%
    filter(rowid %!in%df_ordered_ASVs$ASV_name) %>% distinct() %>% pull()
  
  rowid_order = c(df_ordered_ASVs$ASV_name,rowid_order)
  
  df_plot$rowid <- factor(df_plot$rowid, levels = rowid_order, ordered = T)
  
  df_plot = df_plot %>%
    left_join(df_ordered_ASVs %>% rename('variable'=1,'species1'=2))
  
  vet_limits <- df_ordered_ASVs %>% group_by(Species) %>% 
    summarise(Freq=n()) %>% mutate(Cumm_Freq = cumsum(Freq)) %>% 
    select(Cumm_Freq) %>% pull()
  
  vet_limits= vet_limits+0.5
  if(!specific_species){
  plot_out <- df_plot %>% 
    ggplot(aes(x=rowid, y=variable, fill=simm)) + 
    geom_tile() + 
    theme_minimal() +
    scale_fill_manual(values = pallete10)+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())+
    xlab('ASVs')+
    ylab('ASVs') + 
    geom_hline(yintercept = vet_limits)+
    geom_vline(xintercept = vet_limits[length(vet_limits)])
    #geom_abline(slope = 1,intercept = 0)+
    theme(legend.position = 'bottom')
  }
  
  if(specific_species){
    plot_out <- df_plot %>% 
      ggplot(aes(x=rowid, y=variable, fill=simm)) + 
      geom_tile() + 
      theme_minimal() +
      scale_fill_manual(values = pallete10)+
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
      xlab('ASVs')+
      ylab('ASVs') + 
      geom_vline(xintercept = vet_limits[length(vet_limits)])
    #geom_abline(slope = 1,intercept = 0)+
    theme(legend.position = 'bottom')
  }
  
  return(plot_out)
}
