# 
# length(testing)
# 
# hclust_list <- plyr::llply(testing, function(el){el$hclust})
# pam_list <- plyr::llply(testing, function(el){el$pam})
# 
# test<-Reduce('+', hclust_list)
# test_pam <-Reduce('+', pam_list)
# 
# test = test/250
# test[1:10,1:10]
# 
# test_pam <-Reduce('+', pam_list)
# test_pam = test_pam/250
# 
# library(factoextra)
# 
# factoextra::fviz_dist(as.dist(test),
#                       show_labels = FALSE,
#                       gradient = list(low = "white", high = "#00AFBB"))
# 
# factoextra::fviz_dist(as.dist(test_pam),
#                       show_labels = FALSE,
#                       gradient = list(low = "white", high = "#00AFBB"))
# 
# 
# test_pam[1:10,1:10]
# test[1:10,1:10]
# 
# 
# library(gplots)
# install.packages('gplots')
# library(gplots)
# heatmap.2(my.image, density.info="none", trace="none", dendrogram='none', 
#           Rowv=FALSE, Colv=FALSE)
# 
# install.packages('gclus')
# plot1 <- coldiss(D = as.dist(test_pam),nc = 10)
# 
# # coldiss()
# # Color plots of a dissimilarity matrix, without and with ordering
# #
# # License: GPL-2 
# # Author: Francois Gillet, August 2009
# #
# 
# "coldiss" <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
# {
#   require(gclus)
#   
#   if (max(D)>1) D <- D/max(D)
#   
#   if (byrank) {
#     spe.color = dmat.color(1-D, cm.colors(nc))
#   }
#   else {
#     spe.color = dmat.color(1-D, byrank=FALSE, cm.colors(nc))
#   }
#   
#   spe.o = order.single(1-D)
#   speo.color = spe.color[spe.o,spe.o]
#   
#   op = par(mfrow=c(1,2), pty="s")
#   
#   if (diag) {
#     plotcolors(spe.color, rlabels=attributes(D)$Labels, 
#                main="Dissimilarity Matrix", 
#                dlabels=attributes(D)$Labels)
#     plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
#                main="Ordered Dissimilarity Matrix", 
#                dlabels=attributes(D)$Labels[spe.o])
#   }
#   else {
#     plotcolors(spe.color, rlabels=attributes(D)$Labels, 
#                main="Dissimilarity Matrix")
#     plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o], 
#                main="Ordered Dissimilarity Matrix")
#   }
#   
#   par(op)
# }
# 
# 
# 
# 
# 