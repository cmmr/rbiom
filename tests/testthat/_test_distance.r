# 
# biom <- readRDS("inputs/biom.rds")
# 
# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Distance Matrices
# #   'prev' values confirmed correct using phyloseq for 
# #   all and QIIME for unifrac.
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# context("Beta Diversity")
# 
# prev <- readRDS("outputs/bdiv_distmat.rds")
# 
# for (method in c("manhattan", "euclidean", "bray-curtis", "jaccard", "unifrac")) {
#   
#   wdm <- bdiv_distmat(biom=biom, method=method, weighted=TRUE)
#   udm <- bdiv_distmat(biom=biom, method=method, weighted=FALSE)
#   
#   test_that(method, {
#     expect_equal(as.vector(wdm), as.vector(prev[[method]][['wdm']]))
#     expect_equal(as.vector(udm), as.vector(prev[[method]][['udm']]))
#   })
# }
# 
# remove("biom", "method", "prev", "udm", "wdm")
