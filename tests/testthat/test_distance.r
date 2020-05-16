
biom <- readRDS("inputs/biom.rds")


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Distance Matrices
#   'prev' values confirmed correct using phyloseq for 
#   all and QIIME for unifrac.
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Beta Diversity")

prev <- readRDS("outputs/beta.div.rds")

for (method in c("manhattan", "euclidean", "bray-curtis", "jaccard", "unifrac")) {
  
  wdm <- beta.div(biom=biom, method=method, weighted=TRUE)
  udm <- beta.div(biom=biom, method=method, weighted=FALSE)
  
  test_that(method, {
    expect_equal(as.vector(wdm), as.vector(prev[[method]][['wdm']]))
    expect_equal(as.vector(udm), as.vector(prev[[method]][['udm']]))
  })
}

remove("biom", "method", "prev", "udm", "wdm")
