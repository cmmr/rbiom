library(rbiom)
library(phyloseq)
library(doParallel)

infile <- system.file("extdata", "hmp500.biom", package = "rbiom")

biom <- read.biom(infile)

phy  <- import_biom(infile)
phy_tree(phy) <- biom$phylogeny

cl <- makeCluster(8)
registerDoParallel(cl)




t1 <- Sys.time()
x <- unifrac(biom, weighted=TRUE)
Sys.time() - t1

t1 <- Sys.time()
x <- UniFrac(phy, weighted=TRUE, normalized=FALSE, parallel=TRUE)
Sys.time() - t1



res <- benchmark(unifrac(biom, weighted=TRUE),
                 UniFrac(phy, weighted=TRUE, normalized=FALSE, parallel=TRUE),
                 order="relative")


# rbiom::write.biom(biom, "~/Desktop/hmp500.biom", "json")
# ape::write.tree(biom$phylogeny, "~/Desktop/hmp500.tree")
# time parallel_beta_diversity.py -i hmp500.biom -o beta_div/ -t hmp500.tre -O 48 -m weighted_unifrac
# 1m13.453s

