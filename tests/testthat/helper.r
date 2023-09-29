
rare50 <- sample_rarefy(hmp50)
min50  <- biom_build(counts = otu_matrix(hmp50))

hmp5  <- sample_select(hmp50, 1:5)
rare5 <- sample_rarefy(hmp5)
min5  <- biom_build(counts = otu_matrix(hmp5))

