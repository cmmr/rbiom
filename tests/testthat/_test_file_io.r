# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Biom File Read/Write
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# context("Biom File Read/Write")
# 
# biom <- readRDS("inputs/biom.rds")
# biom$info$id <- "full biom"
# 
# 
# #----------------------------------------------------------
# # Edge cases
# #----------------------------------------------------------
# 
# onesample <- biom
# onesample[['counts']] <- onesample[['counts']][,1]
# onesample$info$id <- "onesample"
# onesample <- repair(onesample)
# 
# onetaxa <- biom
# onetaxa[['counts']] <- onetaxa[['counts']][1,]
# onetaxa$info$id <- "onetaxa"
# onetaxa <- repair(onetaxa)
# 
# 
# lapply(list(biom, onesample, onetaxa), function (x) {
#   
#   
#   #----------------------------------------------------------
#   # Write/Read Classic Tabular Format
#   #----------------------------------------------------------
#   
#   test_that(paste("TSV Intermediate for", x$info$id), {
#     
#     fp <- tempfile(fileext = ".gz")
#     write_biom(x, fp, format="tab")
#     tsv <- read_biom(fp)
#     unlink(fp)
#     
#     expect_equivalent(as.matrix(x$counts), as.matrix(tsv$counts))
#     expect_equivalent(x$taxonomy,          tsv$taxonomy)
#     expect_equivalent(x$shape,             tsv$shape)
#   })
#   
#   
#   #----------------------------------------------------------
#   # Write/Read BIOM Format 1.0 (JSON)
#   #----------------------------------------------------------
#   
#   test_that(paste("JSON Intermediate for", x$info$id), {
#     
#     fp <- tempfile(fileext = ".bz2")
#     write_biom(x, fp, format="json")
#     json <- read_biom(fp)
#     unlink(fp)
#     
#     expect_equivalent(x$counts,    json$counts)
#     expect_equivalent(x$phylogeny, json$phylogeny)
#     expect_equivalent(x$metadata,  json$metadata)
#     expect_equivalent(x$taxonomy,  json$taxonomy)
#     expect_equivalent(x$shape,     json$shape)
#   })
#   
#   
#   #----------------------------------------------------------
#   # Write/Read BIOM Format 2.1 (HDFS)
#   #----------------------------------------------------------
#   if (requireNamespace("rhdf5", quietly = TRUE)) {
#     
#     test_that(paste("HDF5 Intermediate for", x$info$id), {
#       
#       fp <- tempfile()
#       write_biom(x, fp, format="hdf5")
#       hdf5 <- read_biom(fp)
#       unlink(fp)
#       
#       expect_equivalent(x$counts,    hdf5$counts)
#       expect_equivalent(x$phylogeny, hdf5$phylogeny)
#       expect_equivalent(x$metadata,  hdf5$metadata)
#       expect_equivalent(x$taxonomy,  hdf5$taxonomy)
#       expect_equivalent(x$shape,     hdf5$shape)
#     })
#     
#   }
# })
# 
# 
# remove("biom", "onesample", "onetaxa")
