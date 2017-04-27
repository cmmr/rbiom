
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Biom File Read/Write
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Biom File Read/Write")

biom <- readRDS("inputs/biom.rds")
fp   <- tempfile()


#----------------------------------------------------------
# Write/Read Classic Tabular Format
#----------------------------------------------------------
write.biom(biom, fp, format="tab")
tsv <- read.biom(fp)

test_that("TSV Intermediate", {
  expect_equal(as.matrix(biom$counts), as.matrix(tsv$counts))
  expect_equal(biom$taxonomy, tsv$taxonomy)
})


#----------------------------------------------------------
# Write/Read BIOM Format 1.0 (JSON)
#----------------------------------------------------------
write.biom(biom, fp, format="json")
json <- read.biom(fp)

test_that("JSON Intermediate", {
  expect_equal(biom$counts,    json$counts)
  expect_equal(biom$phylogeny, json$phylogeny)
  expect_equal(biom$metadata,  json$metadata)
  expect_equal(biom$taxonomy,  json$taxonomy)
})


# #----------------------------------------------------------
# # Write/Read BIOM Format 2.1 (HDFS)
# #----------------------------------------------------------
# write.biom(biom, fp, format="hdf5")
# hdf5 <- read.biom(fp)
# 
# test_that("HDF5 Intermediate", {
#   expect_equal(biom$counts,    hdf5$counts)
#   expect_equal(biom$phylogeny, hdf5$phylogeny)
#   expect_equal(biom$metadata,  hdf5$metadata)
#   expect_equal(biom$taxonomy,  hdf5$taxonomy)
# })
# 
# 
# unlink(fp)
# remove("fp", "tsv", "json", "hdf5")



