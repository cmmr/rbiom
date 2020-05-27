
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
unlink(fp)

test_that("TSV Intermediate", {
  expect_equivalent(as.matrix(biom$counts), as.matrix(tsv$counts))
  expect_equivalent(biom$taxonomy, tsv$taxonomy)
})


#----------------------------------------------------------
# Write/Read BIOM Format 1.0 (JSON)
#----------------------------------------------------------
write.biom(biom, fp, format="json")
json <- read.biom(fp)
unlink(fp)

test_that("JSON Intermediate", {
  expect_equivalent(biom$counts,    json$counts)
  expect_equivalent(biom$phylogeny, json$phylogeny)
  expect_equivalent(biom$metadata,  json$metadata)
  expect_equivalent(biom$taxonomy,  json$taxonomy)
})


#----------------------------------------------------------
# Write/Read BIOM Format 2.1 (HDFS)
#----------------------------------------------------------
if (requireNamespace("rhdf5", quietly = TRUE)) {
  
  write.biom(biom, fp, format="hdf5")
  hdf5 <- read.biom(fp)
  unlink(fp)
  
  test_that("HDF5 Intermediate", {
    expect_equivalent(biom$counts,    hdf5$counts)
    expect_equivalent(biom$phylogeny, hdf5$phylogeny)
    expect_equivalent(biom$metadata,  hdf5$metadata)
    expect_equivalent(biom$taxonomy,  hdf5$taxonomy)
  })
  
  remove("hdf5")
}


remove("biom", "fp", "tsv", "json")
