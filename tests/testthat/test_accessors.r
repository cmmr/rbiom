
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# All the little accessor functions
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Accessors")

biom <- readRDS("inputs/biom.rds")


test_that("sample.names", expect_equal(sample.names(biom), colnames(biom[['counts']])))
test_that("taxa.names",   expect_equal(taxa.names(biom),   rownames(biom[['counts']])))
test_that("taxa.ranks",   expect_equal(taxa.ranks(biom),   colnames(biom[['taxonomy']])))
test_that("counts",       expect_equal(counts(biom),       as.matrix(biom[['counts']])))
test_that("taxonomy",     expect_equal(taxonomy(biom),     biom[['taxonomy']]))
test_that("phylogeny",    expect_equal(phylogeny(biom),    biom[['phylogeny']]))
test_that("metadata",     expect_equal(metadata(biom),     biom[['metadata']]))
test_that("info",         expect_equal(info(biom),         biom[['info']]))


remove("biom")
