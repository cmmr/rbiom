
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# All the little accessors
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Accessors")

biom <- readRDS("inputs/biom.rds")


test_that("sample_names", expect_equal(sample_names(biom), colnames(biom[['counts']])))
test_that("taxa_names",   expect_equal(taxa_names(biom),   rownames(biom[['counts']])))
test_that("taxa_ranks",   expect_equal(taxa_ranks(biom),   c('OTU', colnames(biom[['taxonomy']]))))
test_that("counts",       expect_equal(counts(biom),       as.matrix(biom[['counts']])))
test_that("taxonomy",     expect_equal(taxonomy(biom),     biom[['taxonomy']]))
test_that("phylogeny",    expect_equal(phylogeny(biom),    biom[['phylogeny']]))
test_that("metadata",     expect_equal(metadata(biom),     biom[['metadata']]))
test_that("info",         expect_equal(info(biom),         biom[['info']]))


remove("biom")
