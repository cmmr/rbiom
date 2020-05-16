
biom <- readRDS("inputs/biom.rds")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Alpha Diversity
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Alpha Diversity")

curr <- alpha.div(biom)
prev <- readRDS("outputs/alpha.div.rds")

test_that("Diversity Metrics", {
  expect_equal(tolerance = 0.0001, curr[['OTUs']],    prev[['OTUs']])
  expect_equal(tolerance = 0.0001, curr[['Shannon']], prev[['Shannon']])
  expect_equal(tolerance = 0.0001, curr[['Simpson']], prev[['Simpson']])
  expect_equal(tolerance = 0.0001, curr[['Chao1']],   prev[['Chao1']])
})

remove("curr", "prev")



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Taxa Rollup
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Taxa Rollup")

tr <- readRDS("outputs/taxa.rollup.rds")

for (txRank in c(taxa.ranks(biom), 'OTU')) {
  test_that(txRank, {
    expect_equal(tr[[txRank]], taxa.rollup(biom, txRank))
  })
}

remove("tr", "txRank")

remove("biom")
