
biom <- readRDS("inputs/biom.rds")
biom <- rarefy(biom, 100)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Alpha Diversity
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Alpha Diversity")

curr <- alpha.div(biom)
prev <- readRDS("outputs/alpha.div.rds")

test_that("Diversity Metrics", {
  expect_equal(curr[,'OTUs'],    prev[,'OTUs'])
  expect_equal(curr[,'Shannon'], prev[,'Shannon'])
  expect_equal(curr[,'Simpson'], prev[,'Simpson'])
  expect_equal(curr[,'Chao1'],   prev[,'Chao1'])
})

remove("curr", "prev")



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Taxa Rollup
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Taxa Rollup")

tr <- readRDS("outputs/taxa.rollup.rds")

for (txRank in c(colnames(biom$taxonomy), 'OTU')) {
  test_that(txRank, {
    expect_equal(tr[[txRank]], taxa.rollup(biom, txRank))
  })
}

remove("tr", "txRank")






