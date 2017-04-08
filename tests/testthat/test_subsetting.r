
biom <- readRDS("inputs/biom.rds")


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Rarefaction
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Rarefaction")

biom <- rarefy(biom, 100)

test_that("Rarefaction", {
  expect_equal(biom, readRDS("outputs/rarefy.rds"))
})



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Subset/Select
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Subset/Select")


ex1 <- subset(biom, Age > 30)
ex2 <- subset(biom, Body.Site %in% c("Saliva", "Stool"))
ex3 <- subset(biom, Age < 25 & BMI > 22)

test_that("Subsetting", {
  expect_equal(ex1, readRDS("outputs/subset_ex1.rds"))
  expect_equal(ex2, readRDS("outputs/subset_ex2.rds"))
  expect_equal(ex3, readRDS("outputs/subset_ex3.rds"))
})
remove("ex1", "ex2", "ex3")


ex1 <- select(biom, c("1S/~5B.", "DV%sB[?=S", ")p*PCn[="))
ex2 <- select(biom, c(32, 11, 28, 16, 46, 5))
ex3 <- select(biom, 1:50 %% 6 == 0)

test_that("Selecting", {
  expect_equal(ex1, readRDS("outputs/select_ex1.rds"))
  expect_equal(ex2, readRDS("outputs/select_ex2.rds"))
  expect_equal(ex3, readRDS("outputs/select_ex3.rds"))
})
remove("ex1", "ex2", "ex3")




#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Beta Diversity
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

context("Beta Diversity")

dm <- readRDS("outputs/beta.div.rds")

for (metric in c("manhattan", "euclidean", "bray", "jaccard", "unifrac")) {
  test_that(metric, {
    expect_equal(dm[['wdm']][[metric]], beta.div(biom, metric, weighted=TRUE))
    expect_equal(dm[['udm']][[metric]], beta.div(biom, metric, weighted=FALSE))
  })
}

remove("dm", "metric")
