# 
# library(slam)
# 
# biom <- readRDS("inputs/biom.rds")
# 
# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Rarefaction
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# context("Rarefaction")
# 
# biom <- rarefy(biom, 100)
# 
# test_that("Rarefaction", {
#   expect_equal(biom, readRDS("outputs/rarefy.rds"))
#   
#   # matrix in, matrix out
#   x  <- matrix(runif(100, 0, 1000), ncol=10)
#   y  <- rarefy(x, depth=100)
#   z  <- merge(reshape2::melt(x), reshape2::melt(y), by=c('Var1', 'Var2'))
#   r2 <- summary(lm(value.x ~ value.y, data=z))$r.squared
#   expect_is(y, "matrix")
#   expect_gt(r2, .3)
#   
#   # simple_triplet_matrix in, simple_triplet_matrix out
#   x  <- as.simple_triplet_matrix(matrix(runif(100, 0, 1000), ncol=10))
#   y  <- rarefy(x, depth=100)
#   z  <- merge(
#     reshape2::melt(as.matrix(x)), 
#     reshape2::melt(as.matrix(y)), 
#     by=c('Var1', 'Var2'))
#   r2 <- summary(lm(value.x ~ value.y, data=z))$r.squared
#   expect_is(y, "simple_triplet_matrix")
#   expect_gt(r2, .3)
#   
#   remove("x", "y", "z", "r2")
#   # ggplot(z, aes(x=value.x, y=value.y)) + geom_smooth(method = "lm") + geom_point()
# })
# 
# 
# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Subset/Select
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# context("Subset/Select")
# 
# 
# ex1 <- subset(biom, Age > 30)
# ex2 <- subset(biom, Body.Site %in% c("Saliva", "Stool"))
# ex3 <- subset(biom, Age < 25 & BMI > 22)
# 
# test_that("Subsetting", {
#   expect_equal(ex1, readRDS("outputs/subset_ex1.rds"))
#   expect_equal(ex2, readRDS("outputs/subset_ex2.rds"))
#   expect_equal(ex3, readRDS("outputs/subset_ex3.rds"))
# })
# remove("ex1", "ex2", "ex3")
# 
# 
# ex1 <- select(biom, c("1S/~5B.", "DV%sB[?=S", ")p*PCn[="))
# ex2 <- select(biom, c(32, 11, 28, 16, 46, 5))
# ex3 <- select(biom, 1:50 %% 6 == 0)
# 
# test_that("Selecting", {
#   expect_equal(ex1, readRDS("outputs/select_ex1.rds"))
#   expect_equal(ex2, readRDS("outputs/select_ex2.rds"))
#   expect_equal(ex3, readRDS("outputs/select_ex3.rds"))
# })
# remove("ex1", "ex2", "ex3")
# 
# 
# remove("biom")
#   
