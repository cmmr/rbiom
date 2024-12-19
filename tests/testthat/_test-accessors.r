# 
# test_that("biom_info", {
#   expect_snapshot(biom_info(hmp50)) })
# 
# test_that("sample_sums", {
#   expect_snapshot(sample_sums(hmp5))
#   expect_snapshot(sample_sums(hmp5, long=TRUE))
#   expect_snapshot(sample_sums(hmp5, long=TRUE, md=TRUE))
#   expect_snapshot(sample_sums(min5, long=TRUE, md=TRUE)) })
# 
# test_that("rarefaction_level", {
#   expect_snapshot(rarefaction_level(hmp5))
#   expect_snapshot(rarefaction_level(rare5))
#   expect_snapshot(rarefaction_level(as_percent(rare5))) })
# 
# test_that("has_tree", {
#   expect_snapshot(has_tree(hmp5))
#   expect_snapshot(has_tree(min5)) })
# 
# test_that("has_sequences", {
#   expect_snapshot(has_sequences(hmp5))
#   expect_snapshot(has_sequences(min5)) })
# 
# test_that("is_rarefied", {
#   expect_snapshot(is_rarefied(hmp5))
#   expect_snapshot(is_rarefied(rare5)) })
# 
# test_that("n_samples", {
#   expect_snapshot(n_samples(hmp5)) })
# 
# test_that("n_otus", {
#   expect_snapshot(n_otus(hmp5)) })
# 
# test_that("n_ranks", {
#   expect_snapshot(n_ranks(hmp5))
#   expect_snapshot(n_ranks(min5)) })
