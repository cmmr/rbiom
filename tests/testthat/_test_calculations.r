# 
# biom <- readRDS(test_path("inputs/biom.rds"))
# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Alpha Diversity
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# context("Alpha Diversity")
# 
# curr <- adiv_matrix(biom)
# prev <- readRDS(test_path("outputs/adiv_table.rds"))
# 
# test_that("Diversity Metrics", {
#   expect_equal(tolerance = 0.0001, unname(curr[,'OTUs']),    prev[['OTUs']])
#   expect_equal(tolerance = 0.0001, unname(curr[,'Shannon']), prev[['Shannon']])
#   expect_equal(tolerance = 0.0001, unname(curr[,'Simpson']), prev[['Simpson']])
#   expect_equal(tolerance = 0.0001, unname(curr[,'Chao1']),   prev[['Chao1']])
# })
# 
# remove("curr", "prev")
# 
# 
# 
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Taxa Rollup
# #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# context("Taxa Rollup")
# 
# tr <- readRDS("outputs/taxa_rollup.rds")
# 
# # turn off locale-specific sorting
# orig_LCC <- Sys.getlocale("LC_COLLATE")
# on.exit(Sys.setlocale("LC_COLLATE", orig_LCC))
# Sys.setlocale("LC_COLLATE", "C")
# 
# for (txRank in taxa_ranks(biom)) {
#   test_that(txRank, {
#     expect_equal(tr[[txRank]], taxa_rollup(biom, txRank))
#   })
# }
# 
# remove("tr", "txRank")
# 
# remove("biom")
