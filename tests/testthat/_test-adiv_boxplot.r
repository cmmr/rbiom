# 
# test_that("adiv_boxplot_1", {
#   p <- adiv_boxplot(rare50, x = "Body Site", adiv = "Shannon")
#   vdiffr::expect_doppelganger("adiv_boxplot_1_plot", p)
#   expect_snapshot(p[['stats']])
# })
# 
# test_that("adiv_boxplot_2", {
#   p <- adiv_boxplot(rare50, x = "Sex", adiv = c("OTUs", "Shannon"), layers="b", stat.by="Body Site", scales="free")
#   vdiffr::expect_doppelganger("adiv_boxplot_2_plot", p)
#   expect_snapshot(p[['stats']])
# })
# 
# test_that("adiv_boxplot_3", {
#   p <- adiv_boxplot(rare50, x = "Body Site", adiv = "Simpson", layers="p", stat.by="Sex", xlab.angle=30)
#   vdiffr::expect_doppelganger("adiv_boxplot_3_plot", p)
#   expect_snapshot(p[['stats']])
# })
