
test_that("adiv_boxplot_1", {
  p <- adiv_boxplot(rare50, x = "Body Site", metric = "Shannon")
  vdiffr::expect_doppelganger("adiv_boxplot_1_plot", p)
  expect_snapshot(p[['stats']])
})

test_that("adiv_boxplot_2", {
  p <- adiv_boxplot(rare50, x = "Sex", metric = c("OTUs", "Shannon"), layers="b", color.by="Body Site", scales="free")
  vdiffr::expect_doppelganger("adiv_boxplot_2_plot", p)
  expect_snapshot(p[['stats']])
})

test_that("adiv_boxplot_3", {
  p <- adiv_boxplot(rare50, x = "Body Site", metric = "Simpson", layers="p", color.by="Sex", xlab.angle=30)
  vdiffr::expect_doppelganger("adiv_boxplot_3_plot", p)
  expect_snapshot(p[['stats']])
})
