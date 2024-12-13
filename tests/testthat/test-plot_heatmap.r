test_that("plot_heatmap", {
  
  mtx           <- matrix(1:25, ncol = 5)
  dimnames(mtx) <- list(letters[1:5], LETTERS[1:5])
  cat_vals      <- c("Male", "Female")[c(1,2,1,2,1)]
  num_vals      <- c(25, 30, 30, 35, 40)
  tracks        <- list(
    list(label = "Sex", values = cat_vals, side = "left", colors   = "bright"),
    list(label = "Age", values = num_vals, side = "left", range    = c(20,40), bins = 3),
    list(label = "Sex", values = cat_vals, side = "top",  na.color = "gray"),
    list(label = "Age", values = num_vals, side = "top",  colors   = c('red', 'blue')) )
  
  expect_silent(plot_heatmap(mtx = mtx, tracks = tracks))
  
  skip_on_cran()
  
  expect_silent(plot_heatmap(mtx = mtx, grid = NULL, rescale = 'rows', clust = FALSE))
  expect_silent(plot_heatmap(mtx = mtx, grid = 'greens', rescale = 'cols'))
  expect_silent(plot_heatmap(mtx = mtx, grid = c('green', 'blue')))
  expect_silent(plot_heatmap(mtx = mtx, tracks = list(list(label = "Sex", values = cat_vals, side = "left"))))
  
})
