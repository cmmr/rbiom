test_that("plot_palettes", {
  
  expect_true(is_palette('-okabe'))
  expect_false(is_palette('dne'))
  
  skip_on_cran()
  
  
  expect_silent(get_palette('-okabe'))
  expect_silent(get_n_colors(12, TRUE))
  expect_silent(get_n_colors(20))
  expect_silent(get_n_colors(25))
  expect_silent(get_n_colors(5, c('red', 'blue')))
  expect_silent(get_n_shapes(12, TRUE))
  expect_silent(get_n_shapes(25))
  expect_silent(get_n_patterns(2, TRUE))
  
  expect_silent(color_palette(pal  = c(green = 'green')))
  expect_silent(color_palette(keys = c(green = 'green')))
  expect_silent(color_palette(keys = factor(LETTERS[1:3])))
  expect_silent(color_palette(keys = LETTERS[1:3]))
  expect_silent(color_palette(keys = FALSE))
  expect_silent(color_palette(n = 8))
  expect_silent(color_palette(n = 12))
  expect_silent(color_palette(n = 20))
  expect_silent(color_palette(n = 21))
  expect_silent(color_palette(pal =  1, n = 4))
  expect_silent(color_palette(pal = -1, n = 4))
  expect_error(color_palette(pal = FALSE, n = 4))
  
})
