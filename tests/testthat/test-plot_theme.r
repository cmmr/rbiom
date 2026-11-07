
test_that("theme_rbiom", {

  th <- theme_rbiom()
  expect_s3_class(th, 'theme')
  expect_true(attr(th, 'complete'))

  # Open axes, no gridlines or panel border, black text.
  el <- function (name, t = th) ggplot2::calc_element(name, t)
  expect_s3_class(el('axis.line.x.bottom'), 'element_line')
  expect_s3_class(el('panel.border'),       'element_blank')
  expect_s3_class(el('panel.grid.major.x'), 'element_blank')
  expect_s3_class(el('panel.grid.major.y'), 'element_blank')
  expect_s3_class(el('strip.background'),   'element_blank')
  expect_identical(el('axis.text.x.bottom')$colour, "black")
  expect_equal(el('axis.text.x.bottom')$size, 11 * 0.85)
  expect_equal(el('plot.title')$size, 11 * 1.15)
  expect_identical(th$plot.title.position, "plot")

  # Every text size follows `text`, so one setting resizes a plot.
  th8 <- th + ggplot2::theme(text = ggplot2::element_text(size = 8))
  expect_equal(el('axis.text.x.bottom', th8)$size, 8 * 0.85)
  expect_equal(el('legend.title', th8)$size, 8)

  # Options.
  expect_equal(el('axis.text.x.bottom', theme_rbiom(base_size = 8))$size, 8 * 0.85)
  expect_identical(el('text', theme_rbiom(base_family = "serif"))$family, "serif")

  box <- theme_rbiom(frame = "box")
  expect_s3_class(el('panel.border', box),       'element_rect')
  expect_s3_class(el('axis.line.x.bottom', box), 'element_blank')
  expect_s3_class(el('axis.ticks.x.bottom', box), 'element_line')

  void <- theme_rbiom(frame = "none")
  expect_s3_class(el('axis.line.x.bottom', void),  'element_blank')
  expect_s3_class(el('axis.ticks.x.bottom', void), 'element_blank')
  expect_s3_class(el('axis.title.x.bottom', void), 'element_blank')
  expect_s3_class(el('axis.text.x.bottom', void),  'element_text')

  expect_s3_class(el('panel.grid.major.y', theme_rbiom(grid = "y")),  'element_line')
  expect_s3_class(el('panel.grid.major.x', theme_rbiom(grid = "y")),  'element_blank')
  expect_s3_class(el('panel.grid.major.x', theme_rbiom(grid = TRUE)), 'element_line')
  expect_s3_class(el('panel.grid.minor',   theme_rbiom(grid = TRUE)), 'element_blank')

  expect_error(theme_rbiom(frame = "x"),       "frame")
  expect_error(theme_rbiom(grid = "z"),        "grid")
  expect_error(theme_rbiom(base_size = -1),    "base_size")
  expect_error(theme_rbiom(base_family = NA),  "base_family")

  # Usable on any ggplot.
  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point() + theme_rbiom()
  expect_silent(ggplot2::ggplot_gtable(ggplot2::ggplot_build(p)))
})


test_that("rbiom plots use theme_rbiom()", {

  skip_on_cran()

  p <- adiv_boxplot(rare50, x = "Sex", stat.by = "Sex")
  expect_true(grepl("theme_rbiom()", p$code, fixed = TRUE))
  expect_false(grepl("theme_bw", p$code, fixed = TRUE))
  expect_false(grepl("theme(", p$code, fixed = TRUE))

  # Theme options can be passed to any plot function.
  p <- adiv_boxplot(rare50, x = "Body Site", stat.by = "Body Site", base_size = 8, grid = "y")
  expect_true(grepl('theme_rbiom\\(\\s*base_size = 8,\\s*grid\\s*= "y" \\)', p$code))
  expect_true(grepl('fill\\s*= "white"', p$code))  # gridlines hidden behind p-values
  expect_equal(ggplot2::calc_element('axis.text.x.bottom', p$theme)$size, 8 * 0.85)

  p <- bdiv_ord_plot(rare50, stat.by = "Sex", layers = "p", permutations = 99)
  expect_true(grepl('theme_rbiom(frame = "box")', p$code, fixed = TRUE))

  p <- taxa_heatmap(rare50, taxa = 5)
  expect_true(grepl('theme_rbiom(frame = "none")', p$code, fixed = TRUE))
  expect_false(grepl("theme_void", p$code, fixed = TRUE))
})


test_that("every facet panel has axis lines", {

  skip_on_cran()

  has_axes <- function (p) grepl('axes\\s*= "all",\\s*axis.labels\\s*= "margins"', p$code)

  # One column: only the bottom panel would have an x-axis line.
  p <- adiv_corrplot(rare50, x = "Age", adiv = c("otu", "shan"), layers = "p", facet.ncol = 1)
  expect_true(has_axes(p))
  expect_identical(p$facet$params$draw_axes, list(x = TRUE, y = TRUE))
  expect_false(p$facet$params$axis_labels$x)  # tick labels on the outer panels only

  # facet_grid() always draws axes on the outer panels only.
  p <- adiv_boxplot(rare50, adiv = c("otu", "shan"), stat.by = "Body Site", facet.by = "Sex")
  expect_s3_class(p$facet, 'FacetGrid')
  expect_true(has_axes(p))

  # One row with free y scales: every panel has its axes already.
  p <- adiv_boxplot(rare50, x = "Sex", adiv = c("otu", "shan"), stat.by = "Sex")
  expect_false(grepl("axes", p$code, fixed = TRUE))

  # Free scales in facet_wrap(): same.
  p <- bdiv_ord_plot(rare50, stat.by = "Body Site", facet.by = "Sex", layers = "p", permutations = 99)
  expect_false(grepl("axes", p$code, fixed = TRUE))

  # The user's choice wins.
  p <- adiv_corrplot(rare50, x = "Age", adiv = c("otu", "shan"), layers = "p", facet.ncol = 1, facet.axes = "margins")
  expect_false(has_axes(p))
  expect_false(p$facet$params$draw_axes$x)
})


test_that("boxplots without x have no x-axis tick", {

  skip_on_cran()

  p <- adiv_boxplot(rare50, stat.by = "Body Site", flip = TRUE)
  expect_true(grepl("scale_x_discrete(breaks = NULL)", p$code, fixed = TRUE))
  expect_null(p$scales$get_scales("x")$breaks)

  p <- bdiv_boxplot(rare50, stat.by = "Sex", layers = "b")
  expect_null(p$scales$get_scales("x")$breaks)

  p <- adiv_boxplot(rare50, x = "Sex", stat.by = "Body Site")
  expect_false(grepl("breaks = NULL", p$code, fixed = TRUE))
})


test_that("default colors and shapes", {

  expect_identical(get_n_colors(3), c("#0072B2", "#D55E00", "#009E73"))
  expect_identical(get_n_colors(8), get_palette("okabe"))
  expect_identical(get_n_colors(12), get_palette("tol"))
  expect_identical(get_n_shapes(4), c(16, 17, 15, 18))

  skip_on_cran()

  # The "Other" taxa group is gray.
  p <- taxa_stacked(rare50, rank = "Phylum", taxa = 3)
  expect_identical(tail(levels(p$data$.taxa), 1), "Other")
  expect_true(grepl('values = c("#0072B2", "#D55E00", "#009E73", "#BBBBBB")', p$code, fixed = TRUE))
})


test_that("p-values in figures", {

  expect_identical(p_text(c(0.0738, 0.05, 0.00345, 2e-4, 0.001, 1)),
    c("p = 0.074", "p = 0.05", "p = 0.0034", "p < 0.001", "p = 0.001", "p = 1") )
  expect_identical(p_math(c(0.0738, 2e-4)),
    c('italic(p) == "0.074"', 'italic(p) < "0.001"') )

  # Every expression parses.
  for (x in p_math(c(0.5, 0.0123, 1e-10))) expect_silent(parse(text = x))

  skip_on_cran()

  p <- adiv_boxplot(rare50, x = "Body Site", stat.by = "Body Site", test = "wilcox")
  labels <- attr(p$data, 'stat_labels')$.label
  expect_true(length(labels) > 0)
  expect_false(any(grepl("e-0", labels, fixed = TRUE)))
  expect_true(all(grepl("^italic\\(p\\) (==|<) \"", labels)))
})


test_that("heatmap labels", {

  skip_on_cran()

  # label = "rows" hides the column labels (no numeric ticks).
  p <- taxa_heatmap(rare50, taxa = 5, tracks = "Sex", label = "rows")
  expect_null(p$scales$get_scales("x")$breaks)
  expect_false(is.null(p$scales$get_scales("y")$breaks))

  p <- taxa_heatmap(rare50, taxa = 5, tracks = "Sex", label = "cols")
  expect_null(p$scales$get_scales("y")$breaks)

  # base_size goes to theme_rbiom(); other arguments to theme().
  p <- taxa_heatmap(rare50, taxa = 5, base_size = 8, legend.position = "bottom")
  expect_true(grepl('theme_rbiom(8, frame = "none")', p$code, fixed = TRUE))
  expect_error(taxa_heatmap(rare50, taxa = 5, base_sz = 8), "theme_rbiom")
})
