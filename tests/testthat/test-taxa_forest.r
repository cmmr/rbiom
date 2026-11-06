
saliva_stool <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))
three_sites  <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool", "Buccal mucosa"))


test_that("taxa_forest", {

  p <- taxa_forest(saliva_stool, stat.by = "Body Site", rank = "Genus", p.top = 8)

  expect_s3_class(p, 'ggplot')
  expect_s3_class(p, 'rbiom_plot')
  expect_s3_class(p$stats, 'rbiom_tbl')
  expect_equal(nrow(p$data), 8L)

  # The most significant taxa, ordered by fold change.
  expect_setequal(p$data$.taxa, head(p$stats$.taxa[order(p$stats$.adj.p, p$stats$.p.val)], 8))
  expect_false(is.unsorted(p$data$.lfc))
  expect_identical(levels(p$data$.y), paste(1:8, p$data$.taxa, sep = "-"))

  # Two groups: the legend names them.
  expect_true(all(levels(p$data$.change) %in% c("Higher in Stool", "Higher in Saliva", "Not significant")))
  expect_identical(p$labels$x, "Log2 fold change")
  expect_identical(p$labels$y, "Genus")

  b <- ggplot2::ggplot_build(p)
  expect_identical(as.vector(b$layout$panel_scales_y[[1]]$get_labels()), as.character(p$data$.taxa))
})


test_that("taxa_forest panels, significance, and colors", {

  p <- taxa_forest(three_sites, stat.by = "Body Site", rank = "Genus", p.top = 5)
  expect_equal(nrow(p$data), 15L)
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 3L)
  expect_match(p$labels$caption, "\"A - B\"", fixed = TRUE)

  # p.top < 1 keeps the taxa below that adjusted p-value.
  p <- taxa_forest(saliva_stool, stat.by = "Body Site", rank = "Genus", p.top = 0.001)
  expect_true(all(p$data$.adj.p <= 0.001))
  expect_equal(nrow(p$data), sum(p$stats$.adj.p <= 0.001))

  # No taxa pass: a placeholder plot.
  p <- taxa_forest(saliva_stool, stat.by = "Sex", rank = "Genus", p.top = 1e-12)
  expect_s3_class(p, 'ggplot')
  expect_match(p$code, "No Data to Display", fixed = TRUE)

  # Significance threshold.
  p <- taxa_forest(saliva_stool, stat.by = "Body Site", rank = "Genus", p.top = Inf, p.sig = 1e-6)
  expect_identical(p$data$.change == "Not significant", p$data$.adj.p > 1e-6)

  # Colors.
  p <- taxa_forest(saliva_stool, stat.by = "Body Site", rank = "Genus", p.top = Inf, colors = c("red", "blue", "green"))
  cols <- unique(ggplot2::ggplot_build(p)$data[[3]][, c('colour', 'group')])
  expect_setequal(cols$colour, c("red", "blue", "green"))

  p <- taxa_forest(saliva_stool, stat.by = "Body Site", rank = "Genus", colors = FALSE)
  expect_null(p$labels$colour)

  expect_error(taxa_forest(saliva_stool, stat.by = "Body Site", colors = 1:4), "colors")
  expect_error(taxa_forest(saliva_stool, stat.by = "Body Site", bogus = 1), "bogus")
})


test_that("taxa_forest with regr and several ranks", {

  p <- taxa_forest(hmp50, regr = "Age", rank = c("Phylum", "Genus"), p.top = 4)

  expect_equal(nrow(p$data), 8L)
  expect_identical(p$labels$x, "Log2 fold change per unit of Age")
  expect_null(p$labels$y)
  expect_true(all(levels(p$data$.change) %in% c("Increases with Age", "Decreases with Age", "Not significant")))
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 2L)
})


test_that("taxa_forest with ANCOM-BC2", {

  skip_if_not_installed("ANCOMBC")
  skip_on_cran()

  p <- suppressWarnings(taxa_forest(saliva_stool, stat.by = "Body Site", rank = "Phylum", method = "ancombc2", p.top = Inf))

  sig <- p$data$.adj.p <= 0.05 & p$data$.robust %in% TRUE
  expect_identical(p$data$.change != "Not significant", sig)
  expect_match(p$labels$caption, "ANCOM-BC2")
  expect_match(p$code, ".robust %in% TRUE", fixed = TRUE)
})
