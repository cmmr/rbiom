
# Replays `p$code` in a clean environment and compares what ggplot2 draws.
expect_replays <- function (p) {
  code <- p$code
  expect_s3_class(code, 'rbiom_code')
  expect_false(any(grepl("not reproducible", code, fixed = TRUE)))
  env   <- new.env(parent = globalenv())
  r     <- eval(parse(text = code), envir = env)
  built <- function (x) lapply(ggplot2::ggplot_build(x)$data, function (d) { d <- as.data.frame(d); rownames(d) <- NULL; d })
  expect_equal(built(r), built(p), tolerance = 1e-6)
  invisible(r)
}

b20 <- rare50[1:20]


test_that("bdiv_ord_plot smoke", {

  expect_silent(bdiv_ord_plot(b20, permutations = 99))
  expect_silent(bdiv_ord_plot(b20, stat.by = 'Sex', permutations = 99))

  skip_on_cran()

  expect_silent(bdiv_ord_plot(b20, layers = 'sp', stat.by = 'Sex', facet.by = 'Body Site', test = 'none'))
  expect_silent(bdiv_ord_plot(b20, layers = 'pa', stat.by = 'Sex', facet.by = 'Sex', test = 'none', permutations = 99))
  expect_silent(bdiv_ord_plot(b20, layers = 'tm', taxa = c('Bacteroides', 'Prevotella')))
  expect_silent(bdiv_ord_plot(b20, layers = 'mn', bdiv = 'Manhattan'))
  expect_silent(bdiv_ord_plot(b20, layers = 'pt', bdiv = 'Jaccard'))

  expect_silent(bdiv_ord_plot(b20, rank = NULL))
  expect_silent(bdiv_ord_plot(b20, rank = c('Phylum', 'Genus')))
  expect_silent(bdiv_ord_plot(b20, bdiv = c('Manhattan', 'Jaccard')))
  expect_silent(bdiv_ord_plot(b20, bdiv = c('Manhattan', 'Jaccard'), facet.by = 'Sex'))
  expect_silent(bdiv_ord_plot(b20, bdiv = c('Manhattan', 'Jaccard'), rank = c('Phylum', 'Genus')))
  expect_silent(bdiv_ord_plot(b20, bdiv = c('Manhattan', 'Jaccard'), rank = c('Phylum', 'Genus'), facet.by = 'Sex'))
  expect_s3_class(suppressWarnings(bdiv_ord_plot(b20, ord = c('PCoA', 'NMDS'), layers = 'p')), 'rbiom_plot')
  expect_silent(bdiv_ord_plot(b20, stat.by = 'Sex', test = 'mrpp', permutations = 99))
  expect_silent(bdiv_ord_plot(b20, stat.by = 'Sex', test = 'none'))
  expect_silent(bdiv_ord_plot(b20, stat.by = 'Sex', strata = 'Body Site', permutations = 99))

  expect_error(bdiv_ord_plot(rare50[1:3]), "four samples")
  expect_error(bdiv_ord_plot(b20, shape.by = 'Sex'), "shape.by")
  expect_error(bdiv_ord_plot(b20, size = 3), "not an argument")
})


test_that("bdiv_ord_plot code replays", {

  skip_on_cran()

  expect_replays(bdiv_ord_plot(rare50, stat.by = "Body Site", permutations = 99))
  expect_replays(bdiv_ord_plot(rare50, stat.by = "Body Site", facet.by = "Sex", rank = "Genus", taxa = 3, permutations = 99))
  expect_replays(bdiv_ord_plot(b20, layers = "psa", bdiv = c("bray", "jaccard"), rank = c("Phylum", "Genus"), stat.by = "Sex", permutations = 99))
  expect_replays(bdiv_ord_plot(b20, layers = "pn", test = "none"))
})


test_that("bdiv_ord_plot labels and statistics", {

  skip_on_cran()

  p <- bdiv_ord_plot(rare50, stat.by = "Body Site", permutations = 99)
  expect_true(grepl('x +=  *"PCoA 1 [(]', p$code))
  expect_true(grepl('y +=  *"PCoA 2 [(]', p$code))
  expect_true(grepl('Bray-Curtis',    p$code, fixed = TRUE))
  expect_true(grepl('PERMANOVA F = ', p$code, fixed = TRUE))
  expect_true(grepl('type += "norm"', p$code))
  expect_s3_class(p$stats, 'data.frame')
  expect_s3_class(p$taxa_stats, 'data.frame')
  expect_s3_class(p$data$taxa_coords, 'data.frame')
  expect_identical(as.character(p$cmd), 'bdiv_ord_plot(biom, stat.by = "Body Site", permutations = 99)')
  expect_identical(as.character(p$stats$history), as.character(p$history))

  # Per-facet statistics go in the strips; generic axis titles.
  p <- bdiv_ord_plot(rare50, stat.by = "Body Site", facet.by = "Sex", rank = NULL, permutations = 99)
  expect_true(grepl("facets <- tibble::tibble", p$code, fixed = TRUE))
  expect_true(grepl('x += "PCoA 1"', p$code))
  expect_equal(nlevels(p$data$.facet), 2L)

  # Ellipses are skipped for groups with fewer than four samples.
  p <- bdiv_ord_plot(rare50[1:12], stat.by = "Body Site", test = "none")
  expect_true(grepl("ellipse_data <-", p$code, fixed = TRUE))
  expect_true(all(table(droplevels(attr(p$data, 'ellipse')[['Body Site']])) >= 4))
  expect_silent(ggplot2::ggplot_build(p))
  expect_replays(p)
})


test_that("bdiv_ord_plot density layer", {
  
  skip_on_cran()
  skip_if_not_installed("ggdensity")
  
  p <- expect_silent(bdiv_ord_plot(b20, layers = "pd", stat.by = "Sex", test = "none"))
  expect_true(grepl("ggdensity::geom_hdr", p$code, fixed = TRUE))
  expect_replays(p)
})


test_that("bdiv_ord_plot taxon labels without ggrepel", {

  skip_on_cran()

  old <- options(rbiom.disable_pkgs = "ggrepel")
  on.exit(options(old), add = TRUE)

  p <- bdiv_ord_plot(b20, layers = "ptm", test = "none")
  expect_false(grepl("ggrepel", p$code, fixed = TRUE))
  expect_true(grepl("geom_label(", p$code, fixed = TRUE))
  expect_replays(p)
})


test_that("bdiv_ord_plot with several metrics and ordinations", {

  skip_on_cran()
  skip_if_not_installed("uwot")

  # Used to fail with "argument is of length zero" (no common title).
  p <- bdiv_ord_plot(b20, stat.by = "Sex", bdiv = c("jaccard", "bray"), ord = c("PCoA", "UMAP"), layers = "pe", permutations = 99)
  expect_equal(nlevels(p$data$.facet), 4L)
  expect_true(all(grepl("PERMANOVA", levels(p$data$.facet), fixed = TRUE)))
  expect_replays(p)

  # With a biplot rank, the title is the rank.
  p <- bdiv_ord_plot(b20, bdiv = c("jaccard", "bray"), ord = c("PCoA", "UMAP"), layers = "pt", rank = "Phylum", taxa = 2, test = "none")
  expect_true(grepl('title = "Phylum"', p$code, fixed = TRUE))
})
