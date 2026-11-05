
# Replays `p$code` in a clean environment and compares what ggplot2 draws.
expect_replays <- function (p) {
  code <- p$code
  expect_s3_class(code, 'rbiom_code')
  expect_false(any(grepl("not reproducible", code, fixed = TRUE)))
  env   <- new.env(parent = globalenv())
  r     <- eval(parse(text = code), envir = env)
  built <- function (x) lapply(ggplot2::ggplot_build(x)$data, function (d) { d <- as.data.frame(d); rownames(d) <- NULL; d })
  if (inherits(p, 'patchwork')) {
    expect_true(inherits(r, 'patchwork'))
    expect_equal(length(r$patches$plots), length(p$patches$plots))
    for (i in seq_along(p$patches$plots))
      expect_equal(built(r$patches$plots[[i]]), built(p$patches$plots[[i]]), tolerance = 1e-6)
  }
  expect_equal(built(r), built(p), tolerance = 1e-6)
  invisible(r)
}

# A fresh (unmodified) object: other test files edit the shared fixtures in
# place, which adds non-reproducible entries to their history.
fresh5 <- as_rbiom(system.file("extdata", "hmp50.bz2", package = "rbiom"))[1:5]


test_that("taxa_stacked", {

  expect_silent(taxa_stacked(rare5))

  skip_on_cran()

  expect_silent(taxa_stacked(rare5, rank = c("Genus", "Phylum")))
  expect_silent(taxa_stacked(rare5, label.by = 'Body Site'))
  expect_silent(taxa_stacked(rare5, label.by = '.sample', xlab.angle = 30))
  expect_silent(taxa_stacked(rare5, label.by = 'Body Site', order.by = 'Sex'))
  expect_silent(taxa_stacked(rare5, order.by = c('Sex', 'Age'), facet.by = 'Sex'))
  expect_silent(taxa_stacked(rare5, colors = FALSE))
  expect_error(taxa_stacked(rare5, size = 2), "not an argument")

  skip_if_not_installed("fillpattern")
  expect_silent(taxa_stacked(hmp5, patterns = TRUE, order.by = 'Sex'))
})


test_that("taxa_stacked abundances and labels", {

  skip_on_cran()

  # Rarefied: relative abundance; raw counts otherwise.
  p <- taxa_stacked(rare5, rank = "Phylum")
  expect_true(grepl("Relative Abundance", p$code, fixed = TRUE))
  expect_true(all(p$data$.abundance <= 1))
  expect_identical(as.character(p$cmd), 'taxa_stacked(biom, "Phylum")')
  expect_true(grepl("mtx <- taxa_matrix(biom", p$code, fixed = TRUE))

  p <- taxa_stacked(hmp5, rank = "Phylum")
  expect_true(grepl("Raw Abundance", p$code, fixed = TRUE))
  expect_true(any(p$data$.abundance > 1))

  # Labels other than sample names are kept unique with a prefix.
  p <- taxa_stacked(rare5, rank = "Phylum", label.by = "Body Site")
  expect_true(all(grepl("^[0-9]+-", levels(p$data$.sample))))
  expect_true(grepl("Labeled by Body Site", p$code, fixed = TRUE))

  # Samples with NA metadata are dropped and recorded.
  b <- rare50[1:8]
  md <- b$metadata; md$Sex[[1]] <- NA; b$metadata <- md
  p <- taxa_stacked(b, rank = "Phylum", facet.by = "Sex")
  expect_equal(nlevels(p$data$.sample), 7L)
  expect_true(grepl("na.omit(biom", p$history, fixed = TRUE))
})


test_that("taxa_stacked code replays", {

  skip_on_cran()

  b <- rare50[1:12]
  expect_replays(taxa_stacked(b, rank = "Phylum", taxa = 5))
  expect_replays(taxa_stacked(b, rank = c("Phylum", "Genus"), taxa = 4))
  expect_replays(taxa_stacked(b, rank = "Phylum", label.by = "Body Site", order.by = "Sex", facet.by = "Sex"))
  expect_replays(taxa_stacked(fresh5, rank = "Phylum", label.by = ".sample"))
})
