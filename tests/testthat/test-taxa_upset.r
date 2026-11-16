
# Six samples of 100 reads each, so relative abundances can be read off
# the counts. `Site` has its levels in non-alphabetical order.
#
#       S1  S2  S3 | S4  S5  S6
#       --- Arm --- | --- Gut ---
#   G1  50  40  60 | 10   0   0
#   G2  30  40   0 |  0   0   0
#   G3  20  20  39 | 45  50  50
#   G4   0   0   1 | 45  50  30
#   G5   0   0   0 |  0   0  20
tiny_biom <- function () {
  as_rbiom(list(
    counts = matrix(
      data     = c(50, 40, 60, 10,  0,  0,
                   30, 40,  0,  0,  0,  0,
                   20, 20, 39, 45, 50, 50,
                    0,  0,  1, 45, 50, 30,
                    0,  0,  0,  0,  0, 20),
      nrow     = 5,
      byrow    = TRUE,
      dimnames = list(paste0("OTU", 1:5), paste0("S", 1:6)) ),
    metadata = data.frame(
      .sample = paste0("S", 1:6),
      Site    = factor(rep(c("Arm", "Gut"), each = 3), levels = c("Gut", "Arm")),
      Sex     = factor(c("F", "M", "F", "M", "F", "M")),
      Age     = c(21, 34, 45, 52, 28, 39) ),
    taxonomy = data.frame(
      .otu   = paste0("OTU", 1:5),
      Phylum = c("P1", "P1", "P2", "P2", "P2"),
      Genus  = paste0("G", 1:5) )))
}

# What ggplot2 draws for each layer of a plot.
built <- function (x) {
  lapply(ggplot2::ggplot_build(x)$data, function (d) { d <- as.data.frame(d); rownames(d) <- NULL; d })
}

# Replays `p$code` in a clean environment and compares what ggplot2 draws.
expect_replays <- function (p) {
  code <- p$code
  expect_s3_class(code, 'rbiom_code')
  expect_false(any(grepl("not reproducible", code, fixed = TRUE)))
  env <- new.env(parent = globalenv())
  r   <- eval(parse(text = code), envir = env)
  expect_true(inherits(r, 'patchwork'))
  expect_equal(length(r$patches$plots), length(p$patches$plots))
  for (i in seq_along(p$patches$plots))
    expect_equal(built(r$patches$plots[[i]]), built(p$patches$plots[[i]]), tolerance = 1e-6)
  expect_equal(built(r), built(p), tolerance = 1e-6)
  expect_equal(as.data.frame(r$data), as.data.frame(p$data), ignore_attr = TRUE)
  invisible(r)
}


test_that("taxa_upset: what is drawn", {

  biom <- tiny_biom()
  p    <- taxa_upset(biom, "Site")

  expect_s3_class(p, 'rbiom_plot')
  expect_s3_class(p, 'patchwork')
  expect_silent(ggplot2::ggplot_build(p))

  # Panels: set sizes, group totals, and (the outer plot) the dot matrix.
  expect_identical(length(p$patches$plots), 2L)
  bars   <- built(p$patches$plots[[1]])
  totals <- built(p$patches$plots[[2]])
  dots   <- built(p)

  # Bars: largest set first, with the count written above each.
  expect_equal(bars[[1]]$y, c(3, 1, 1))
  expect_equal(as.numeric(bars[[1]]$x), 1:3)
  expect_identical(as.character(bars[[2]]$label), c("3", "1", "1"))

  # Totals: taxa per group, the first level on the top row.
  expect_equal(totals[[1]]$x, c(-4, -4))   # reversed axis: bars grow leftwards

  # Matrix: one dot per set and group; filled dots joined by a line.
  expect_identical(nrow(dots[[1]]), 6L)
  expect_identical(nrow(dots[[2]]), 4L)
  expect_identical(nrow(dots[[3]]), 4L)

  # Gut (the first level) is the top row: y = 2 on a discrete axis.
  filled <- dots[[3]]
  expect_equal(sort(as.numeric(filled$y[filled$x == 2])), 2)   # set 2: Gut only
  expect_equal(sort(as.numeric(filled$y[filled$x == 3])), 1)   # set 3: Arm only
  expect_equal(sort(as.numeric(filled$y[filled$x == 1])), 1:2) # set 1: both
  expect_identical(
    ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels(),
    c("Arm", "Gut") )

  # Colour-neutral by default.
  expect_true(all(bars[[1]]$fill == "grey20"))
  expect_true(all(filled$colour == "black"))
})


test_that("taxa_upset: the table behind the plot", {

  biom <- tiny_biom()
  p    <- taxa_upset(biom, "Site", prev = 0.5)
  tbl  <- taxa_shared(biom, "Site", prev = 0.5)

  expect_equal(as.data.frame(p$stats), as.data.frame(tbl), ignore_attr = TRUE)
  expect_equal(as.data.frame(p$data),  as.data.frame(tbl), ignore_attr = TRUE)
  expect_equal(attr(p$stats, 'sets'), attr(tbl, 'sets'), ignore_attr = TRUE)

  expect_equal(built(p$patches$plots[[1]])[[1]]$y, attr(tbl, 'sets')$.n)

  # The table has its own code.
  expect_match(p$stats$code, 'taxa_shared(biom, "Site", prev = 0.5)', fixed = TRUE)
  expect_identical(as.character(p$cmd), 'taxa_upset(biom, "Site", prev = 0.5)')
})


test_that("taxa_upset: options", {

  skip_on_cran()

  # At most `sets` bars: the largest.
  full <- taxa_upset(rare50, "Body Site", sets = 99)
  p    <- taxa_upset(rare50, "Body Site", sets = 4)
  n    <- attr(p$stats, 'sets')$.n
  expect_true(length(n) > 4)
  expect_equal(built(p$patches$plots[[1]])[[1]]$y, n[1:4])
  expect_equal(built(full$patches$plots[[1]])[[1]]$y, n)
  expect_identical(nrow(built(p)[[1]]), 4L * 5L)

  # The totals count every taxon of a group, not just those in the bars shown.
  tbl <- p$stats
  expect_equal(
    -built(p$patches$plots[[2]])[[1]]$x,
    unname(vapply(levels(rare50$metadata[["Body Site"]]), function (g) sum(tbl[[g]]), numeric(1))) )

  # Without the totals panel.
  p <- taxa_upset(rare50, "Body Site", totals = FALSE)
  expect_identical(length(p$patches$plots), 1L)
  expect_false(grepl("totals", p$code, fixed = TRUE))

  # Arguments for the layers.
  p <- taxa_upset(rare50, "Sex", bar.fill = "steelblue", dot.size = 5, line.linewidth = 2, total.fill = "red", label.size = 4, empty.color = "white", base_size = 9)
  expect_true(all(built(p$patches$plots[[1]])[[1]]$fill == "steelblue"))
  expect_true(all(built(p$patches$plots[[1]])[[2]]$size == 4))
  expect_true(all(built(p$patches$plots[[2]])[[1]]$fill == "red"))
  expect_true(all(built(p)[[1]]$colour == "white"))
  expect_true(all(built(p)[[2]]$linewidth == 2))
  expect_true(all(built(p)[[3]]$size == 5))
  expect_match(p$code, 'geom_col(fill = "steelblue"', fixed = TRUE)
  expect_match(p$code, 'theme_rbiom(9', fixed = TRUE)

  # Thresholds, ranks, and lineage reach taxa_shared().
  p <- taxa_upset(rare50, "Body Site", rank = "Phylum", detect = 0.001, prev = 0.9, unc = "grouped", lineage = TRUE)
  expect_match(p$code, 'taxa_shared(biom, "Body Site", "Phylum", 0.001, 0.9, "grouped", TRUE)', fixed = TRUE)

  # Two to eight groups.
  b <- mutate(rare50, grp = factor(rep_len(LETTERS[1:8], length(Sex))))
  p <- taxa_upset(b, "grp")
  expect_silent(ggplot2::ggplot_build(p))
  expect_identical(
    ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels(),
    rev(LETTERS[1:8]) )
})


test_that("taxa_upset: errors", {

  biom <- tiny_biom()

  expect_error(taxa_upset(biom), "`stat.by` is required")
  expect_error(taxa_upset(biom, NULL))
  expect_error(taxa_upset(biom, "Age"), "not categorical")
  expect_error(taxa_upset(biom, "nosuchfield"), "nosuchfield")
  expect_error(taxa_upset(biom, "Site", sets = 0), "between 1 and Inf")
  expect_error(taxa_upset(biom, "Site", totals = "yes"), "TRUE or FALSE")
  expect_error(taxa_upset(biom, "Site", detect = 0.99), "No taxa belong to any group")
  expect_error(taxa_upset(biom, "Site", nosuch.size = 2), "not an argument")
  expect_error(taxa_upset(biom, "Site", zzz = 2), "not an argument")

  # A single group.
  one <- biom[1:3]
  expect_error(taxa_upset(one, "Site"), "at least two groups")
})


test_that("taxa_upset: schema", {

  x <- fn_schema("taxa_upset")
  expect_identical(x$fn, "taxa_upset")
  expect_identical(x$returns, list(type = "plot", has = c("data", "code", "methods", "stats")))
  expect_identical(x$counts, "rarefied")
  expect_identical(vapply(x$params, `[[`, "", 'name'), names(formals(taxa_upset)))
  expect_true(x$params[[2]]$required)

  d <- fn_check_data(hmp50)
  expect_identical(
    fn_check("taxa_upset", list(stat.by = "body", rank = -2, sets = 10), d),
    list(stat.by = "Body Site", rank = "Family", sets = 10) )
  expect_error(fn_check("taxa_upset", list(), d), "needs")
  expect_error(fn_check("taxa_upset", list(stat.by = "Age"), d), "not categorical")
  expect_error(fn_check("taxa_upset", list(stat.by = "Sex", sets = 0), d), "between 1 and Inf")
})


test_that("taxa_upset code replays", {

  skip_on_cran()

  fresh <- as_rbiom(system.file("extdata", "hmp50.bz2", package = "rbiom"))
  b     <- rarefy(fresh, depth = 1000, seed = 3)

  expect_replays(taxa_upset(b, "Body Site"))
  expect_replays(taxa_upset(b, "Body Site", rank = "Phylum", detect = 0.001, prev = 0.5, sets = 5))
  expect_replays(taxa_upset(b, "Sex", totals = FALSE, bar.fill = "steelblue", base_size = 9))
  expect_replays(taxa_upset(b, "Sex", unc = "drop", lineage = TRUE, dot.size = 4, total.width = 0.4))

  # Cached and uncached plots agree.
  p <- taxa_upset(b, "Body Site")
  rlang::with_options(rbiom.cache_dir = "FALSE", {
    q <- taxa_upset(b, "Body Site")
  })
  expect_equal(built(q), built(p))
  expect_identical(as.character(q$code), as.character(p$code))
})
