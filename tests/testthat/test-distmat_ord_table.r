
replay_stats <- function (x) {
  expect_false(any(grepl("not reproducible", x$code, fixed = TRUE)))
  env <- new.env(parent = globalenv())
  eval(parse(text = x$code), envir = env)
  env$stats
}

strip <- function (df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (a in c('cmd', 'code', 'history', 'data_cmd', 'tbl_sum', 'stats', 'taxa_stats', 'taxa_coords', 'axes')) attr(df, a) <- NULL
  rownames(df) <- NULL
  df
}


test_that("distmat_ord_table", {

  dm15 <- expect_silent(bdiv_distmat(rare50[1:15]))
  ord  <- expect_silent(distmat_ord_table(dm15, ord = "PCoA"))

  # PCoA is classical multidimensional scaling.
  pc <- stats::cmdscale(dm15, k = 2, eig = TRUE)
  expect_equal(ord$.x, unname(pc$points[, 1]))
  expect_equal(ord$.y, unname(pc$points[, 2]))
  expect_identical(as.character(ord$.sample), attr(dm15, 'Labels'))

  axes <- attr(ord, 'axes')
  expect_equal(axes$.var, 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]))
  expect_match(axes$.label[[1]], "^PCoA 1 \\([0-9.]+%\\)$")

  expect_error(distmat_ord_table(bdiv_distmat(rare50[1:2])), "three samples")

  skip_on_cran()

  dm40 <- expect_silent(bdiv_distmat(rare50[1:40]))
  expect_silent(distmat_ord_table(dm15, ord = "PCoA", k = 3))
  expect_silent(distmat_ord_table(dm15, ord = "UMAP"))
  expect_silent(distmat_ord_table(dm40, ord = "NMDS"))

  # Stochastic ordinations are reproducible through `seed`.
  n1 <- distmat_ord_table(dm40, ord = "NMDS", seed = 1)
  n2 <- distmat_ord_table(dm40, ord = "NMDS", seed = 1)
  expect_equal(n1$.x, n2$.x)
  set.seed(99); before <- .Random.seed
  distmat_ord_table(dm40, ord = "NMDS", seed = 2)
  expect_identical(.Random.seed, before)
})


test_that("distmat_ord_table tSNE", {

  skip_on_cran()
  skip_if_not_installed("Rtsne")

  dm15 <- bdiv_distmat(rare50[1:15])

  # Rtsne under `seed`, with perplexity capped at (n - 1) / 3.
  ord <- expect_silent(distmat_ord_table(dm15, ord = "tSNE", seed = 3))
  ref <- with_seed(3, Rtsne::Rtsne(dm15, dims = 2, perplexity = 14 / 3)$Y)
  expect_equal(ord$.x, ref[, 1])
  expect_equal(ord$.y, ref[, 2])
  expect_identical(as.character(ord$.sample), attr(dm15, 'Labels'))
  expect_identical(attr(ord, 'tbl_sum')[['Ordination']], "Rtsne::Rtsne(dm, dims = 2)")

  # `k` becomes `dims`; `...` reaches Rtsne.
  ord <- distmat_ord_table(dm15, ord = "tSNE", k = 3, seed = 3, perplexity = 2, theta = 0)
  ref <- with_seed(3, Rtsne::Rtsne(dm15, dims = 3, perplexity = 2, theta = 0)$Y)
  expect_equal(ord$.z, ref[, 3])

  # Works down to the three-sample minimum.
  expect_silent(distmat_ord_table(bdiv_distmat(rare50[1:3]), ord = "tSNE"))
})


test_that("bdiv_ord_table", {

  skip_on_cran()

  b   <- rarefy(hmp50, depth = 1000, seed = 7)
  ord <- bdiv_ord_table(b, bdiv = "bray", ord = "PCoA", stat.by = "Body Site", rank = "Genus", taxa = 3, permutations = 199)

  # Distances are not relativized first.
  dm <- bdiv_distmat(b, "bray")
  pc <- stats::cmdscale(dm, k = 2, eig = TRUE)
  expect_equal(ord$.x, unname(pc$points[, 1]))
  expect_equal(ord$stats$.stat, distmat_stats(dm, pull(b, "Body Site"), permutations = 199)$.stat)

  # Sample statistics replay.
  expect_equal(strip(replay_stats(ord$stats)), strip(ord$stats), tolerance = 1e-6)

  # Taxa biplot: envfit r2 / p-values, adjusted over all taxa, with arrows.
  ts <- ord$taxa_stats
  expect_identical(names(ts), c('.bdiv', '.rank', '.ord', '.taxa', '.r.sqr', '.p.val', '.adj.p'))
  expect_equal(nrow(ts), 3L)
  expect_equal(ts$.adj.p, p.adjust(ts$.p.val, 'fdr'))
  expect_equal(strip(replay_stats(ts)), strip(ts), tolerance = 1e-6)

  tc <- ord$taxa_coords
  expect_true(all(c('.taxa', '.size', '.x', '.y', '.x0', '.y0') %in% names(tc)))
  expect_equal(tc$.x0, rep(mean(ord$.x), 3), tolerance = 1e-5)

  # p.top keeps the most significant taxa; .adj.p was computed before filtering.
  o2 <- bdiv_ord_table(b, bdiv = "bray", rank = "Genus", taxa = 6, p.top = 2, permutations = 99)
  expect_equal(nrow(o2$taxa_stats), 2L)
  expect_equal(nrow(o2$taxa_coords), 2L)
  o6 <- bdiv_ord_table(b, bdiv = "bray", rank = "Genus", taxa = 6, permutations = 99)
  expect_equal(o2$taxa_stats$.adj.p, o6$taxa_stats$.adj.p[match(o2$taxa_stats$.taxa, o6$taxa_stats$.taxa)])

  # Several metrics / ordinations / strata; blocked permutations.
  o3 <- bdiv_ord_table(b, bdiv = c("bray", "jaccard"), ord = c("PCoA", "NMDS"), split.by = "Sex", stat.by = "Body Site", permutations = 49)
  expect_equal(nrow(o3), 4 * b$n_samples)
  expect_equal(nrow(o3$stats), 4L)
  expect_equal(strip(replay_stats(o3$stats)), strip(o3$stats), tolerance = 1e-6)

  expect_error(bdiv_ord_table(b[1:2]), "three samples")

  # Repeated measures: a between-subject stat.by permutes whole subjects.
  bb <- rarefy(subset(babies, `Age (days)` <= 10))
  o4 <- bdiv_ord_table(bb, stat.by = "Delivery mode", strata = "Subject ID")
  expect_match(attr(o4$stats, 'tbl_sum')[['Test']], "whole `Subject ID` strata permuted .*all 220 arrangements, exact")
  expect_lt(o4$stats$.p.val, 1)
  expect_equal(strip(replay_stats(o4$stats)), strip(o4$stats), tolerance = 1e-6)
})
