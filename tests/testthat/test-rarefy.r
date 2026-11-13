
replay_biom <- function (b) {
  env <- new.env(parent = globalenv())
  eval(parse(text = b$code), envir = env)
  env$biom
}


test_that("rarefy", {

  skip_on_cran()

  mtx <- min5$counts

  expect_error(suppressWarnings(mtx_rarefy(mtx, margin = 2L, depth = 1:2)))
  expect_error(suppressWarnings(mtx_rarefy(mtx, margin = 2L, depth = 1.5)))
  expect_error(suppressWarnings(mtx_rarefy(mtx, margin = 2L, seed = 1.5)))

  expect_warning(mtx_rarefy(mtx, margin = 1L))
  expect_warning(mtx_rarefy(mtx, margin = 2L))

  expect_silent(suggest_rarefy_depth(mtx))

  # The deprecated `upsample` path inflates before rarefying.
  rel <- biom_relativize(min5)$counts
  expect_equal(unname(colSums(suppressWarnings(mtx_rarefy(rel, depth = 500, upsample = 1000)))), rep(500, 5))
})


test_that("rarefy is reproducible and records its depth and seed", {

  h5 <- hmp50[1:5]
  r1 <- rarefy(h5, depth = 1000, seed = 7)
  r2 <- rarefy(h5, depth = 1000, seed = 7)
  r3 <- rarefy(h5, depth = 1000, seed = 8)
  expect_identical(r1$hash, r2$hash)
  expect_false(identical(r1$hash, r3$hash))
  expect_equal(unname(sample_sums(r1)), rep(1000, 5))
  expect_equal(r1$depth, 1000)

  expect_true(grepl("rarefy(biom, depth = 1000L, seed = 7)", r1$code, fixed = TRUE))
  expect_equal_rbiom(replay_biom(r1), r1)

  # Auto-selected depths are written into the history.
  auto <- rarefy(hmp50)
  expect_true(grepl(sprintf("depth = %dL", suggest_rarefy_depth(hmp50)), auto$code, fixed = TRUE))
  expect_equal_rbiom(replay_biom(auto), auto)

  expect_error(rarefy(biom_relativize(hmp5)), "inflate")
})


test_that("biom_relativize / biom_rescale / biom_inflate", {

  h5  <- hmp50[1:5]
  rel <- biom_relativize(h5)
  expect_true(rel$relative)
  expect_null(rel$depth)
  expect_equal(unname(sample_sums(rel)), rep(1, 5))
  expect_equal_rbiom(replay_biom(rel), rel)

  rs <- biom_rescale(h5, range = c(0, 100))
  expect_equal(unname(apply(as.matrix(rs), 2, max)), rep(100, 5))
  expect_true(grepl("biom_rescale(biom, range = c(0, 100))", rs$code, fixed = TRUE))
  expect_equal_rbiom(replay_biom(rs), rs)

  # Integer counts are left alone when no depth is requested.
  inf <- biom_inflate(h5)
  expect_equal(as.matrix(inf), as.matrix(h5))

  inf <- biom_inflate(rel, depth = 1000)
  expect_equal(unname(sample_sums(inf)), rep(1000, 5))
  expect_equal_rbiom(replay_biom(inf), inf)

  inf <- biom_inflate(rel)
  expect_false(any(inf$counts@x %% 1 != 0))
  expect_true(all(sample_sums(inf) > 1))

  inf <- rarefy(rel, depth = 500, inflate = TRUE)
  expect_equal(unname(sample_sums(inf)), rep(500, 5))
  expect_equal_rbiom(replay_biom(inf), inf)
})


test_that("biom_inflate recovers counts divided by a constant, and proportions", {

  # Whole multiples of each sample's smallest value give exact depths, on
  # any scale (the atima report: counts / 4 inflated to a median depth of 3).
  m     <- as.matrix(hmp50$counts)
  depth <- stats::setNames(as.integer(sample_sums(hmp50)), hmp50$samples)

  for (x in list(m / 4, m * 2.5, sweep(m, 2, colSums(m), "/"))) {
    b <- as_rbiom(list(counts = x))
    expect_identical(suggest_inflate_depths(b), depth)
    expect_equal(as.matrix(biom_inflate(b)$counts), m)
  }
  expect_equal(as.matrix(biom_inflate(biom_relativize(hmp50))$counts), m)

  # Otherwise the heuristic estimates depths from proportions, so the scale
  # of the values doesn't matter.
  set.seed(1)
  noisy <- m * (1 + 1e-4 * matrix(stats::runif(length(m)), nrow(m)))
  est   <- suggest_inflate_depths(as_rbiom(list(counts = noisy)))
  expect_identical(suggest_inflate_depths(as_rbiom(list(counts = noisy / 4))), est)
  expect_true(all(est > 10))
  expect_false(identical(est, depth))
})


test_that("suggest_inflate_depths handles uninformative samples", {

  # Integer counts: the current sums.
  expect_equal(unname(suggest_inflate_depths(hmp5)), unname(sample_sums(hmp5)))

  # A sample with a single feature gets a depth too.
  counts <- as.matrix(biom_relativize(hmp5)$counts)
  counts <- cbind(counts, Solo = c(1, rep(0, nrow(counts) - 1)))
  depths <- suggest_inflate_depths(as_rbiom(counts))
  expect_false(anyNA(depths))
  expect_true(all(depths >= 1))

  expect_equal(unname(sample_sums(rarefy(as_rbiom(counts), depth = 100, inflate = TRUE))), rep(100, 6))
})
