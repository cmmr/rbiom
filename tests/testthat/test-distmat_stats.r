
replay_stats <- function (x) {
  expect_false(any(grepl("not reproducible", x$code, fixed = TRUE)))
  env <- new.env(parent = globalenv())
  eval(parse(text = x$code), envir = env)
  env$stats
}

strip <- function (df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (a in c('cmd', 'code', 'history', 'data_cmd', 'tbl_sum')) attr(df, a) <- NULL
  rownames(df) <- NULL
  df
}


test_that("distmat_stats", {

  dm <- bdiv_distmat(rare50)
  expect_null(distmat_stats(dm, pull(rare50, 'Sex'), test = 'none'))

  x <- distmat_stats(dm, pull(rare50, 'Sex'), permutations = 199)
  x <- distmat_stats(dm, pull(rare50, 'Sex'), permutations = 199) # cached
  expect_identical(names(x), c('.n', '.stat', '.r.sqr', '.p.val'))
  expect_equal(x$.n, 49L)

  # Matches vegan::adonis2 with the same seed and permutations.
  grp <- pull(rare50, 'Sex')[attr(dm, 'Labels')]
  set.seed(0); a <- vegan::adonis2(dm ~ grp, permutations = 199)
  expect_equal(x$.stat,  a$F[[1]])
  expect_equal(x$.r.sqr, a$R2[[1]])
  expect_equal(x$.p.val, a$`Pr(>F)`[[1]])

  # The embedded code reproduces the result.
  expect_equal(strip(replay_stats(x)), strip(x))
  expect_true(grepl("vegan::adonis2(dm ~ grouping, permutations = 199)", x$code, fixed = TRUE))

  # Unnamed groups in matrix order are accepted.
  y <- distmat_stats(dm, unname(grp), permutations = 199)
  expect_equal(strip(x), strip(y))
})


test_that("mrpp, betadisper, and mantel", {

  skip_on_cran()

  dm  <- bdiv_distmat(rare50)
  grp <- pull(rare50, 'Body Site')

  x <- distmat_stats(dm, grp, test = "mrpp", permutations = 99)
  expect_identical(names(x), c('.n', '.stat', '.z', '.p.val'))
  set.seed(0); m <- vegan::mrpp(dm, grp[attr(dm, 'Labels')], permutations = 99)
  expect_equal(x$.stat, m$delta)
  expect_equal(x$.p.val, m$Pvalue)
  expect_equal(strip(replay_stats(x)), strip(x))

  x <- suppressWarnings(distmat_stats(dm, grp, test = "betadisper", permutations = 99))
  expect_identical(names(x), c('.n', '.stat', '.p.val'))
  expect_equal(strip(suppressWarnings(replay_stats(x))), strip(x))

  age <- pull(rare50, 'Age')
  x <- distmat_stats(dm, age, test = "mantel", permutations = 99)
  set.seed(0); m <- vegan::mantel(dm, dist(age[attr(dm, 'Labels')]), method = 'spearman', permutations = 99)
  expect_equal(x$.stat, m$statistic)
  expect_equal(x$.p.val, m$signif)
  expect_equal(strip(replay_stats(x)), strip(x))

  # Numeric covariates with PERMANOVA; wrong types are rejected.
  expect_s3_class(distmat_stats(dm, age, permutations = 49), 'rbiom_tbl')
  expect_error(distmat_stats(dm, age, test = "mrpp"), "numeric")
  expect_error(distmat_stats(dm, grp, test = "mantel"), "numeric")
  expect_error(distmat_stats(dm, grp[1:10]), "missing value")
  expect_error(distmat_stats(dm, unname(grp)[1:10]), "one value per sample")
})


#________________________________________________________
# Repeated-measures data: 12 infants (3 Cesarean, 9 vaginal), 10 daily
# samples each.
#________________________________________________________
babies_subset <- function () {
  rarefy(subset(babies, `Age (days)` <= 10))
}


test_that("a within-strata variable permutes within blocks", {

  skip_on_cran()

  b   <- babies_subset()
  dm  <- bdiv_distmat(b)
  lab <- attr(dm, 'Labels')
  sid <- pull(b, 'Subject ID')[lab]
  age <- pull(b, 'Age (days)')[lab]

  early <- ifelse(age <= stats::median(age), "early", "late")
  expect_identical(strata_scheme(early, sid), "within")

  free    <- distmat_stats(dm, early, permutations = 199)
  blocked <- distmat_stats(dm, early, strata = sid, permutations = 199)

  expect_equal(free$.stat, blocked$.stat)
  expect_identical(perm_plan(early, sid, 199)$scheme, "within")
  expect_match(attr(blocked, 'tbl_sum')[['Test']], "199 permutations within strata")
  expect_true(grepl("permute::how(nperm = 199, blocks = blocks)", blocked$code, fixed = TRUE))
  expect_equal(strip(replay_stats(blocked)), strip(blocked))

  set.seed(0)
  v <- vegan::adonis2(dm ~ early, permutations = permute::how(nperm = 199, blocks = factor(sid)))
  expect_equal(blocked$.p.val, v$`Pr(>F)`[[1]])
})


test_that("a between-strata variable permutes whole strata", {

  skip_on_cran()

  b   <- babies_subset()
  dm  <- bdiv_distmat(b)
  lab <- attr(dm, 'Labels')
  sid <- as.character(pull(b, 'Subject ID')[lab])
  del <- as.character(pull(b, 'Delivery mode')[lab])
  expect_identical(strata_scheme(del, sid), "between")

  # Observed statistics are identical to vegan's.
  x <- distmat_stats(dm, del, strata = sid)
  v <- vegan::adonis2(dm ~ del, permutations = 9)
  expect_equal(x$.stat,  v$F[[1]],  tolerance = 1e-10)
  expect_equal(x$.r.sqr, v$R2[[1]], tolerance = 1e-10)

  subj <- stats::setNames(seq_along(unique(sid)) * 1.7, unique(sid))
  cov  <- stats::setNames(unname(subj[sid]), lab)
  x2   <- distmat_stats(dm, cov, strata = sid, permutations = 199)
  v2   <- vegan::adonis2(dm ~ cov, permutations = 9)
  expect_equal(x2$.stat,  v2$F[[1]],  tolerance = 1e-10)
  expect_equal(x2$.r.sqr, v2$R2[[1]], tolerance = 1e-10)

  x3 <- distmat_stats(dm, cov, test = "mantel", strata = sid, permutations = 199)
  m3 <- vegan::mantel(dm, stats::dist(cov), method = "spearman", permutations = 9)
  expect_equal(x3$.stat, unname(m3$statistic), tolerance = 1e-10)

  # 3 Cesarean vs 9 vaginal infants: all choose(12, 3) = 220 arrangements
  # are enumerated, so the p-value is exact. Check it against a brute-force
  # computation on the full design, and that it isn't the degenerate 1
  # that within-subject blocks give.
  plan <- perm_plan(del, sid, 999)
  expect_identical(plan$scheme, "between")
  expect_true(plan$exact)
  expect_identical(plan$n, 220L)
  expect_match(attr(x, 'tbl_sum')[['Test']], "whole strata permuted .*all 220 arrangements, exact")

  # Pseudo-F straight from the squared distances (Anderson 2001):
  # SS_T = sum of d^2 over pairs / n, SS_W = the same within each group /
  # its size. No matrix products, so the check is independent of
  # distmat_stats()'s algebra and of the platform's BLAS (OpenBLAS 0.3.20
  # on an ubuntu-22.04 runner once got the earlier, product-based version
  # wrong).
  D2   <- as.matrix(dm)^2; n <- nrow(D2)
  ss_t <- sum(D2) / (2 * n)
  pseudo_F <- function (g) {
    ss_w <- sum(vapply(split(seq_len(n), g), FUN.VALUE = 0, function (i) sum(D2[i, i]) / (2 * length(i))))
    (ss_t - ss_w) / (ss_w / (n - 2))
  }
  expect_equal(pseudo_F(del), x$.stat, tolerance = 1e-10)
  subjects <- sort(unique(sid))
  n_ces    <- length(unique(sid[del == "Cesarean"]))
  splits   <- utils::combn(subjects, n_ces, simplify = FALSE)
  all_F    <- vapply(splits, FUN.VALUE = 0, function (s) pseudo_F(sid %in% s))
  expect_equal(x$.p.val, mean(all_F >= x$.stat - 1e-8))
  expect_lt(x$.p.val, 1)

  # A between-subject factor with a strong effect is detected: the most
  # extreme split of the infants has the smallest attainable p, 1/220.
  extreme <- ifelse(sid %in% splits[[which.max(all_F)]], "A", "B")
  xe      <- distmat_stats(dm, extreme, strata = sid)
  expect_equal(xe$.p.val, 1 / 220)
  expect_lt(xe$.p.val, 0.05)

  # $code calls distmat_stats() itself, and replays exactly.
  expect_true(grepl("distmat_stats(dm, groups, test = \"adonis2\", seed = 0, permutations = 999, strata = strata)", x$code, fixed = TRUE))
  expect_false(grepl("permute::how", x$code, fixed = TRUE))
  expect_equal(strip(replay_stats(x)), strip(x))
  expect_equal(strip(replay_stats(x3)), strip(x3))

  # Random (not exhaustive) permutations when arrangements outnumber them.
  x4 <- distmat_stats(dm, cov, strata = sid, permutations = 99, seed = 1)
  expect_false(perm_plan(cov, sid, 99)$exact)
  expect_match(attr(x4, 'tbl_sum')[['Test']], "between-strata variable; 99 permutations")
  expect_identical(x4, distmat_stats(dm, cov, strata = sid, permutations = 99, seed = 1))
  expect_true(x4$.p.val >= 1 / 100)

  # MRPP and betadisper can't permute whole strata.
  expect_error(distmat_stats(dm, del, strata = sid, test = "mrpp"),       "constant within each stratum")
  expect_error(distmat_stats(dm, del, strata = sid, test = "betadisper"), "constant within each stratum")

  # Coarse null distributions are flagged.
  few <- sid %in% subjects[1:4]
  dm4 <- stats::as.dist(as.matrix(dm)[few, few])
  expect_warning(
    distmat_stats(dm4, ifelse(sid[few] %in% subjects[1:2], "A", "B"), strata = sid[few], permutations = 99),
    "distinct arrangement" )
})


test_that("whole-strata permutation has nominal type-I error", {

  skip_on_cran()

  b   <- babies_subset()
  dm  <- bdiv_distmat(b)
  lab <- attr(dm, 'Labels')
  sid <- as.character(pull(b, 'Subject ID')[lab])
  us  <- sort(unique(sid))

  rlang::with_options(rbiom.cache_dir = "FALSE", {
    set.seed(2)
    p <- t(vapply(1:100, FUN.VALUE = c(0, 0), function (i) {
      f <- stats::setNames(sample(rep(c("A", "B"), length.out = length(us))), us)
      g <- unname(f[sid])
      c(distmat_stats(dm, g, permutations = 199, seed = i)$.p.val,
        distmat_stats(dm, g, strata = sid, permutations = 199, seed = i)$.p.val )
    }))
  })

  # Ignoring the repeated measures is badly anti-conservative ...
  expect_gt(mean(p[, 1] < 0.05), 0.3)
  # ... while permuting whole subjects holds the nominal level.
  expect_lte(mean(p[, 2] < 0.05), 0.12)
})


test_that("NA groups are dropped", {

  skip_on_cran()

  # NA groups are dropped with a message.
  grp <- pull(rare50, 'Sex'); grp[1:2] <- NA
  dm  <- bdiv_distmat(rare50)
  expect_message(x <- distmat_stats(dm, grp, permutations = 49), "Dropped 2")
  expect_equal(x$.n, 47L)
  expect_equal(strip(replay_stats(x)), strip(x))
})


test_that("distmat_stats leaves the RNG alone", {
  set.seed(5)
  before <- .Random.seed
  distmat_stats(bdiv_distmat(rare5), pull(rare5, 'Sex'), permutations = 19)
  expect_identical(.Random.seed, before)
})
