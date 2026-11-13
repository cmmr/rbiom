#________________________________________________________
# Provenance acceptance tests.
#
# Every rbiom object, table, statistics table, and plot must carry `$code`
# that, when evaluated in a fresh environment, reproduces the object.
#________________________________________________________

# Evaluate `x$code` in a clean environment. Returns a list with the
# environment (`env`) and the value of the final expression (`value`).
replay <- function (x) {
  code <- if (inherits(x, 'rbiom')) x$code else x$code
  expect_s3_class(code, 'rbiom_code')
  expect_false(any(grepl("not reproducible", code, fixed = TRUE)), label = "code is reproducible")
  expect_false(any(grepl("not recorded", code, fixed = TRUE)),     label = "no unrecorded edits")
  expect_false(any(grepl(HEX_RE, code)),                           label = "no hexadecimal numbers in code")
  if (!inherits(x, 'rbiom'))
    expect_false(any(grepl(HEX_RE, attr(x, 'cmd'))),               label = "no hexadecimal numbers in cmd")
  env   <- new.env(parent = globalenv())
  value <- eval(parse(text = code), envir = env)
  list(env = env, value = value)
}

# Compare two data frames ignoring provenance attributes and row names.
strip <- function (df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (a in c('cmd', 'code', 'history', 'data_cmd', 'tbl_sum', 'stats', 'stat_labels'))
    attr(df, a) <- NULL
  rownames(df) <- NULL
  df
}

# Compare the data that ggplot2 would draw for two plots.
built_data <- function (p) {
  lapply(ggplot2::ggplot_build(p)$data, strip)
}

expect_same_plot <- function (replayed, original) {
  if (inherits(original, 'patchwork')) {
    expect_true(inherits(replayed, 'patchwork'))
    expect_equal(length(replayed$patches$plots), length(original$patches$plots))
    for (i in seq_along(original$patches$plots))
      expect_equal(built_data(replayed$patches$plots[[i]]), built_data(original$patches$plots[[i]]), tolerance = 1e-6)
    expect_equal(built_data(replayed), built_data(original), tolerance = 1e-6)
  } else {
    expect_equal(built_data(replayed), built_data(original), tolerance = 1e-6)
  }
}


biom_file <- system.file("extdata", "hmp50.bz2", package = "rbiom")


test_that("rbiom object history replays", {

  b <- read_biom(biom_file)
  b <- subset(b, Age > 25 & `Body Site` != "Stool")
  cutoff <- 22
  b <- subset(b, BMI > cutoff)
  b <- rarefy(b, depth = 1000, seed = 7)
  b <- mutate(b, Old = Age > 35)
  b <- rename(b, Years = "Age")
  b <- b[head(b$samples, 12)]
  b <- subset_taxa(b, Phylum == "Firmicutes" | Phylum == "Bacteroidetes")
  b <- na.omit(b)

  r <- replay(b)
  expect_identical(r$env$biom$hash, b$hash)
  expect_equal(length(b$history), 9L)

  # Direct edits are recorded as non-reproducible comments.
  b2 <- b$clone()
  b2$metadata$Years <- 0
  expect_true(any(grepl("not reproducible", b2$code, fixed = TRUE)))
  expect_false(any(grepl("not reproducible", b$code, fixed = TRUE)))

  # Relativize / rescale / inflate are recorded too.
  b3 <- biom_relativize(b)
  expect_true(b3$relative)
  expect_identical(replay(b3)$env$biom$hash, b3$hash)

  # Merging nests the inputs' histories.
  m <- biom_merge(b[1:4], b[5:8])
  expect_identical(replay(m)$env$biom$hash, m$hash)
})


test_that("cache does not interfere with provenance", {

  b1 <- read_biom(biom_file)
  b2 <- subset(read_biom(biom_file), TRUE)   # same contents, different history

  t1 <- adiv_table(b1)
  t2 <- adiv_table(b2)                       # cache hit on contents

  expect_false(identical(t1$history, t2$history))
  expect_true(grepl("subset(biom, TRUE)", t2$code, fixed = TRUE))
})


test_that("tables replay", {

  b <- rarefy(read_biom(biom_file), depth = 1000, seed = 7)

  check <- function (x) {
    r <- replay(x)
    if (is.data.frame(x))      expect_equal(strip(r$env$data), strip(x))
    else if (inherits(x, 'dist')) expect_equal(as.matrix(r$env$data), as.matrix(x))
    else if (is.matrix(x))     expect_equal(unclass(r$env$data), unclass(x), ignore_attr = TRUE)
    else                       expect_equal(unname(r$env$data), unname(x), ignore_attr = TRUE)
  }

  # A table and a distance matrix everywhere; the rest off CRAN.
  check(adiv_table(b, adiv = c("shannon", "observed"), md = "Sex"))
  check(bdiv_distmat(b, "unweighted_unifrac"))

  skip_on_cran()

  for (x in list(
    adiv_matrix(b),
    adiv_vector(b, "simpson", transform = "rank"),
    bdiv_table(b, bdiv = c("bray", "jaccard"), md = c("Sex", "==Body Site")),
    bdiv_matrix(b, "bray"),
    taxa_table(b, rank = "Genus", taxa = 5, md = "Sex"),
    taxa_matrix(b, rank = "Phylum", taxa = 3),
    taxa_map(b, rank = "Genus", taxa = 4, other = TRUE),
    bdiv_ord_table(b, bdiv = "bray", ord = "PCoA", stat.by = "Body Site", rank = "Genus", taxa = 3) ))
    check(x)
})


test_that("statistics replay", {

  b <- rarefy(read_biom(biom_file), depth = 1000, seed = 7)
  b <- mutate(b, Old = Age > 30)

  check <- function (x) {
    r <- suppressMessages(replay(x))
    s <- if (!is.null(r$env$stats)) r$env$stats else r$value
    expect_equal(strip(s), strip(x), tolerance = 1e-6)
  }

  # One of each kind of test everywhere; the rest off CRAN.
  check(adiv_stats(b, stat.by = "Body Site"))
  check(bdiv_stats(b, stat.by = "Body Site"))

  skip_on_cran()

  for (x in list(
    adiv_stats(b, stat.by = "Sex", test = "wilcox"),
    adiv_stats(b, stat.by = "Body Site", test = "kruskal", split.by = "Sex"),
    adiv_stats(b, regr = "Age"),
    adiv_stats(b, regr = "Age", stat.by = "Sex"),
    adiv_stats(b, regr = "Age", stat.by = "Sex", test = "emtrends", fit = "lm"),
    adiv_stats(b, regr = "Age", stat.by = "Sex", at = c(25, 30)),
    bdiv_stats(b, stat.by = "Sex", split.by = "Old"),
    bdiv_stats(b, regr = "Age"),
    taxa_stats(b, stat.by = "Sex", rank = "Phylum", taxa = 3),
    taxa_stats(b, regr = "BMI", rank = "Phylum", taxa = 3, fit = "lm"),
    distmat_stats(bdiv_distmat(b, "bray"), pull(b, "Body Site")),
    stats_table(adiv_table(b, md = c("Sex", "Age")), stat.by = "Sex", regr = "Age"),
    adiv_stats(b, stat.by = "Body Site", covariates = c("Sex", "Age")),
    adiv_stats(b, regr = "Age", stat.by = "Sex", covariates = "Body Site", test = "emtrends"),
    taxa_stats(b, regr = "BMI", covariates = c("Age", "Sex"), rank = "Phylum", taxa = 3, fit = "log"),
    suppressMessages(taxa_da(b, stat.by = "Sex", rank = "Genus", taxa = 5)),
    suppressMessages(taxa_da(b, regr = "Age", covariates = "Sex", split.by = "Old", rank = c("Phylum", "Genus"))) ))
    check(x)
})


test_that("plots replay", {

  b <- rarefy(read_biom(biom_file), depth = 1000, seed = 7)
  b <- mutate(b, Old = Age > 30)

  # Replays re-emit the same warnings as the original call. The statistics
  # code shown beside a plot is checked for hexadecimal numbers too.
  check <- function (p) {
    r <- suppressWarnings(replay(p))
    expect_same_plot(r$value, p)
    if (!is.null(p$stats))
      expect_false(any(grepl(HEX_RE, p$stats$code)), label = "no hexadecimal numbers in stats code")
  }

  # A plot with statistics and one without everywhere; the rest off CRAN.
  check(adiv_boxplot(b, x = "Body Site", stat.by = "Body Site"))
  check(taxa_stacked(b, rank = "Phylum", taxa = 5))

  skip_on_cran()

  # At least one call of every plotting function. The heatmap's values need
  # 16-17 digits, but no two of its row or column distances tie: the code
  # re-runs hclust(), and tied distances (as in seq_len(20) / 3) make the
  # order hinge on the last bit, which R 4.3 on macOS (Apple silicon) can
  # read back one unit off (see num_literal()).
  mtx <- matrix(seq_len(20)^2 / 7, 4, dimnames = list(letters[1:4], LETTERS[1:5]))

  for (p in list(
    adiv_boxplot(b, "Body Site", layers = c("box", "dot")),
    stats_boxplot(adiv_table(b, md = c("Sex", "Age")), x = "Sex", y = "Age", stat.by = "Sex"),
    stats_corrplot(adiv_table(b, md = c("Sex", "Age")), "Age", stat.by = "Sex", layers = "tc", fit = "lm"),
    plot_heatmap(mtx = mtx, tracks = list(z = c(A = 1/3, B = 2/3, C = 0.1 + 0.2, D = 1, E = 2))),
    adiv_boxplot(b, x = "Sex", stat.by = "Sex", flip = TRUE, stripe = TRUE),
    adiv_boxplot(b, x = "Sex", stat.by = "Sex", facet.by = "Old", adiv = c("shannon", "simpson")),
    # Mid vagina has only one Sex, which is reported in a summary warning.
    local({
      expect_warning(
        p <- adiv_boxplot(b, x = "Body Site", layers = "bdp", stat.by = "Sex", p.label = 1),
        "could not be computed" )
      p }),
    adiv_boxplot(b, flip = TRUE),
    adiv_corrplot(b, x = "Age", stat.by = "Sex", layers = "tcp"),
    adiv_corrplot(b, x = "Age", fit = "lm", layers = "tc"),
    bdiv_boxplot(b, x = "==Body Site", bdiv = "bray"),
    bdiv_corrplot(b, x = "Age", stat.by = "Sex"),
    bdiv_ord_plot(b, stat.by = "Body Site"),
    bdiv_ord_plot(b, stat.by = "Body Site", facet.by = "Sex", rank = "Genus", taxa = 3),
    bdiv_heatmap(b, bdiv = "bray", tracks = "Sex"),
    bdiv_heatmap(b, bdiv = c("bray", "jaccard"), tracks = "Body Site"),
    taxa_boxplot(b, stat.by = "Body Site", stripe = TRUE, taxa = 3),
    taxa_corrplot(b, x = "Age", taxa = 2, layers = "ptc", fit = "lm"),
    taxa_heatmap(b, rank = "Phylum", taxa = 5, tracks = "Sex"),
    taxa_heatmap(b, rank = "Phylum", taxa = 1),
    taxa_stacked(b, rank = c("Phylum", "Genus"), taxa = 4),
    suppressMessages(taxa_forest(b, stat.by = "Body Site", rank = "Phylum", p.top = 3)),
    suppressMessages(taxa_forest(b, regr = "Age", rank = "Genus", p.top = 0.99)),
    rare_corrplot(b, stat.by = "Sex"),
    rare_stacked(b),
    rare_multiplot(b, stat.by = "Body Site") ))
    check(p)
})


test_that("hostile names cannot inject code", {

  # Data-originating strings (sample names, metadata columns and levels, taxa
  # names, ids, paths) must appear in generated code only as quoted data.
  inj <- 'x"); stop("INJECTED"); ("y'

  b <- hmp50[1:12]
  b$samples <- paste0(b$samples, c(inj, ' `\\"#', "\t\\n"))
  b$id      <- inj
  b$comment <- paste0("# ", inj)
  tx <- b$taxonomy
  tx$Genus  <- paste0(as.character(tx$Genus), inj)
  tx$Phylum <- paste0(as.character(tx$Phylum), "`'\\")
  b$taxonomy <- tx
  md <- b$metadata
  md[[inj]] <- rep(c(inj, "ok"), length.out = nrow(md))
  b$metadata <- md

  dir  <- file.path(tempfile(pattern = "it's a dir "), "sub dir")
  dir.create(dir, recursive = TRUE)
  file <- file.path(dir, paste0("biom ", "'quote", ".biom"))
  write_biom(b, file, format = "json")

  b <- read_biom(file)
  b <- eval(bquote(subset(b, .(as.name(inj)) == .(inj) | Sex == "Male")))
  b <- mutate(b, !!inj := paste0(.data[[inj]], "!"))
  b <- rarefy(b, depth = 50, seed = 1)

  outputs <- list(
    b,
    adiv_table(b, md = inj),
    taxa_table(b, rank = "Genus", taxa = 3, md = inj),
    adiv_stats(b, stat.by = inj, test = "wilcox"),
    adiv_stats(b, regr = "Age", covariates = inj, fit = "lm"),
    adiv_boxplot(b, x = inj, stat.by = inj),
    taxa_stacked(b, rank = "Genus", taxa = 3),
    taxa_boxplot(b, rank = "Phylum", taxa = 2, x = inj),
    suppressMessages(taxa_da(b, stat.by = inj, rank = "Genus", prev.cut = 0)),
    suppressMessages(taxa_forest(b, stat.by = inj, rank = "Genus", prev.cut = 0, p.top = Inf)) )

  for (x in outputs) {
    code <- x$code
    expect_true(grepl("INJECTED", code, fixed = TRUE))          # recorded ...
    r <- tryCatch(replay(x), error = function (e) e)
    expect_false(inherits(r, 'error'), info = conditionMessage(r)) # ... but as data
    if (inherits(x, 'rbiom'))         expect_identical(r$env$biom$hash, x$hash)
    else if (is.data.frame(x))        expect_equal(strip(if (!is.null(r$env$stats)) r$env$stats else r$env$data), strip(x))
    else                              expect_same_plot(r$value, x)
  }
})


test_that("provenance helpers", {
  expect_error(P('wmp6CWP7kx::q8jc'))

  pi <- prov_inline(quote(Age > cutoff & Sex == "Male"), env = list2env(list(cutoff = 30)), exclude = c("Age", "Sex"))
  expect_identical(pi$expr, quote(Age > 30 & Sex == "Male"))
  expect_true(pi$reproducible)

  big <- 1:1000
  pi <- prov_inline(quote(Age %in% big), env = environment(), exclude = "Age")
  expect_false(pi$reproducible)
})


test_that("numbers in generated code are decimal and parse back", {

  # R's parser reads 16-17 digit literals exactly only where a long double
  # is wider than a double. On macOS (Apple silicon) 15.299999999999999
  # (9 * 1.7) can come back one unit in the last place off.
  expect_same_double <- function (a, b) {
    if (isTRUE(.Machine$longdouble.digits > 53)) expect_identical(a, b)
    else expect_equal(a, b, tolerance = 4 * .Machine$double.eps)
  }

  # The fewest of 15, 16, or 17 significant digits that read back as the
  # same double, and never hexadecimal: people read this code.
  x   <- c(seq_len(12) * 1.7, 0.1 + 0.2, 1/3, 1e5, 1e-30, 1e300, -2.5, 0, Inf, -Inf, NaN)
  lit <- num_literal(x)
  expect_false(any(grepl(HEX_RE, lit)))
  expect_true(all(nchar(sub("^0+", "", gsub("[^0-9]", "", sub("e.*$", "", lit)))) <= 17))
  expect_same_double(vapply(lit, function (l) eval(str2lang(l)), 0, USE.NAMES = FALSE), x)
  expect_identical(num_literal(9 * 1.7), "15.299999999999999")
  expect_identical(num_literal(c(1.7, 2, 0.5, 0.1 + 0.2)), c("1.7", "2", "0.5", "0.30000000000000004"))
  expect_identical(num_literal(1/3), "0.3333333333333333")
  expect_identical(num_literal(NA_real_), "NA_real_")
  expect_silent(expect_identical(num_literal(c(0, NA)), c("0", "NA")))

  # Embedded vectors and logged calls; names are still quoted by deparse().
  v <- c(`it's` = 9 * 1.7, "a`b" = 1.5, `c(1)` = 1/3, ok = 2)
  expect_false(grepl(HEX_RE, deparse_vec(v)))
  expect_same_double(eval(str2lang(deparse_vec(v))), v)

  cl     <- quote(f(biom, keep = c("a", "b"), n = 3L))
  cl$w   <- c(a = 0.1 + 0.2, b = 1)
  cl$lvl <- 1/3
  txt    <- paste(deparse_exact(cl), collapse = "")
  back   <- str2lang(txt)
  expect_false(grepl(HEX_RE, txt))
  expect_same_double(eval(back$w), cl$w)
  expect_same_double(back$lvl, cl$lvl)
  expect_identical(back$n, 3L)

  # $cmd too.
  expect_identical(as.args(list(1/3)), "0.3333333333333333")
  expect_same_double(eval(str2lang(as.args(list(1/3)))), 1/3)
})


test_that("boxplot p-value positions are short decimals in $code", {

  # The atima report: label heights were written in hexadecimal.
  b <- rarefy(read_biom(biom_file), depth = 1183)
  p <- adiv_boxplot(b, "Body Site", layers = c("box", "dot"))
  expect_false(grepl(HEX_RE, p$code))
  lbl <- Filter(function (l) is.data.frame(l$data) && '.label' %in% names(l$data), p$layers)
  expect_length(lbl, 1L)
  expect_true(all(nchar(sub("^0+", "", gsub("[^0-9]", "", num_literal(lbl[[1]]$data$.y)))) <= 10))

  # pretty() breaks such as 0.6000000000000001 are rounded too.
  p <- bdiv_boxplot(b, x = "==Body Site", bdiv = "bray")
  expect_false(grepl(HEX_RE, p$code))
  expect_false(grepl("0.6000000000000001", p$code, fixed = TRUE))
})
