# `$methods`: manuscript-ready methods text with references.

# The paragraph, and the references as one string.
m_text <- function (x) attr(x$methods, 'text',       exact = TRUE)
m_refs <- function (x) paste(attr(x$methods, 'references', exact = TRUE), collapse = "\n")

expect_mentions <- function (txt, ...) {
  for (pattern in c(...))
    expect_match(txt, pattern, fixed = TRUE, label = "methods text")
}

expect_omits <- function (txt, ...) {
  for (pattern in c(...))
    expect_no_match(txt, pattern, fixed = TRUE, label = "methods text")
}

# Arguments that make each covered function run quickly on the test data.
methods_test_args <- function (fn) {
  switch(
    EXPR              = fn,
    adiv_boxplot      = list(biom = rare50, x = "Sex", stat.by = "Sex"),
    adiv_corrplot     = list(biom = rare50, x = "Age"),
    adiv_matrix       = list(biom = rare50),
    bdiv_distmat      = list(biom = rare50),
    adiv_stats        = list(biom = rare50, stat.by = "Sex", regr = "BMI"),
    adiv_table        = list(biom = rare50),
    bdiv_boxplot      = list(biom = rare50, x = "Sex", stat.by = "Sex", permutations = 9),
    bdiv_corrplot     = list(biom = rare50, x = "Age", permutations = 9),
    bdiv_heatmap      = list(biom = rare50),
    bdiv_ord_plot     = list(biom = rare50, stat.by = "Sex", permutations = 9),
    bdiv_ord_table    = list(biom = rare50, stat.by = "Sex", permutations = 9),
    bdiv_stats        = list(biom = rare50, stat.by = "Sex", permutations = 9),
    bdiv_table        = list(biom = rare5),
    distmat_stats     = list(dm = bdiv_distmat(rare50), groups = pull(rare50, "Sex"), permutations = 9),
    distmat_ord_table = list(dm = bdiv_distmat(rare50)),
    rare_corrplot     = list(biom = hmp50),
    rare_multiplot    = list(biom = hmp50),
    rare_stacked      = list(biom = hmp50),
    rarefy            = list(biom = hmp50),
    taxa_boxplot      = list(biom = rare50, x = "Sex", stat.by = "Sex", taxa = 3),
    taxa_corrplot     = list(biom = rare50, x = "Age", taxa = 2),
    taxa_heatmap      = list(biom = rare50),
    taxa_matrix       = list(biom = rare50),
    taxa_stacked      = list(biom = rare50),
    taxa_forest       = list(biom = hmp50, stat.by = "Sex"),
    taxa_da           = list(biom = hmp50, stat.by = "Sex"),
    taxa_stats        = list(biom = rare50, stat.by = "Sex", regr = "BMI", taxa = 2),
    taxa_table        = list(biom = rare50),
    taxa_top          = list(biom = rare50),
    taxa_shared       = list(biom = rare50),
    taxa_upset        = list(biom = rare50, stat.by = "Body Site"),
    NULL )
}


test_that("return value: one string, with the text and references as attributes", {

  m <- adiv_table(rare50)$methods

  expect_s3_class(m, 'rbiom_methods')
  expect_type(m, 'character')
  expect_length(m, 1)

  text <- attr(m, 'text',       exact = TRUE)
  refs <- attr(m, 'references', exact = TRUE)

  expect_type(text, 'character')
  expect_length(text, 1)
  expect_type(refs, 'character')
  expect_true(length(refs) >= 3)

  # as.character() is the full text: paragraph, then the references.
  full <- as.character(m)
  expect_null(attributes(full))
  expect_identical(full, paste0(text, "\n\nReferences\n", paste(refs, collapse = "\n")))

  # Alphabetical, each once.
  expect_identical(refs, refs[order(tolower(refs), method = "radix")])
  expect_false(anyDuplicated(refs) > 0)

  # No markdown, no first person.
  expect_no_match(text, "[`*#]")
  expect_no_match(text, "\\b(we|our|I)\\b")
})


test_that("print() wraps the text and lists the references", {

  skip_on_cran()

  m   <- adiv_table(rare50)$methods
  rlang::local_options(width = 60)
  out <- capture.output(res <- print(m))

  expect_identical(res, m)   # returned invisibly
  expect_true(all(nchar(out) <= 60))
  expect_true("References" %in% out)
  expect_true(any(startsWith(out, "Shannon CE (1948)")))
  expect_false(any(grepl('^\\[1\\]|"', out)))   # cat(), not a quoted string
})


test_that("data preparation comes from the history", {

  skip_on_cran()

  txt <- m_text(rare50)
  expect_mentions(txt, sprintf("rarefied to %s reads per sample", formatC(rare50$depth, format = "d", big.mark = ",")))
  expect_mentions(txt, "random seed 0", sprintf("comprised %i samples", rare50$n_samples))
  expect_mentions(txt, sprintf("rbiom %s", packageVersion("rbiom")), "ecodive")
  expect_mentions(m_refs(rare50), "R Core Team", "Smith DP. rbiom:", "Smith DP. ecodive:")

  # Subsetting, with the expression that was used.
  b <- subset(rare50, Age > 25)
  expect_mentions(m_text(b), "Samples were restricted to those satisfying Age > 25.")
  expect_mentions(m_text(b), sprintf("comprised %i samples", b$n_samples))

  # A seed other than the default.
  expect_mentions(m_text(rarefy(hmp50, depth = 500, seed = 3)), "to 500 reads per sample", "random seed 3")

  # Unrarefied input is described as such; file names lose their directory.
  tf <- tempfile(fileext = ".biom")
  on.exit(unlink(tf), add = TRUE)
  write_biom(hmp5, tf)
  b   <- read_biom(tf)
  txt <- m_text(adiv_table(b))
  expect_mentions(txt, sprintf("Data were imported from %s.", basename(tf)), "unrarefied counts")
  expect_omits(txt, dirname(tf), normalizePath(dirname(tf), winslash = "/"))
  expect_no_match(as.character(b$methods), "[A-Za-z]:[/\\\\]|/tmp/|/var/")

  # Relative abundances and direct edits.
  expect_mentions(m_text(biom_relativize(hmp5)), "Counts were converted to relative abundances.")
  b <- hmp5$clone()
  b$counts <- b$counts[, 1:3]
  expect_mentions(m_text(b), "replaced by direct assignment")
})


test_that("alpha diversity: metrics, test, adjustment, citations", {

  skip_on_cran()

  p   <- adiv_boxplot(rare50, x = "Body Site", adiv = c("shannon", "faith"))
  txt <- m_text(p)
  ref <- m_refs(p)

  expect_mentions(txt, "Shannon diversity index (Shannon 1948)", "Faith's phylogenetic diversity (Faith 1992)")
  expect_mentions(txt, "Wilcoxon rank-sum (Mann-Whitney U) tests", "each pair of Body Site groups")
  expect_mentions(txt, "Benjamini-Hochberg false discovery rate procedure (Benjamini and Hochberg 1995)")
  expect_mentions(txt, "each alpha diversity metric", "95% level", "ggplot2")
  expect_mentions(ref, "Shannon CE (1948)", "Faith DP (1992)", "Mann HB, Whitney DR (1947)", "Wilcoxon F (1945)", "Benjamini Y, Hochberg Y (1995)", "Wickham H (2016)")
  expect_omits(ref, "Anderson", "Bray", "Kruskal", "Oksanen", "Zhou", "Lozupone")

  # The plot's statistics table describes the same analysis, without the graphics.
  txt <- m_text(p$stats)
  expect_mentions(txt, "Shannon diversity index", "Wilcoxon rank-sum")
  expect_omits(txt, "ggplot2")

  # Kruskal-Wallis; no statistics at all.
  p <- adiv_boxplot(rare50, x = "Body Site", test = "kruskal")
  expect_mentions(m_text(p), "Kruskal-Wallis rank-sum test (Kruskal and Wallis 1952)")
  expect_omits(m_refs(p), "Wilcoxon", "Mann HB")

  p <- adiv_boxplot(rare50)
  expect_omits(m_text(p), "Wilcoxon", "Kruskal", "P-values")

  # A table: the metric and the software, nothing else.
  t   <- adiv_table(rare50, adiv = c("chao1", "observed"))
  txt <- m_text(t)
  expect_mentions(txt, "Chao1 (Chao 1984)", "observed features", "rarefied to")
  expect_omits(txt, "ggplot2", "P-values", "tested")
  expect_omits(m_refs(t), "Shannon", "Wickham")
})


test_that("regression models: emmeans, GAM, strata", {

  skip_on_cran()

  t   <- adiv_stats(rare50, regr = "Age", stat.by = "Sex", test = "emtrends", fit = "lm")
  txt <- m_text(t)
  expect_mentions(txt, "a linear model", "slope of alpha diversity against Age", "estimated marginal trends from the emmeans package")
  expect_mentions(m_refs(t), "Lenth RV. emmeans:")
  expect_omits(m_refs(t), "Wood SN")

  p   <- adiv_corrplot(rare50, x = "Age", stat.by = "Sex", test = "emmeans", fit = "gam")
  txt <- m_text(p)
  expect_mentions(txt, "generalized additive model (GAM; Wood 2011)", "smooth term for Age", "mgcv", "estimated marginal means")
  expect_mentions(m_refs(p), "Wood SN (2011)")

  t   <- adiv_stats(rare50, stat.by = "Sex", strata = "Body Site", covariates = "Age")
  txt <- m_text(t)
  expect_mentions(txt, "linear mixed model with a random intercept for each Body Site", "as a function of Sex and Age")
})


test_that("beta diversity: UniFrac variants, PERMANOVA, permutations", {

  skip_on_cran()

  p   <- bdiv_boxplot(rare50, x = "==Body Site", bdiv = c("unweighted_unifrac", "generalized_unifrac"), stat.by = "Body Site", permutations = 99, alpha = 0.25)
  txt <- m_text(p)
  ref <- m_refs(p)

  expect_mentions(txt, "unweighted UniFrac (Lozupone and Knight 2005)", "generalized UniFrac (alpha = 0.25; Chen et al. 2012)")
  expect_mentions(txt, "Only pairs of samples within Body Site were compared.")
  expect_mentions(txt, "PERMANOVA; Anderson 2001", "adonis2 from vegan", "99 permutations", "random seed 0")
  expect_mentions(txt, "each beta diversity metric", "Benjamini-Hochberg")
  expect_mentions(ref, "Lozupone C, Knight R (2005)", "Chen J, Bittinger K", "Anderson MJ (2001)", "Oksanen J")
  expect_omits(ref, "Bray JR", "Lozupone CA, Hamady M", "Shannon", "Gower", "Mielke")

  # Weighted UniFrac cites the 2007 paper; a table has no test.
  t <- bdiv_table(rare5, bdiv = "weighted_unifrac", within = "Sex")
  expect_mentions(m_text(t), "weighted UniFrac (Lozupone et al. 2007)", "within Sex")
  expect_omits(m_text(t), "PERMANOVA", "permutations")
  expect_omits(m_refs(t), "Anderson", "Oksanen")

  # MRPP, dispersion, Mantel.
  t <- bdiv_stats(rare50, stat.by = "Sex", test = "mrpp", permutations = 19)
  expect_mentions(m_text(t), "Bray-Curtis dissimilarity (Bray and Curtis 1957)", "MRPP; Mielke and Berry 2001", "19 permutations")
  expect_omits(m_refs(t), "Anderson")

  t <- bdiv_stats(rare50, stat.by = "Sex", test = "betadisper", permutations = 19)
  expect_mentions(m_text(t), "multivariate dispersion", "Anderson 2006")
  expect_omits(m_refs(t), "Anderson MJ (2001)")

  p <- bdiv_corrplot(rare50, x = "Age", test = "mantel", permutations = 19)
  expect_mentions(m_text(p), "Mantel test (Mantel 1967)", "differences in Age", "19 permutations")

  # Strata restrict the permutations.
  t <- bdiv_stats(rare50, stat.by = "Sex", strata = "Body Site", permutations = 19)
  expect_mentions(m_text(t), "19 permutations within Body Site")

  # Heat map: the metric and the ordering.
  p <- bdiv_heatmap(rare50, bdiv = "jaccard")
  expect_mentions(m_text(p), "Jaccard distance (Jaccard 1912)", "hierarchical clustering")
  expect_omits(m_text(p), "PERMANOVA")
})


test_that("ordination: PCoA and NMDS with stat.by", {

  skip_on_cran()

  p   <- bdiv_ord_plot(rare50, bdiv = "bray", ord = "PCoA", stat.by = "Body Site", permutations = 99)
  txt <- m_text(p)
  ref <- m_refs(p)

  expect_mentions(txt, "Bray-Curtis dissimilarity (Bray and Curtis 1957)")
  expect_mentions(txt, "principal coordinates analysis (PCoA; Gower 1966)")
  expect_mentions(txt, "effect of Body Site on community composition", "PERMANOVA; Anderson 2001", "99 permutations")
  expect_mentions(txt, "envfit from vegan")
  expect_mentions(ref, "Gower JC (1966)", "Anderson MJ (2001)", "Bray JR, Curtis JT (1957)", "Oksanen J")
  expect_omits(txt, "NMDS", "stress")
  expect_omits(ref, "Kruskal JB", "Lozupone")

  p   <- bdiv_ord_plot(rare50, bdiv = "bray", ord = "NMDS", stat.by = "Body Site", permutations = 99, seed = 2)
  txt <- m_text(p)

  stress <- unique(attr(p$data, 'axes')[['.stress']])
  expect_length(stress, 1)

  expect_mentions(txt, "non-metric multidimensional scaling (NMDS; Kruskal 1964", sprintf("stress = %.3f", stress))
  expect_mentions(txt, "random seed 2", "PERMANOVA", "vegan")
  expect_mentions(m_refs(p), "Kruskal JB (1964)")
  expect_omits(m_refs(p), "Gower")

  # The table and its `$stats`.
  t <- bdiv_ord_table(rare50, stat.by = "Sex", test = "mrpp", permutations = 19)
  expect_mentions(m_text(t),       "PCoA; Gower 1966", "MRPP")
  expect_mentions(m_text(t$stats), "MRPP; Mielke and Berry 2001", "19 permutations", "Bray-Curtis")
  expect_omits(m_text(t$stats), "ggplot2")

  # Without stat.by there is no test.
  p <- bdiv_ord_plot(rare50, layers = "p")
  expect_omits(m_text(p), "PERMANOVA", "permutations")
})


test_that("taxa abundance: rank, transform, selection, tests", {

  skip_on_cran()

  p   <- taxa_stacked(rare50, rank = "Phylum", taxa = 4)
  txt <- m_text(p)
  expect_mentions(txt, "summarized at the Phylum rank", "relative abundances", "The 4 most abundant taxa were retained", "hierarchical clustering of euclidean distances")
  expect_omits(txt, "P-values", "tested", "Shannon")
  expect_omits(m_refs(p), "Benjamini", "Bray")

  p   <- taxa_corrplot(rare50, x = "Age", rank = "Phylum", taxa = 2, stat.by = "Sex", test = "emtrends", fit = "lm")
  txt <- m_text(p)
  expect_mentions(txt, "the Phylum rank", "The 2 most abundant taxa", "slope of taxon abundance against Age", "each taxon", "emmeans")
  expect_mentions(m_refs(p), "Lenth RV", "Benjamini Y, Hochberg Y (1995)")

  t   <- taxa_table(rare50, rank = "Genus", taxa = 0.01, transform = "log1p")
  txt <- m_text(t)
  expect_mentions(txt, "the Genus rank", "log(1 + x)-transformed", "at least 1%")
  expect_omits(txt, "tested", "ggplot2")

  t   <- taxa_stats(rare50, stat.by = "Sex", rank = "Phylum", taxa = 3, test = "wilcox")
  expect_mentions(m_text(t), "the Phylum rank", "Taxon abundance was compared", "Wilcoxon rank-sum", "each taxon")
})


test_that("taxa_top, taxa_shared, taxa_upset", {

  skip_on_cran()

  t <- taxa_top(rare50, rank = "Phylum", n = 3, by = "Body Site")
  expect_mentions(m_text(t), "the Phylum rank", "within each group of Body Site", "the 3 most abundant")
  expect_omits(m_text(t), "tested", "P-values")

  expect_mentions(m_text(taxa_top(rare50, n = 1, by = ".sample")), "within each sample")
  expect_mentions(m_text(taxa_top(rare50)), "across all samples, and the 5 most abundant were listed.")

  t <- taxa_shared(rare50, "Body Site", rank = "Phylum", detect = 0.001, prev = 0.9)
  expect_mentions(m_text(t), "the Phylum rank", "relative abundance exceeded 0.1%", "a Body Site group", "at least 90% of its samples")
  expect_omits(m_refs(t), "Lex A")

  t <- taxa_shared(rare50)
  expect_mentions(m_text(t), "when it was observed", "at least one of its samples")

  p <- taxa_upset(rare50, stat.by = "Body Site")
  expect_mentions(m_text(p), "UpSet plot (Lex et al. 2014)", "a Body Site group")
  expect_mentions(m_refs(p), "Lex A, Gehlenborg N")
})


test_that("differential abundance: LinDA", {

  skip_on_cran()

  t   <- taxa_da(hmp50, stat.by = "Sex", rank = "Phylum", prev.cut = 0.2)
  txt <- m_text(t)
  ref <- m_refs(t)

  expect_mentions(txt, "unrarefied counts", "LinDA (Zhou et al. 2022)", "at the Phylum rank", "log abundance ~ Sex")
  expect_mentions(txt, "fewer than 20% of the samples", "95% level")
  expect_mentions(txt, "Benjamini-Hochberg", "across taxa within each comparison")
  expect_mentions(ref, "Zhou H, He K, Chen J, Zhang X (2022)", "Benjamini Y, Hochberg Y (1995)")
  expect_omits(ref, "Oksanen", "Shannon", "Anderson", "Lenth", "Wickham")
  expect_omits(txt, "rarefied to")

  p   <- taxa_forest(hmp50, stat.by = "Sex", rank = "Phylum", p.sig = 0.01)
  txt <- m_text(p)
  expect_mentions(txt, "LinDA (Zhou et al. 2022)", "adjusted p-value of at most 0.01", "ggplot2")
  expect_mentions(m_refs(p), "Zhou H", "Wickham H (2016)")
  expect_mentions(m_text(p$stats), "LinDA (Zhou et al. 2022)")
})


test_that("rarefaction curves", {

  skip_on_cran()

  p      <- rare_corrplot(hmp50, adiv = "observed")
  txt    <- m_text(p)
  depths <- unique(p$data[['.depth']])

  expect_mentions(txt, "unrarefied counts", sprintf("Each sample was rarefied to %i depths", length(depths)))
  expect_mentions(txt, formatC(max(depths), format = "d", big.mark = ","), "observed features at each depth")
  expect_mentions(txt, "ecodive")
  expect_omits(m_refs(p), "Shannon", "Anderson")

  expect_mentions(m_text(rare_stacked(hmp50, rline = 500)), "rarefying to 500 reads per sample")
})


test_that("results rbiom can't describe get the software sentence", {

  skip_on_cran()

  generic <- sprintf("Analyses were performed with rbiom %s", packageVersion("rbiom"))

  m <- rbiom:::methods_text(NULL)
  expect_s3_class(m, 'rbiom_methods')
  expect_match(attr(m, 'text'), generic, fixed = TRUE)

  # A data.frame statistics table: no history, but the test is named.
  df <- data.frame(y = c(1, 3, 2, 5, 4, 6, 8, 7), g = rep(c("a", "b"), 4))
  t  <- stats_table(df, resp = "y", stat.by = "g", test = "wilcox")
  expect_mentions(m_text(t), "The response was compared", "Wilcoxon rank-sum", generic)
  expect_omits(m_text(t), "imported", "rarefied")

  # A command that can't be parsed, or names an unknown function.
  t <- adiv_table(rare50)
  attr(t, 'cmd') <- "not ( valid R"
  expect_match(m_text(t), generic, fixed = TRUE)
  attr(t, 'cmd') <- "some_other_function(biom, 1)"
  expect_match(m_text(t), generic, fixed = TRUE)

  # Recorded commands are parsed, never run.
  attr(t, 'cmd') <- "adiv_table(biom, stop('evaluated'))"
  expect_no_error(m_text(t))

  # A column named `methods` wins, as for `$code`.
  t <- adiv_table(rare50)
  t$methods <- "x"
  expect_identical(t$methods, rep("x", nrow(t)))
})


test_that("every function that advertises `methods` returns non-empty text", {

  skip_on_cran()

  schema <- fn_schema()
  fns    <- names(schema)[vapply(schema, function (e) "methods" %in% unlist(e$returns$has), logical(1))]

  # Every table and plot (and rbiom object) is covered.
  kinds <- vapply(schema, function (e) e$returns$type, character(1))
  expect_setequal(fns, names(kinds)[kinds %in% c("table", "plot", "biom", "matrix", "dist")])
  expect_true(all(c('taxa_top', 'taxa_shared', 'taxa_upset', 'rarefy') %in% fns))

  for (fn in fns) {

    args <- methods_test_args(fn)
    expect(!is.null(args), sprintf("no test arguments for %s()", fn))
    if (is.null(args)) next

    res <- suppressWarnings(suppressMessages(do.call(fn, args)))
    m   <- tryCatch(res$methods, error = function (e) conditionMessage(e))

    expect(inherits(m, 'rbiom_methods'), sprintf("%s()$methods failed: %s", fn, paste(m, collapse = " ")))
    if (!inherits(m, 'rbiom_methods')) next

    text <- attr(m, 'text', exact = TRUE)

    expect(nchar(text) > 80,                                             sprintf("%s()$methods is too short", fn))
    expect(grepl("Analyses were performed with rbiom", text, fixed = TRUE), sprintf("%s()$methods has no software sentence", fn))
    expect(!grepl("NULL|NA\\b|character\\(0\\)|  ", text),                sprintf("%s()$methods has a formatting slip: %s", fn, text))
    expect(!grepl("[A-Za-z]:[/\\\\]", as.character(m)),                   sprintf("%s()$methods has an absolute path", fn))
    expect(length(attr(m, 'references')) >= 2,                           sprintf("%s()$methods has no references", fn))

    # Each in-text citation of the form (Author Year) has a reference.
    refs  <- paste(attr(m, 'references'), collapse = "\n")
    years <- unique(unlist(regmatches(text, gregexpr("[A-Z][A-Za-z]+( et al\\.| and [A-Z][a-z]+)? (19|20)[0-9]{2}", text))))
    for (y in years) {
      author <- sub(" .*$", "", sub("^van der ", "", y))
      year   <- sub("^.* ", "", y)
      expect(
        any(grepl(author, strsplit(refs, "\n")[[1]], fixed = TRUE) & grepl(sprintf("(%s)", year), strsplit(refs, "\n")[[1]], fixed = TRUE)),
        sprintf("%s()$methods cites '%s' without a reference", fn, y) )
    }

    # The statistics tables reachable from the result.
    for (part in intersect(c('stats', 'taxa_stats'), unlist(schema[[fn]]$returns$has))) {
      s <- attr(res, part, exact = TRUE)
      if (is.null(s)) next
      sm <- tryCatch(s$methods, error = function (e) conditionMessage(e))
      expect(inherits(sm, 'rbiom_methods'), sprintf("%s()$%s$methods failed: %s", fn, part, paste(sm, collapse = " ")))
    }
  }
})


test_that("the text is deterministic and leaves the result's code alone", {

  skip_on_cran()

  a <- bdiv_ord_plot(rare50, stat.by = "Sex", permutations = 9)
  b <- bdiv_ord_plot(rare50, stat.by = "Sex", permutations = 9)

  expect_identical(a$methods, b$methods)
  expect_identical(a$methods, a$methods)
  expect_identical(a$stats$methods, b$stats$methods)

  # Asking for the methods changes nothing about the result.
  code <- a$code
  cmd  <- a$cmd
  invisible(a$methods)
  expect_identical(a$code, code)
  expect_identical(a$cmd,  cmd)

  # The description of the input rides on the history lines; it is not
  # part of the code or of `$history`.
  t <- adiv_table(rare50)
  expect_no_match(as.character(t$code), "biom_info|methods")
  expect_named(attr(attr(t, 'history'), 'biom_info'), c('samples', 'otus', 'table'))
  expect_null(attr(t, 'biom_info'))
  expect_null(attr(t$history, 'biom_info'))

  # An object built the same way twice gives the same text.
  expect_identical(rarefy(hmp50)$methods, rarefy(hmp50)$methods)
})


test_that("metric references match ecodive's manual pages", {

  skip_on_cran()

  db <- tools::Rd_db("ecodive")

  # The first reference in a page's references section, as plain text.
  first_ref <- function (id) {
    rd <- db[[paste0(id, ".Rd")]]
    if (is.null(rd)) return (NULL)
    tags <- vapply(rd, function (x) attr(x, "Rd_tag"), "")
    i    <- which(tags == "\\references")
    if (!length(i)) return (NULL)
    txt <- paste(utils::capture.output(tools::Rd2txt(
      rd[i], fragment = TRUE,
      options = list(underline_titles = FALSE, width = 100000L) )), collapse = "\n")
    txt <- sub("^\\s*References:?\\s*", "", txt)
    gsub("\\s+", " ", trimws(strsplit(txt, "\n\\s*\n")[[1]][[1]]))
  }

  for (id in ecodive::list_metrics()$id) {

    eco <- first_ref(id)
    ref <- METRIC_REFS[[id]]

    if (is.null(eco)) {
      expect(is.null(ref), sprintf("%s: ecodive has no reference, rbiom has one", id))
      next
    }
    expect(!is.null(ref), sprintf("%s: ecodive's reference is missing from METRIC_REFS", id))
    if (is.null(ref)) next

    year   <- regmatches(eco, regexpr("\\(\\d{4}\\)", eco))
    author <- sub(",.*$", "", eco)
    doi    <- regmatches(eco, regexpr("10\\.\\d{4,}/[^ <>]+[^ <>.]", eco))

    expect(grepl(year, ref$ref, fixed = TRUE),          sprintf("%s: year %s", id, year))
    expect(startsWith(ref$ref, author),                 sprintf("%s: first author %s", id, author))
    expect(startsWith(ref$cite, author),                sprintf("%s: in-text citation %s", id, ref$cite))
    if (length(doi))
      expect(grepl(doi, ref$ref, fixed = TRUE),         sprintf("%s: doi %s", id, doi))
  }
})
