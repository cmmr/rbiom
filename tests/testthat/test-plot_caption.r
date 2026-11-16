
test_that("input_table_caption describes the table", {

  ss <- sample_sums(hmp50)

  # Rarefied: one depth.
  expect_identical(
    input_table_caption(rare50),
    sprintf("Rarefied to %s reads per sample.", formatC(rare50$depth, format = "d", big.mark = ",")) )

  # Unrarefied: the range of sample sums.
  expect_identical(
    input_table_caption(hmp50),
    sprintf(
      "Unrarefied counts (%s to %s reads per sample).",
      formatC(min(ss), format = "d", big.mark = ","),
      formatC(max(ss), format = "d", big.mark = ",") ))

  # Relative abundances.
  expect_identical(input_table_caption(biom_relativize(hmp50)), "Relative abundances.")

  # Thousands separators.
  expect_match(input_table_caption(hmp50), "[0-9],[0-9]{3} reads per sample")

  # A function that rarefies internally says what it analysed.
  expect_identical(input_table_caption(hmp50, depth = 1000), "Rarefied to 1,000 reads per sample.")
  expect_identical(
    input_table_caption(hmp50, depth = c(5, 500, 2500)),
    "Rarefied to 3 depths, from 5 to 2,500 reads per sample." )

  # Always one string; a pure function of its arguments.
  for (b in list(hmp50, rare50, biom_relativize(hmp50), hmp5)) {
    txt <- input_table_caption(b)
    expect_true(is.character(txt) && length(txt) == 1 && !is.na(txt))
    expect_identical(txt, input_table_caption(b))
  }

  # Non-integer values that aren't relative abundances.
  half <- as_rbiom(list(counts = as.matrix(hmp5$counts) / 2))
  expect_match(input_table_caption(half), "^Unnormalized abundances \\(sample sums from ")

  pct <- as_rbiom(list(counts = 100 * as.matrix(biom_relativize(hmp5)$counts)))
  expect_identical(input_table_caption(pct), "Normalized abundances (each sample sums to 100).")
})


test_that("caption helpers", {

  # Only generalized UniFrac has a parameter its name doesn't state.
  expect_null(bdiv_caption("bray"))
  expect_null(bdiv_caption(c("bray", "weighted_unifrac", "unweighted_unifrac")))
  expect_identical(bdiv_caption("generalized_unifrac", 0.25), "Generalized UniFrac with alpha = 0.25.")
  expect_identical(bdiv_caption(c("bray", "generalized_unifrac")), "Generalized UniFrac with alpha = 0.5.")

  expect_null(within_between_caption())
  expect_identical(within_between_caption("Sex"), "Comparing samples within Sex.")
  expect_identical(within_between_caption(c("Sex", "Sex")), "Comparing samples within Sex.")
  expect_identical(
    within_between_caption("Sex", c("Body Site", "Age")),
    "Comparing samples within Sex and between Body Site, Age." )

  expect_identical(clust_caption("complete", "euclidean"), "'complete' clustering of euclidean distances")
  expect_identical(clust_caption("average", stats::dist(1:3), "the plotted distances"), "'average' clustering of the plotted distances")
  expect_identical(clust_caption(stats::hclust(stats::dist(1:3)), "euclidean"), "the given clustering")

  # One line per argument; sentences in one argument share a line.
  expect_null(caption_lines())
  expect_null(caption_lines(NULL, character(), ""))
  expect_identical(caption_lines("A.", NULL, c("B.", NULL, "C.")), "A.\nB. C.")
})


test_that("every plot's caption states the input table", {

  skip_on_cran()

  rare  <- "^Rarefied to [0-9,]+ reads per sample\\."
  unrar <- "^Unrarefied counts \\([0-9,]+ to [0-9,]+ reads per sample\\)\\."

  cap <- function (p) p$labels$caption

  expect_match(cap(adiv_boxplot(rare50, x = "Sex")),                 rare)
  expect_match(cap(adiv_boxplot(hmp50)),                             unrar)
  expect_match(cap(bdiv_boxplot(rare50, x = "Sex")),                 rare)
  expect_match(cap(taxa_boxplot(rare50, taxa = 2)),                  rare)
  expect_match(cap(taxa_boxplot(biom_relativize(hmp50), taxa = 2)),  "^Relative abundances\\.")
  expect_match(cap(adiv_corrplot(rare50, "Age")),                    rare)
  expect_match(cap(bdiv_corrplot(rare50, "Age", permutations = 9)),  rare)
  expect_match(cap(taxa_corrplot(rare50, "Age", taxa = 2)),          rare)
  expect_match(cap(bdiv_ord_plot(rare50, rank = NULL)),              rare)
  expect_match(cap(bdiv_heatmap(rare50[1:8])),                       rare)
  expect_match(cap(taxa_heatmap(rare50[1:8], taxa = 5)),             rare)
  expect_match(cap(taxa_stacked(hmp50)),                             unrar)
  expect_match(cap(rare_corrplot(hmp50)),   "^Rarefied to 10 depths, from 5 to [0-9,]+ reads per sample\\. Dashed line: [0-9,]+ reads per sample\\.$")
  expect_match(cap(rare_stacked(hmp50)),    "^Rarefaction depth \\(dashed line\\): [0-9,]+ reads per sample\\.$")

  # The statistics follow on their own line(s).
  expect_match(cap(adiv_boxplot(rare50, x = "Sex")),                "\nPairwise Mann-Whitney U tests")
  expect_match(cap(bdiv_corrplot(rare50, "Age", permutations = 9)), "\nMantel test of the distances against differences in Age \\(9 permutations\\)\\.")
  expect_match(cap(adiv_corrplot(rare50, "Age")),                   "\nGAM; p-values test for a non-flat trend")
  expect_match(cap(adiv_corrplot(rare50, "Age", stat.by = "Sex", fit = "lm", test = "emtrends")), "\nLinear model; p-values compare the groups' slopes\\.")

  # Generalized UniFrac's alpha.
  expect_match(cap(bdiv_boxplot(rare50, x = "Sex", bdiv = "generalized_unifrac", alpha = 0.25, test = "none")), "Generalized UniFrac with alpha = 0\\.25\\.")
  expect_no_match(cap(bdiv_boxplot(rare50, x = "Sex", bdiv = "weighted_unifrac", test = "none")), "alpha")

  # Heat maps: how the rows and columns were ordered.
  expect_match(cap(bdiv_heatmap(rare50[1:8])),           "Samples ordered by 'complete' clustering of the plotted distances\\.")
  expect_match(cap(taxa_heatmap(rare50[1:8], taxa = 5)), "Taxa and samples ordered by 'complete' clustering of euclidean distances\\.")
  expect_match(cap(taxa_heatmap(rare50[1:8], taxa = 5, clust = c("average", NA), rescale = "rows")), "rescaled to 0 - 1 within each taxon\\.\nTaxa ordered by 'average' clustering")

  m <- matrix(c(1, 5, 2, 8, 3, 1, 9, 4, 6), nrow = 3, dimnames = list(letters[1:3], LETTERS[1:3]))
  expect_identical(cap(plot_heatmap(m)), "Rows and columns ordered by 'complete' clustering of euclidean distances.")
  expect_null(cap(plot_heatmap(m, clust = NA)))

  # Captions can be turned off, or replaced.
  expect_null(cap(adiv_boxplot(rare50, x = "Sex", caption = FALSE)))
  expect_null(cap(bdiv_ord_plot(rare50, rank = NULL, caption = FALSE)))
  expect_null(cap(rare_corrplot(hmp50, caption = FALSE)))
  expect_null(cap(taxa_heatmap(rare50[1:8], taxa = 5, labs.caption = NULL)))
  expect_null(cap(taxa_stacked(hmp50, labs.caption = NULL)))
  expect_identical(cap(bdiv_heatmap(rare50[1:8], labs.caption = "Mine.")), "Mine.")

  # At most a few lines.
  for (p in list(
    bdiv_boxplot(rare50, x = "==Sex", stat.by = "!=Body Site", bdiv = "generalized_unifrac", permutations = 9),
    taxa_forest(hmp50, regr = "Age", rank = "Phylum"),
    bdiv_ord_plot(rare50, ord = "NMDS", stat.by = "Sex", rank = NULL, permutations = 9) ))
    expect_lte(length(strsplit(cap(p), "\n", fixed = TRUE)[[1]]), 4)
})


test_that("bdiv plots caption `within` once", {

  skip_on_cran()

  # `stat.by` defaults to `x`, so "==Sex" used to be listed twice.
  p <- bdiv_boxplot(rare50, x = "==Sex", permutations = 9)
  expect_match(p$labels$caption, "Comparing samples within Sex.", fixed = TRUE)
  expect_no_match(p$labels$caption, "Sex, Sex", fixed = TRUE)

  p <- bdiv_corrplot(rare50, "Age", stat.by = "==Sex", facet.by = "!=Body Site", permutations = 9)
  expect_match(p$labels$caption, "Comparing samples within Sex and between Body Site.", fixed = TRUE)
})


test_that("taxa_forest's caption names the p-value adjustment", {

  skip_on_cran()

  p <- taxa_forest(hmp50, regr = "Age", rank = "Phylum")
  expect_match(p$labels$caption, "smallest FDR (Benjamini-Hochberg) adjusted p-values", fixed = TRUE)

  p <- taxa_forest(hmp50, regr = "Age", rank = "Phylum", p.adj = "holm", p.top = Inf)
  expect_match(p$labels$caption, "p-values adjusted with the Holm method.", fixed = TRUE)

  p <- taxa_forest(hmp50, regr = "Age", rank = "Phylum", p.adj = "none", p.top = 5)
  expect_match(p$labels$caption, "the 5 taxa with the smallest p-values", fixed = TRUE)
  expect_no_match(p$labels$caption, "adjusted", fixed = TRUE)
})


test_that("deprecation warnings name the function the user called", {

  skip_on_cran()

  old <- options(lifecycle_verbosity = "warning")
  on.exit(options(old), add = TRUE)

  expect_warning(bdiv_distmat(rare50, bdiv = "unifrac", weighted = TRUE), "argument of `bdiv_distmat()`", fixed = TRUE)

  # Validated inside with(params, ...), where the caller is eval().
  expect_warning(
    bdiv_boxplot(rare50, x = "Sex", bdiv = "unifrac", weighted = FALSE, test = "none"),
    "`weighted` argument of `bdiv_boxplot()`", fixed = TRUE )

  expect_warning(
    bdiv_corrplot(rare50, "Age", bdiv = "unifrac", weighted = FALSE, test = "none"),
    "`weighted` argument of `bdiv_corrplot()`", fixed = TRUE )
})


test_that("captions are literal strings in `$code` and replay", {

  skip_on_cran()

  plots <- list(
    adiv_boxplot(rare50, x = "Sex"),
    bdiv_boxplot(rare50, x = "==Sex", bdiv = "generalized_unifrac", alpha = 0.25, permutations = 9),
    taxa_corrplot(rare50, "Age", taxa = 2),
    rare_corrplot(hmp50),
    rare_stacked(hmp50),
    bdiv_ord_plot(hmp50, ord = "NMDS", stat.by = "Sex", rank = NULL, permutations = 9),
    bdiv_heatmap(rare50[1:8]),
    taxa_heatmap(rare50[1:8], taxa = 5),
    taxa_stacked(hmp50, label.by = "Body Site"),
    taxa_forest(hmp50, regr = "Age", rank = "Phylum") )

  for (p in plots) {

    caption <- p$labels$caption
    expect_true(is.character(caption) && length(caption) == 1)

    # No helper calls in the code: the text itself.
    expect_false(grepl("input_table_caption|caption_lines", p$code))
    expect_true(grepl(strsplit(caption, "\n", fixed = TRUE)[[1]][[1]], p$code, fixed = TRUE))

    r <- eval(parse(text = p$code), envir = new.env(parent = globalenv()))
    expect_identical(r$labels$caption, caption)
  }
})


test_that("caption = FALSE removes the caption from heatmaps and stacked plots", {

  biom <- rarefy(hmp50)
  mtx  <- taxa_matrix(biom, rank = "Phylum", taxa = 4, sparse = FALSE)

  cap <- function (p) p$labels$caption

  expect_false(is.null(cap(taxa_stacked(biom))))
  expect_false(is.null(cap(rare_stacked(hmp50))))
  expect_false(is.null(cap(bdiv_heatmap(biom))))
  expect_false(is.null(cap(taxa_heatmap(biom))))
  expect_false(is.null(cap(plot_heatmap(mtx))))

  expect_null(cap(taxa_stacked(biom,  caption = FALSE)))
  expect_null(cap(rare_stacked(hmp50, caption = FALSE)))
  expect_null(cap(bdiv_heatmap(biom,  caption = FALSE)))
  expect_null(cap(taxa_heatmap(biom,  caption = FALSE)))
  expect_null(cap(plot_heatmap(mtx,   caption = FALSE)))

  expect_error(taxa_stacked(biom, caption = "no"))
  expect_error(plot_heatmap(mtx,  caption = NA))
})
