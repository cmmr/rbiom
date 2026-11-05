
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


test_that("boxplot smoke", {

  df <- expect_silent(adiv_table(rare5))
  expect_silent(stats_boxplot(df))
  expect_silent(adiv_boxplot(rare5))
  expect_silent(bdiv_boxplot(rare5))
  expect_silent(taxa_boxplot(rare5))

  # Empty data gives a placeholder plot.
  p <- expect_silent(stats_boxplot(df[0,], x = "Body Site"))
  expect_s3_class(p, 'rbiom_plot')
  expect_true(grepl("No Data to Display", p$code, fixed = TRUE))

  skip_on_cran()

  expect_silent(adiv_boxplot(rare5, x = 'Body Site', flip = TRUE, transform = 'log1p'))
  expect_silent(adiv_boxplot(rare5, x = 'Body Site', flip = TRUE, transform = 'sqrt'))
  expect_silent(adiv_boxplot(rare5, flip = TRUE))             # x = NULL with flip
  expect_silent(adiv_boxplot(rare5, stripe = TRUE))           # x = NULL with stripe

  expect_silent(bdiv_boxplot(rare5, x = 'Body Site', flip = TRUE, test = 'none'))
  expect_silent(bdiv_boxplot(rare5, x = 'Body Site', within = 'Body Site', test = 'none'))

  expect_silent(taxa_boxplot(rare5, x = "Body Site", rank = c("Phylum", "Genus")))
  expect_silent(taxa_boxplot(rare50, layers = "bed", stat.by = 'Body Site', p.top = 3))
  expect_silent(taxa_boxplot(rare50, taxa = 3, transform = 'percent', y.transform = 'none'))

  expect_silent(adiv_boxplot(rare50, x = 'Body Site', stat.by = 'Body Site'))
  expect_silent(adiv_boxplot(rare50, stat.by = 'Sex'))
  expect_silent(adiv_boxplot(rare50, x = 'Sex', facet.by = 'Body Site', scales = 'free', test = 'none'))
  expect_silent(adiv_boxplot(rare50, x = 'Body Site', stat.by = 'Body Site', facet.by = 'Sex', adiv = c('shannon', 'simpson'), test = 'none'))
})


test_that("boxplot code replays", {

  skip_on_cran()

  expect_replays(adiv_boxplot(rare50, x = "Body Site", stat.by = "Body Site"))
  expect_replays(adiv_boxplot(rare50, x = "Sex", flip = TRUE, stripe = TRUE))
  expect_replays(adiv_boxplot(rare50, x = "Sex", layers = "bdp", stat.by = "Sex", p.label = 1))
  expect_replays(adiv_boxplot(rare50, flip = TRUE))
  expect_replays(bdiv_boxplot(rare50[1:20], x = "==Sex", bdiv = "bray", permutations = 99))
  expect_replays(taxa_boxplot(rare50, stat.by = "Body Site", stripe = TRUE, taxa = 3))
  expect_replays(taxa_boxplot(rare50, x = "Body Site", stat.by = "Body Site", rank = c("Phylum", "Genus"), taxa = 4, p.top = 2))
  expect_replays(stats_boxplot(adiv_table(rare50, md = c("Sex", "Age")), x = "Sex", y = "Age", layers = "lv", ci = "sd"))
})


test_that("boxplot statistics routing", {

  skip_on_cran()

  # x == stat.by: pairwise Mann-Whitney U.
  p <- adiv_boxplot(rare50, x = "Body Site", stat.by = "Body Site")
  expect_equal(nrow(p$stats), choose(5, 2))
  expect_true(grepl("Mann-Whitney", p$code, fixed = TRUE))
  expect_true(all(c("labels", "brackets") %in% ls(envir = (function () { e <- new.env(parent = globalenv()); eval(parse(text = p$code), e); e })())))

  # x != stat.by: Kruskal-Wallis per x category.
  p <- adiv_boxplot(rare50, x = "Sex", stat.by = "Body Site", layers = "b")
  expect_true(grepl("Kruskal", p$code, fixed = TRUE))
  expect_equal(nrow(p$stats), 2L)

  # Explicit tests.
  p <- adiv_boxplot(rare50, x = "Sex", test = "kruskal")
  expect_equal(nrow(p$stats), 1L)
  expect_true(grepl("Kruskal", p$code, fixed = TRUE))

  p <- adiv_boxplot(rare50, x = "Sex", test = "emmeans")
  expect_true(grepl("marginal means", p$code, fixed = TRUE))

  expect_null(adiv_boxplot(rare50, x = "Sex", test = "none")$stats)
  expect_null(adiv_boxplot(rare50, x = "Sex", p.label = 0)$stats)
  expect_error(adiv_boxplot(rare50, x = "Sex", test = "mrpp"), "not available")

  # Beta diversity: PERMANOVA on the samples, and the statistics table
  # carries the plot's history.
  p <- bdiv_boxplot(rare50[1:20], x = "==Sex", permutations = 99)
  expect_true(grepl("PERMANOVA", p$code, fixed = TRUE))
  expect_true(grepl("distmat_stats", p$stats$code, fixed = TRUE))
  expect_true(grepl("rarefy(biom", p$stats$code, fixed = TRUE))
  expect_equal(p$stats$.n[[1]], 20L)
  expect_silent(bdiv_boxplot(rare50[1:20], x = "==Sex", test = "mrpp", permutations = 99))
  expect_silent(bdiv_boxplot(rare50[1:20], x = "==Sex", test = "betadisper", permutations = 99))
  expect_error(bdiv_boxplot(rare50[1:20], x = "==Sex", test = "wilcox"), "not valid")

  # Strata: mixed model for pairwise comparisons, blocks for permutations.
  p <- adiv_boxplot(babies, x = "Delivery", strata = "Subject")
  expect_true(grepl("random intercept per Subject", p$code, fixed = TRUE))
  expect_true(grepl("Subject", p$stats$code, fixed = TRUE))
  p <- bdiv_boxplot(babies[1:30], x = "==Delivery", strata = "Subject", permutations = 99)
  expect_true(grepl("within Subject", p$code, fixed = TRUE))

  # Every comparison fails (each facet has one group): NA p-values are
  # dropped before labeling, so there is no "italic(p)== " parse error.
  df <- adiv_table(rare50)
  df <- df[(df$`Body Site` == "Saliva" & df$Sex == "Female") | (df$`Body Site` == "Stool" & df$Sex == "Male"), ]
  p  <- suppressWarnings(stats_boxplot(df, x = "Sex", facet.by = "Body Site", p.label = 1))
  expect_true(all(is.na(p$stats$.p.val)))
  expect_null(attr(p$data, 'stat_labels'))
  expect_silent(ggplot2::ggplot_build(p))

  # Empty statistics table: no label layer, no parse error.
  df <- adiv_table(rare50)
  df$.diversity[df$Sex == "Male"] <- NA
  df <- df[!is.na(df$.diversity), ]
  expect_error(stats_boxplot(df, x = "Sex"), "only one group")
  expect_silent(stats_boxplot(df, x = "Sex", test = "none"))
})


test_that("boxplot provenance", {

  skip_on_cran()

  p <- adiv_boxplot(rare50, x = "Body Site", stat.by = "Body Site", adiv = c("shannon", "simpson"))
  expect_identical(as.character(p$cmd), 'adiv_boxplot(biom, "Body Site", c("shannon", "simpson"))')
  expect_true(grepl("rarefy(biom", p$code, fixed = TRUE))
  expect_true(grepl("data <- adiv_table(biom", p$code, fixed = TRUE))
  expect_identical(as.character(p$stats$history), as.character(p$history))

  p <- bdiv_boxplot(rare5, x = "==Body Site", bdiv = "jaccard")
  expect_identical(as.character(p$cmd), 'bdiv_boxplot(biom, "==Body Site", "jaccard")')

  p <- taxa_boxplot(rare5, rank = "Phylum", taxa = 2)
  expect_identical(as.character(p$cmd), 'taxa_boxplot(biom, rank = "Phylum", taxa = 2)')

  # The same contents with a different history: the code follows the object.
  b2 <- subset(rare50, TRUE)
  p  <- adiv_boxplot(b2, x = "Sex")
  expect_true(grepl("subset(biom, TRUE)", p$code,       fixed = TRUE))
  expect_true(grepl("subset(biom, TRUE)", p$stats$code, fixed = TRUE))
})


test_that("boxplot caching", {

  skip_on_cran()

  old <- options(rbiom.cache_dir = "FALSE")
  on.exit(options(old), add = TRUE)
  p1 <- adiv_boxplot(rare50, x = "Sex")

  options(rbiom.cache_dir = tempfile())
  p2 <- adiv_boxplot(rare50, x = "Sex")
  p3 <- adiv_boxplot(rare50, x = "Sex")

  expect_identical(p1$code, p2$code)
  expect_identical(p2$code, p3$code)
  expect_equal(p2$stats, p3$stats)
})


test_that("boxplot errors", {

  expect_error(adiv_boxplot(rare50, x = "Sex", size = 2),        "not an argument")
  expect_error(adiv_boxplot(rare50, x = "Sex", colour.by = "Sex"), "not an argument")
  expect_error(adiv_boxplot(rare50, x = "Sex", color.by = "Sex"),  "color.by")
  expect_error(adiv_boxplot(rare50, x = "Sex", colors = c(Male = "red", Zzz = "blue")), "Zzz")
  expect_error(adiv_boxplot(rare50, x = "Sex", colors = c(Male = "red")), "Missing")
  expect_error(adiv_boxplot(rare50, x = "Sex", test = "nope"))
  expect_error(stats_boxplot(1:10), "data.frame")

  # Arguments reach the right layers; `b.` means the bar layer.
  p <- adiv_boxplot(rare50, x = "Sex", layers = "b", b.alpha = 0.2, test = "none")
  expect_true(grepl("alpha += 0.2", p$code))
  p <- adiv_boxplot(rare50, x = "Sex", facet.by = "Body Site", scales = "free", test = "none")
  expect_true(grepl('scales += "free"', p$code))
})
