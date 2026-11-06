
# hmp50 is unrarefied; taxa_da() expects counts.
saliva_stool <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool"))
three_sites  <- subset(hmp50, `Body Site` %in% c("Saliva", "Stool", "Buccal mucosa"))


test_that("mean-shift mode matches modeest::mlv(method = 'meanshift')", {

  # Reference values from modeest 2.4.0.
  expect_equal(da_mode(c(-3.2, -0.5, -0.1, 0, 0.05, 0.1, 0.2, 0.3, 2.5, 4)), 0.0774484077671091, tolerance = 1e-12)
  expect_equal(da_mode(seq(-1, 1, length.out = 21)^3 + 0.25), 0.25, tolerance = 1e-12)
  expect_equal(da_mode(c(1, 1, 1, 2, 2, 2, 5, 9, -4, 1.5)), 1.5, tolerance = 1e-12)
  expect_identical(da_mode(3), 3)
  expect_identical(da_mode(numeric()), NA_real_)
})


test_that("LinDA matches MicrobiomeStat::linda()", {

  da <- taxa_da(saliva_stool, stat.by = "Body Site", rank = "Genus")

  expect_s3_class(da, 'rbiom_tbl')
  expect_identical(
    names(da),
    c("Body Site", ".taxa", ".lfc", ".se", ".lower", ".upper", ".t.ratio", ".df", ".p.val", ".adj.p", ".n") )
  expect_identical(levels(da$`Body Site`), "Stool - Saliva")
  expect_false(is.unsorted(da$.p.val))
  expect_true(all(da$.n == 20L))

  # Reference values from MicrobiomeStat 1.4:
  # linda(taxa_matrix(saliva_stool, "Genus"), md, "~grp", prev.filter = 0.1)
  ref <- data.frame(
    row.names = c("Bacteroides", "Streptococcus", "Neisseria"),
    lfc  = c(8.71707819361206, -8.72967330054562, -7.82223620183642),
    se   = c(0.791807615124557, 0.518397714750094, 0.507108520156645),
    p    = c(1.99466645427255e-09, 1.83294115604018e-12, 8.06775478726673e-12),
    padj = c(1.66146515641917e-08, 4.49070583229844e-11, 1.07814541248019e-10) )

  x <- da[match(rownames(ref), da$.taxa), ]
  expect_equal(x$.lfc,   ref$lfc,  tolerance = 1e-10)
  expect_equal(x$.se,    ref$se,   tolerance = 1e-10)
  expect_equal(x$.p.val, ref$p,    tolerance = 1e-8)
  expect_equal(x$.adj.p, ref$padj, tolerance = 1e-8)
  expect_true(all(x$.df == 18))

  # Bias-corrected: the fold changes' mode is zero.
  expect_equal(da_mode(sqrt(20) * da$.lfc), 0, tolerance = 1e-6)

  # Confidence intervals.
  q <- qt(0.975, 18)
  expect_equal(da$.lower, da$.lfc - q * da$.se)
  expect_equal(da$.upper, da$.lfc + q * da$.se)
})


test_that("LinDA with a numeric variable, covariates, and zero imputation", {

  # Sequencing depth is associated with Age, so zeros are imputed.
  # Reference values from MicrobiomeStat 1.4 with adaptive = FALSE and
  # zero.handling = "imputation", converted from per-SD to per-year.
  da <- taxa_da(hmp50, regr = "Age", covariates = "Sex", rank = "Genus")

  expect_false(hasName(da, 'Age'))
  x <- da[match(c("Fusobacterium", "Neisseria"), da$.taxa), ]
  expect_equal(x$.lfc,   c(0.0738758338608936, -0.0279397021892972), tolerance = 1e-10)
  expect_equal(x$.se,    c(0.104191041009793,   0.118243083263858),  tolerance = 1e-10)
  expect_equal(x$.p.val, c(0.48180086079324,    0.814234107882185),  tolerance = 1e-8)
  expect_equal(x$.adj.p, c(0.987763959611515,   0.987763959611515),  tolerance = 1e-8)

  # Fold changes are per unit of `regr`.
  b2 <- mutate(hmp50, Age10 = Age * 10)
  da10 <- taxa_da(b2, regr = "Age10", covariates = "Sex", rank = "Genus")
  da10 <- da10[match(da$.taxa, da10$.taxa), ]
  expect_equal(da10$.lfc * 10, da$.lfc, tolerance = 1e-8)
  expect_equal(da10$.p.val, da$.p.val, tolerance = 1e-8)
})


test_that("LinDA pairwise comparisons among three groups", {

  da <- taxa_da(three_sites, stat.by = "Body Site", rank = "Genus")

  expect_identical(
    levels(da$`Body Site`),
    c("Saliva - Buccal mucosa", "Stool - Buccal mucosa", "Stool - Saliva") )

  # Each comparison is corrected by its own mode.
  for (cmp in levels(da$`Body Site`)) {
    x <- da[da$`Body Site` == cmp, ]
    expect_equal(da_mode(sqrt(30) * x$.lfc), 0, tolerance = 1e-6)
    # p-values are adjusted within each comparison.
    expect_equal(x$.adj.p, p.adjust(x$.p.val, "fdr"))
  }

  # The comparisons don't depend on which group is the reference level:
  # reversing the order flips the signs.
  lfc <- function (d, cmp) { x <- d[d$`Body Site` == cmp, ]; stats::setNames(x$.lfc, x$.taxa) }
  rev <- mutate(three_sites, `Body Site` = factor(`Body Site`, levels = c("Stool", "Saliva", "Buccal mucosa")))
  dr  <- taxa_da(rev, stat.by = "Body Site", rank = "Genus")
  for (pair in list(c("Saliva", "Buccal mucosa"), c("Stool", "Buccal mucosa"), c("Stool", "Saliva"))) {
    a  <- lfc(da, paste(pair[[1]], "-", pair[[2]]))
    b  <- lfc(dr, paste(pair[[2]], "-", pair[[1]]))
    expect_equal(unname(b[names(a)]), -unname(a), tolerance = 1e-8)
  }
})


test_that("taxa, ranks, split.by, and NA handling", {

  all <- taxa_da(saliva_stool, stat.by = "Body Site", rank = "Genus")
  top <- taxa_da(saliva_stool, stat.by = "Body Site", rank = "Genus", taxa = 5)
  expect_equal(nrow(top), 5L)
  # Adjusted for every tested taxon, not just the reported ones.
  expect_equal(top$.adj.p, all$.adj.p[match(top$.taxa, all$.taxa)])

  named <- taxa_da(saliva_stool, stat.by = "Body Site", rank = "Genus", taxa = c("Bacteroides", "Nope"))
  expect_identical(named$.taxa, "Bacteroides")

  multi <- taxa_da(saliva_stool, stat.by = "Body Site", rank = c("Phylum", "Genus"))
  expect_identical(levels(multi$.rank), c("Phylum", "Genus"))
  expect_equal(
    multi[multi$.rank == "Genus", ][['.lfc']],
    all$.lfc[match(multi[multi$.rank == "Genus", ][['.taxa']], all$.taxa)] )

  split <- taxa_da(three_sites, stat.by = "Sex", split.by = "Body Site", rank = "Phylum")
  expect_true(hasName(split, 'Body Site'))
  expect_identical(levels(split$Sex), "Male - Female")

  b <- saliva_stool$clone()
  md <- b$metadata
  md$Sex[1:3] <- NA
  b$metadata <- md
  expect_message(
    x <- taxa_da(b, stat.by = "Body Site", covariates = "Sex", rank = "Phylum"),
    "Dropped 3 of 20 samples" )
  expect_true(all(x$.n == 17L))
})


test_that("argument checks", {

  expect_error(taxa_da(saliva_stool, rank = "Genus"), "neither")
  expect_error(taxa_da(saliva_stool, stat.by = "Sex", regr = "Age"), "both")
  expect_error(taxa_da(saliva_stool, stat.by = "Sex", covariates = "Sex"), "different fields")
  expect_error(taxa_da(saliva_stool, stat.by = "Age"), "not categorical")
  expect_error(taxa_da(saliva_stool, regr = "Sex"), "not numeric")
  expect_error(taxa_da(saliva_stool, stat.by = "Sex", method = "deseq"), "method")
  expect_error(taxa_da(biom_relativize(saliva_stool), stat.by = "Sex"), "relative abundances")
  expect_error(taxa_da(subset(hmp50, Sex == "Male"), stat.by = "Sex"), "fewer than two groups")
})


test_that("LinDA with repeated measures (lmerTest)", {

  skip_if_not_installed("lmerTest")
  skip_on_cran()

  bb <- subset(babies, `Age (days)` %in% c(10, 11, 12, 100, 101, 102))
  bb <- mutate(bb, Phase = ifelse(`Age (days)` < 50, "early", "late"))

  da <- suppressMessages(taxa_da(bb, stat.by = "Phase", strata = "Subject ID", rank = "Phylum"))

  expect_identical(levels(da$Phase), "late - early")
  expect_true(all(is.finite(da$.p.val)))
  expect_true(all(da$.df > 0))
  expect_match(attr(da, 'tbl_sum')[['Model']], "(1 | `Subject ID`)", fixed = TRUE)
})


test_that("ANCOM-BC2 matches ANCOMBC::ancombc2()", {

  skip_if_not_installed("ANCOMBC")
  skip_on_cran()

  # ANCOM-BC2's own warnings (here: few taxa) are passed on.
  expect_warning(
    da <- taxa_da(saliva_stool, stat.by = "Body Site", rank = "Phylum", method = "ancombc2", seed = 1),
    "large number of taxa" )

  expect_true(hasName(da, '.robust'))
  expect_identical(levels(da$`Body Site`), "Stool - Saliva")
  expect_match(da$code, "ANCOMBC::ancombc2(", fixed = TRUE)

  counts <- taxa_matrix(saliva_stool, "Phylum")
  md <- data.frame(
    row.names = saliva_stool$samples,
    grp = droplevels(pull(saliva_stool, "Body Site"))[saliva_stool$samples],
    id  = saliva_stool$samples )
  set.seed(1)
  ref <- suppressWarnings(suppressMessages(ANCOMBC::ancombc2(
    data = counts, meta_data = md, fix_formula = "grp",
    p_adj_method = "fdr", prv_cut = 0.1, verbose = FALSE )))$res
  ref <- ref[!is.na(ref$lfc_grpStool), ]

  x <- da[match(ref$taxon, da$.taxa), ]
  expect_equal(x$.lfc * log(2), ref$lfc_grpStool, tolerance = 1e-10)
  expect_equal(x$.se  * log(2), ref$se_grpStool,  tolerance = 1e-10)
  expect_equal(x$.p.val,  ref$p_grpStool)
  expect_equal(x$.adj.p,  ref$q_grpStool)
  expect_equal(x$.robust, ref$passed_ss_grpStool)

  # The recovered degrees of freedom reproduce the p-values.
  ok <- is.finite(x$.df) & x$.t.ratio != 0
  expect_equal(2 * pt(-abs(x$.t.ratio[ok]), x$.df[ok]), x$.p.val[ok], tolerance = 1e-6)
})


test_that("degrees of freedom recovered from t and p", {
  for (d in c(3, 17.5, 250))
    expect_equal(da_t_dof(2.1, 2 * pt(-2.1, d)), d, tolerance = 1e-6)
  expect_true(is.na(da_t_dof(0, 1)))
  expect_true(is.na(da_t_dof(2, 0)))
  expect_true(is.na(da_t_dof(2, 1)))
  expect_identical(da_t_dof(2, 0.01), Inf) # smaller than the normal tail
})
