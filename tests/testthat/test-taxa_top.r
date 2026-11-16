
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

replay_data <- function (x) {
  env <- new.env(parent = globalenv())
  eval(parse(text = x$code), envir = env)
  env$data
}


test_that("taxa_top: all samples together", {

  biom <- tiny_biom()

  tbl <- taxa_top(biom, n = 2)

  expect_s3_class(tbl, 'rbiom_tbl')
  expect_named(tbl, c('.rank', '.top', '.taxa', '.mean', '.prevalence'))
  expect_identical(as.character(tbl$.rank), c("Genus", "Genus"))
  expect_identical(tbl$.top, 1:2)
  expect_identical(as.character(tbl$.taxa), c("G3", "G1"))
  expect_equal(tbl$.mean, c(224, 160) / 600)
  expect_equal(tbl$.prevalence, c(1, 4/6))

  # Counts rather than fractions.
  tbl <- taxa_top(biom, n = 2, transform = "none")
  expect_equal(tbl$.mean, c(224, 160) / 6)

  # Other transformations change the units, not the taxa listed.
  tbl <- taxa_top(biom, n = 3, transform = "log1p")
  expect_identical(as.character(tbl$.taxa), c("G3", "G1", "G4"))
  expect_equal(tbl$.mean[[1]], mean(log1p(c(20, 20, 39, 45, 50, 50))))

  # More than there are: all the taxa.
  expect_identical(nrow(taxa_top(biom, n = 99)), 5L)

  # The same taxa that taxa_table() and taxa_matrix() select.
  expect_setequal(
    as.character(taxa_top(rare50, rank = "Genus", n = 6)$.taxa),
    rownames(taxa_matrix(rare50, rank = "Genus", taxa = 6)) )
})


test_that("taxa_top: per group", {

  biom <- tiny_biom()

  tbl <- taxa_top(biom, n = 5, by = "Site")

  expect_named(tbl, c('.rank', 'Site', '.top', '.taxa', '.mean', '.prevalence'))

  # Groups in factor level order; absent taxa aren't listed.
  expect_identical(levels(tbl$Site), c("Gut", "Arm"))
  expect_identical(as.character(tbl$Site), rep(c("Gut", "Arm"), each = 4))
  expect_identical(tbl$.top, c(1:4, 1:4))
  expect_identical(as.character(tbl$.taxa), c("G3", "G4", "G5", "G1",   "G1", "G3", "G2", "G4"))
  expect_equal(tbl$.mean, c(145, 125, 20, 10,   150, 79, 70, 1) / 300)
  expect_equal(tbl$.prevalence, c(3, 3, 1, 1,   3, 3, 2, 1) / 3)

  # The top taxa differ between groups.
  tbl <- taxa_top(biom, n = 1, by = "Site")
  expect_identical(as.character(tbl$.taxa), c("G3", "G1"))

  # Abbreviated field names, as elsewhere.
  expect_identical(names(taxa_top(biom, n = 1, by = "si"))[[2]], "Site")

  # Two fields: one group per combination that occurs.
  tbl <- taxa_top(biom, n = 1, by = c("Site", "Sex"))
  expect_named(tbl, c('.rank', 'Site', 'Sex', '.top', '.taxa', '.mean', '.prevalence'))
  expect_identical(as.character(tbl$Site), c("Gut", "Gut", "Arm", "Arm"))
  expect_identical(as.character(tbl$Sex),  c("F", "M", "F", "M"))
  expect_identical(as.character(tbl$.taxa), c("G3", "G3", "G1", "G1"))
  expect_equal(tbl$.mean, c(50, 95 / 2, 110 / 2, 40) / 100)

  # Samples with NA in `by` are left out.
  b  <- tiny_biom()
  md <- b$metadata
  md$Site[[6]] <- NA
  b$metadata <- md
  tbl <- taxa_top(b, n = 5, by = "Site")
  expect_false("G5" %in% tbl$.taxa)
  expect_false(anyNA(tbl$Site))
})


test_that("taxa_top: per sample, with ties broken by name", {

  biom <- tiny_biom()

  tbl <- taxa_top(biom, n = 2, by = ".sample")

  expect_named(tbl, c('.rank', '.sample', '.top', '.taxa', '.abundance'))
  expect_identical(as.character(tbl$.sample), rep(paste0("S", 1:6), each = 2))
  expect_identical(
    as.character(tbl$.taxa),
    c("G1", "G2",   "G1", "G2",   "G1", "G3",   "G3", "G4",   "G3", "G4",   "G3", "G4") )
  expect_equal(tbl$.abundance, c(50, 30,   40, 40,   60, 39,   45, 45,   50, 50,   50, 30) / 100)

  # Other fields come along as columns.
  tbl <- taxa_top(biom, n = 1, by = c("Site", ".sample"))
  expect_named(tbl, c('.rank', 'Site', '.sample', '.top', '.taxa', '.abundance'))
  expect_identical(as.character(tbl$Site), rep(c("Arm", "Gut"), each = 3))

  # No randomness: the same answer every time, without a seed.
  rlang::with_options(rbiom.cache_dir = "FALSE", {
    expect_identical(taxa_top(biom, n = 2, by = ".sample"), taxa_top(biom, n = 2, by = ".sample"))
  })
})


test_that("taxa_top: several ranks, lineage, and unclassified taxa", {

  biom <- tiny_biom()

  tbl <- taxa_top(biom, rank = c("Phylum", "Genus"), n = 1, by = "Site")
  expect_identical(levels(tbl$.rank), c("Phylum", "Genus"))
  expect_identical(as.character(tbl$.rank), c("Phylum", "Phylum", "Genus", "Genus"))
  expect_identical(as.character(tbl$.taxa), c("P2", "P1", "G3", "G1"))
  expect_equal(tbl$.mean, c(290, 220, 145, 150) / 300)

  # Rank numbers.
  expect_identical(as.character(taxa_top(biom, rank = 1, n = 1)$.rank), "Phylum")

  tbl <- taxa_top(biom, n = 1, lineage = TRUE)
  expect_identical(as.character(tbl$.taxa), "P2; G3")

  expect_no_error(taxa_top(rare50, unc = "drop", by = "Body Site"))
  expect_no_error(taxa_top(rare50, rank = "otu", n = 3))
})


test_that("taxa_top: errors", {

  expect_error(taxa_top(hmp50, by = "Age"), "not categorical")
  expect_error(taxa_top(hmp50, by = "nosuchfield"), "nosuchfield")
  expect_error(taxa_top(hmp50, n = 0), "between 1 and Inf")
  expect_error(taxa_top(hmp50, n = 1.5), "whole numbers")
  expect_error(taxa_top(hmp50, n = c(1, 2)), "length 1")
  expect_error(taxa_top(hmp50, rank = "nosuchrank"), "nosuchrank")
  expect_error(taxa_top(hmp50, transform = "zzzz"), "Invalid `transform` option")
})


test_that("taxa_top: schema and provenance", {

  x <- fn_schema("taxa_top")
  expect_identical(x$fn, "taxa_top")
  expect_identical(x$returns$type, "table")
  expect_identical(x$counts, "rarefied")
  expect_identical(vapply(x$params, `[[`, "", 'name'), names(formals(taxa_top)))

  expect_identical(
    fn_check("taxa_top", list(rank = 2, by = "body", transform = "perc"), hmp50),
    list(rank = "Phylum", by = "Body Site", transform = "percent") )
  expect_identical(fn_check("taxa_top", list(by = ".sample"), fn_check_data(hmp50)), list(by = ".sample"))
  expect_error(fn_check("taxa_top", list(by = "Age"), fn_check_data(hmp50)), "not categorical")

  skip_on_cran()

  fresh <- as_rbiom(system.file("extdata", "hmp50.bz2", package = "rbiom"))
  b     <- rarefy(fresh, depth = 1000, seed = 3)

  tbl <- taxa_top(b, rank = "Genus", n = 3, by = "Body Site")
  expect_identical(as.character(tbl$cmd), 'taxa_top(biom, "Genus", 3, "Body Site")')
  expect_false(grepl("not reproducible", tbl$code, fixed = TRUE))
  expect_equal(as.data.frame(replay_data(tbl)), as.data.frame(tbl), ignore_attr = TRUE)

  tbl <- taxa_top(b, n = 1, by = ".sample", transform = "none")
  expect_equal(as.data.frame(replay_data(tbl)), as.data.frame(tbl), ignore_attr = TRUE)

  # Cached and uncached results agree.
  rlang::with_options(rbiom.cache_dir = "FALSE", {
    expect_equal(taxa_top(b, n = 1, by = ".sample", transform = "none"), tbl)
  })
})
