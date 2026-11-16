
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

# The taxa in each set, as a named list.
set_taxa <- function (tbl) {
  lapply(split(as.character(tbl$.taxa), tbl$.set), sort)
}


test_that("taxa_shared: any observation", {

  biom <- tiny_biom()
  tbl  <- taxa_shared(biom, "Site")

  expect_s3_class(tbl, 'rbiom_tbl')

  # Group columns in factor level order, not alphabetical.
  expect_named(tbl, c('.taxa', 'Gut', 'Arm', '.n_groups', '.set', '.prev.Gut', '.prev.Arm'))

  # Largest set first, then by name.
  expect_identical(as.character(tbl$.taxa), c("G1", "G3", "G4", "G5", "G2"))
  expect_identical(tbl$Gut, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(tbl$Arm, c(TRUE, TRUE, TRUE, FALSE, TRUE))
  expect_identical(tbl$.n_groups, c(2L, 2L, 2L, 1L, 1L))
  expect_identical(as.character(tbl$.set), c("Gut & Arm", "Gut & Arm", "Gut & Arm", "Gut", "Arm"))
  expect_identical(levels(tbl$.set), c("Gut & Arm", "Gut", "Arm"))
  expect_equal(tbl$.prev.Gut, c(1, 3, 3, 1, 0) / 3)
  expect_equal(tbl$.prev.Arm, c(3, 3, 1, 0, 2) / 3)

  # Set sizes, shared with the plot.
  sets <- attr(tbl, 'sets')
  expect_s3_class(sets, 'tbl_df')
  expect_named(sets, c('.set', '.n', '.n_groups', 'Gut', 'Arm'))
  expect_identical(as.character(sets$.set), c("Gut & Arm", "Gut", "Arm"))
  expect_identical(sets$.n, c(3L, 1L, 1L))
  expect_identical(sets$.n_groups, c(2L, 1L, 1L))
  expect_identical(sets$Gut, c(TRUE, TRUE, FALSE))
  expect_identical(sets$Arm, c(TRUE, FALSE, TRUE))
  expect_identical(sum(sets$.n), nrow(tbl))

  # Abbreviations and rank numbers.
  expect_equal(taxa_shared(biom, "si", rank = -1), tbl, ignore_attr = TRUE)

  # Another rank.
  tbl <- taxa_shared(biom, "Site", rank = "Phylum")
  expect_identical(as.character(tbl$.taxa), c("P1", "P2"))
  expect_identical(as.character(tbl$.set),  c("Gut & Arm", "Gut & Arm"))
})


test_that("taxa_shared: prevalence threshold", {

  biom <- tiny_biom()

  # At least half of a group's samples.
  tbl <- taxa_shared(biom, "Site", prev = 0.5)
  expect_identical(set_taxa(tbl), list('Arm' = c("G1", "G2"), 'Gut & Arm' = "G3", 'Gut' = "G4"))
  expect_identical(levels(tbl$.set), c("Arm", "Gut & Arm", "Gut"))
  expect_identical(attr(tbl, 'sets')$.n, c(2L, 1L, 1L))

  # The threshold is inclusive: 2 of 3 samples meets prev = 2/3.
  tbl <- taxa_shared(biom, "Site", prev = 2/3)
  expect_true(tbl$Arm[tbl$.taxa == "G2"])
  tbl <- taxa_shared(biom, "Site", prev = 0.67)
  expect_false("G2" %in% tbl$.taxa)

  # prev = 1: every sample of the group.
  tbl <- taxa_shared(biom, "Site", prev = 1)
  expect_identical(as.character(tbl$.taxa), c("G3", "G4", "G1"))
  expect_identical(as.character(tbl$.set),  c("Gut & Arm", "Gut", "Arm"))
  expect_identical(attr(tbl, 'sets')$.n, c(1L, 1L, 1L))

  # Prevalence is reported for groups the taxon doesn't belong to, too.
  expect_equal(tbl$.prev.Gut[tbl$.taxa == "G1"], 1/3)

  # prev = 0 still needs one sample: taxa in no group aren't listed.
  expect_true(all(taxa_shared(biom, "Site")$.n_groups >= 1))
})


test_that("taxa_shared: detection threshold", {

  biom <- tiny_biom()

  # G4 is 1% of S3, its only Arm sample.
  expect_true(taxa_shared(biom, "Site")$Arm[[3]])

  tbl <- taxa_shared(biom, "Site", detect = 0.05)
  expect_identical(
    set_taxa(tbl),
    list('Gut & Arm' = c("G1", "G3"), 'Gut' = c("G4", "G5"), 'Arm' = "G2") )
  expect_equal(tbl$.prev.Arm[tbl$.taxa == "G4"], 0)

  # Strictly greater than `detect`: G1 is exactly 10% of S4.
  tbl <- taxa_shared(biom, "Site", detect = 0.099)
  expect_true(tbl$Gut[tbl$.taxa == "G1"])
  tbl <- taxa_shared(biom, "Site", detect = 0.10)
  expect_false(tbl$Gut[tbl$.taxa == "G1"])

  # Both thresholds together.
  tbl <- taxa_shared(biom, "Site", detect = 0.25, prev = 1)
  expect_identical(set_taxa(tbl), list('Gut' = c("G3", "G4"), 'Arm' = "G1"))

  # Nothing passes: an empty table, not an error.
  tbl <- taxa_shared(biom, "Site", detect = 0.99)
  expect_identical(nrow(tbl), 0L)
  expect_identical(nrow(attr(tbl, 'sets')), 0L)
})


test_that("taxa_shared: one group of all samples (core microbiome)", {

  biom <- tiny_biom()

  tbl <- taxa_shared(biom)
  expect_named(tbl, c('.taxa', 'All', '.n_groups', '.set', '.prev.All'))
  expect_identical(as.character(tbl$.taxa), paste0("G", 1:5))
  expect_equal(tbl$.prev.All, c(4, 2, 6, 4, 1) / 6)

  tbl <- taxa_shared(biom, prev = 0.6)
  expect_identical(as.character(tbl$.taxa), c("G1", "G3", "G4"))
  expect_true(all(tbl$All))

  expect_identical(as.character(taxa_shared(biom, stat.by = NULL, prev = 1)$.taxa), "G3")
})


test_that("taxa_shared: groups", {

  biom <- tiny_biom()

  # Samples with NA are left out, and empty levels dropped.
  md <- biom$metadata
  md$Site[[6]] <- NA
  biom$metadata <- md

  tbl <- taxa_shared(biom, "Site")
  expect_false("G5" %in% tbl$.taxa)
  expect_equal(tbl$.prev.Gut[tbl$.taxa == "G1"], 1/2)

  # More than two groups; every combination is labelled in level order.
  tbl <- taxa_shared(rare50, "Body Site")
  lvls <- levels(rare50$metadata[["Body Site"]])
  expect_identical(names(tbl)[2:6], lvls)
  expect_true("Anterior nares & Buccal mucosa & Mid vagina & Saliva & Stool" %in% levels(tbl$.set))
  expect_identical(tbl$.n_groups, unname(as.integer(rowSums(as.matrix(tbl[, lvls])))))

  sets <- attr(tbl, 'sets')
  expect_false(is.unsorted(rev(sets$.n)))
  expect_identical(as.integer(table(tbl$.set)[as.character(sets$.set)]), sets$.n)

  # Agrees with a direct computation.
  mtx <- taxa_matrix(rare50, "Genus")
  for (lvl in lvls) {
    ids <- rare50$samples[pull(rare50, "Body Site") == lvl]
    expect_setequal(
      as.character(tbl$.taxa[tbl[[lvl]]]),
      rownames(mtx)[rowSums(mtx[, ids, drop = FALSE] > 0) > 0] )
  }

  # A level that is also a column name can't be used.
  b <- tiny_biom()
  md <- b$metadata
  md$Site <- factor(c(".set", ".set", ".set", "Gut", "Gut", "Gut"))
  b$metadata <- md
  expect_error(taxa_shared(b, "Site"), "column name")
})


test_that("taxa_shared: errors", {

  expect_error(taxa_shared(hmp50, "Age"), "not categorical")
  expect_error(taxa_shared(hmp50, "nosuchfield"), "nosuchfield")
  expect_error(taxa_shared(hmp50, c("Sex", "Body Site")), "length 2")
  expect_error(taxa_shared(hmp50, "Sex", rank = c("Phylum", "Genus")), "length 2")
  expect_error(taxa_shared(hmp50, "Sex", prev = 1.5), "between 0 and 1")
  expect_error(taxa_shared(hmp50, "Sex", prev = -0.1), "between 0 and 1")
  expect_error(taxa_shared(hmp50, "Sex", detect = 2), "between 0 and 1")
  expect_error(taxa_shared(hmp50, "Sex", detect = "a"), "numeric")
  expect_error(taxa_shared(hmp50, "Sex", unc = "zzzz"), "Invalid `unc` option")
})


test_that("taxa_shared: schema and provenance", {

  x <- fn_schema("taxa_shared")
  expect_identical(x$fn, "taxa_shared")
  expect_identical(x$returns$type, "table")
  expect_identical(x$counts, "rarefied")
  expect_identical(vapply(x$params, `[[`, "", 'name'), names(formals(taxa_shared)))

  expect_identical(
    fn_check("taxa_shared", list(stat.by = "body", rank = 2, prev = 0.5), fn_check_data(hmp50)),
    list(stat.by = "Body Site", rank = "Phylum", prev = 0.5) )
  expect_error(fn_check("taxa_shared", list(stat.by = "Age"), fn_check_data(hmp50)), "not categorical")
  expect_error(fn_check("taxa_shared", list(prev = 2), fn_check_data(hmp50)), "between 0 and 1")

  skip_on_cran()

  fresh <- as_rbiom(system.file("extdata", "hmp50.bz2", package = "rbiom"))
  b     <- rarefy(fresh, depth = 1000, seed = 3)

  tbl <- taxa_shared(b, "Body Site", rank = "Genus", detect = 0.001, prev = 0.5)
  expect_identical(as.character(tbl$cmd), 'taxa_shared(biom, "Body Site", "Genus", 0.001, 0.5)')
  expect_false(grepl("not reproducible", tbl$code, fixed = TRUE))

  r <- replay_data(tbl)
  expect_equal(as.data.frame(r), as.data.frame(tbl), ignore_attr = TRUE)
  expect_equal(as.data.frame(attr(r, 'sets')), as.data.frame(attr(tbl, 'sets')))

  tbl <- taxa_shared(b, prev = 0.9)
  expect_equal(as.data.frame(replay_data(tbl)), as.data.frame(tbl), ignore_attr = TRUE)

  # Cached and uncached results agree.
  rlang::with_options(rbiom.cache_dir = "FALSE", {
    expect_equal(taxa_shared(b, prev = 0.9), tbl)
  })
})
