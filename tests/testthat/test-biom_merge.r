test_that("biom_merge", {

  actual        <- biom_merge(hmp5[1:2], hmp5[3:5], tree = hmp5$tree)
  actual$counts <- actual$counts[sort(actual$otus), sort(actual$samples)]

  expected        <- hmp5$clone()
  expected$counts <- expected$counts[sort(expected$otus), sort(expected$samples)]

  expect_equal_rbiom(actual, expected)

  skip_on_cran()

  expect_error(biom_merge())
  expect_error(biom_merge(NA))
  expect_error(biom_merge(list()))
  expect_error(biom_merge(list(NA, NA)))

  # Overlapping samples: fine when identical, an error otherwise.
  expect_equal_rbiom(biom_merge(min5, min5), min5)
  x <- min5$clone(); x$counts <- x$counts * 2
  expect_error(biom_merge(min5, x), "different counts")

  expect_warning(biom_merge(min5[1:2], {x <- min5[3:5]; x$otus %<>% paste0('_'); x}))

  expect_equal_rbiom(biom_merge(hmp5), hmp5)
  expect_s3_class(biom_merge(list(min5)), 'rbiom')
  expect_s3_class(biom_merge(system.file("extdata", "hmp50.bz2", package = "rbiom")), 'rbiom')

  # Matrices, lists, and file names are accepted.
  m1 <- as.matrix(hmp5$counts[, 1:2])
  m2 <- hmp5$counts[, 3:5]
  x  <- expect_silent(biom_merge(m1, m2))
  expect_identical(sort(x$samples), sort(hmp5$samples))
  x  <- expect_silent(biom_merge(list(m1, list(counts = m2))))
  expect_identical(sort(x$samples), sort(hmp5$samples))

  # Explicit components.
  x <- biom_merge(hmp5[1:2], hmp5[3:5], metadata = NULL, taxonomy = NULL, sequences = NULL, id = "x", comment = "y")
  expect_identical(x$fields, ".sample")
  expect_identical(x$ranks,  ".otu")
  expect_null(x$sequences)
  expect_identical(x$id, "x")
  expect_identical(x$comment, "y")
})


test_that("biom_merge unions metadata, taxonomy, sequences and trees", {

  skip_on_cran()

  b1 <- hmp5[1:2]
  b2 <- hmp5[3:5]

  # Metadata type conflicts: one summarizing warning, converted to categorical.
  # (Logical and numeric are compatible and stay numeric.)
  b2 <- mutate(b2, Age = as.character(Age), BMI = as.Date(BMI, origin = "1970-01-01"), Obese = BMI > 22)
  b1 <- mutate(b1, Obese = 1)
  expect_warning(x <- biom_merge(b1, b2, tree = NULL), "Age.*BMI")
  expect_true(is.factor(x$metadata$Age))
  expect_true(is.factor(x$metadata$BMI))
  expect_true(is.numeric(x$metadata$Obese))
  b1 <- hmp5[1:2]

  # New fields on one side are NA on the other.
  b2 <- mutate(hmp5[3:5], extra = 1)
  x  <- expect_silent(biom_merge(b1, b2, tree = NULL))
  expect_identical(x$fields, c(hmp5$fields, "extra"))
  expect_identical(sum(is.na(x$metadata$extra)), 2L)

  # Taxonomy conflicts: one warning, most detailed lineage kept.
  b2 <- hmp5[3:5]
  tx <- b2$taxonomy
  shared <- intersect(b1$otus, b2$otus)
  tx$Genus[tx$.otu %in% shared] <- NA
  b2$taxonomy <- tx
  expect_warning(x <- biom_merge(b1, b2, tree = NULL), "taxonomy")
  expect_identical(
    as.character(x$taxonomy$Genus[match(shared, x$taxonomy$.otu)]),
    as.character(b1$taxonomy$Genus[match(shared, b1$taxonomy$.otu)]) )

  # Sequences in only one input: kept, NA elsewhere, one warning.
  b2 <- hmp5[3:5]; b2$sequences <- NULL
  expect_warning(x <- biom_merge(b1, b2, tree = NULL), "sequence")
  expect_identical(x$sequences[b1$otus], b1$sequences)
  expect_true(anyNA(x$sequences))

  # Trees: subsets of one dataset keep its tree (the stored, unpruned tree
  # covers every merged OTU), pruned to the merged OTUs.
  x <- expect_silent(biom_merge(hmp50[1:4], hmp50[5:8]))
  expect_s3_class(x$tree, 'phylo')
  expect_setequal(x$tree$tip.label, x$otus)
  expect_identical(x$tree, tree_subset(hmp50$tree, x$otus))
  dm <- bdiv_distmat(x, 'unweighted_unifrac')
  expect_s3_class(dm, 'dist')
  expect_identical(attr(dm, 'Size'), 8L)
  expect_true(all(is.finite(dm)))

  # Pruning is compositional, so the merged tree doesn't depend on whether
  # an input's tree had already been pruned (as long as it still covers).
  s1 <- hmp50$otus[1:200]; s2 <- s1[seq(1, 200, 3)]
  a  <- tree_subset(tree_subset(hmp50$tree, s1), s2)
  b  <- tree_subset(hmp50$tree, s2)
  expect_identical(a[c('edge', 'tip.label', 'node.label')], b[c('edge', 'tip.label', 'node.label')])
  expect_equal(a$edge.length, b$edge.length)   # summed lengths; floating point

  # Trees already pruned to each input's own OTUs: none covers the union,
  # so the merged object has no tree and a message says so.
  b1 <- hmp5[1:2]; b1$tree <- b1$tree
  b2 <- hmp5[3:5]; b2$tree <- b2$tree
  expect_false(all(b2$otus %in% b1$tree$tip.label))
  expect_message(x <- biom_merge(b1, b2), "tree")
  expect_null(x$tree)

  # Same OTUs but different samples: kept.
  b1 <- hmp5$clone(); b1$samples <- paste0("a", seq_len(5))
  b2 <- hmp5$clone(); b2$samples <- paste0("b", seq_len(5))
  x  <- expect_silent(biom_merge(b1, b2))
  expect_identical(x$tree$tip.label, hmp5$tree$tip.label)
  expect_identical(x$n_samples, 10L)
})


test_that("biom_merge provenance", {

  skip_on_cran()

  f  <- system.file("extdata", "hmp50.bz2", package = "rbiom")
  b1 <- subset(read_biom(f), Age > 30)
  b2 <- subset(read_biom(f), Age <= 30)
  m  <- suppressMessages(biom_merge(b1, b2))

  expect_length(m$history, 1L)
  expect_identical(m$history[[1]], quote(biom_merge(biom1, biom2)))

  code <- m$code
  expect_true(grepl("biom1 <- read_biom(", code, fixed = TRUE))
  expect_true(grepl("biom2 <- subset(biom2, Age <= 30)", code, fixed = TRUE))
  expect_true(grepl("biom <- biom_merge(biom1, biom2)", code, fixed = TRUE))
  expect_false(grepl("not reproducible", code, fixed = TRUE))

  env <- new.env(parent = globalenv())
  suppressMessages(eval(parse(text = code), envir = env))
  expect_identical(env$biom$hash, m$hash)

  # Extra arguments are recorded; in-memory values are flagged.
  m <- suppressMessages(biom_merge(b1, b2, id = "merged", tree = hmp50$tree))
  expect_identical(m$history[[1]][["id"]], "merged")
  expect_true(grepl("not reproducible", m$code, fixed = TRUE))
})
