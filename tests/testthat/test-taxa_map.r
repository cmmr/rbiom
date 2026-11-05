test_that("taxa_map", {

  skip_on_cran()

  expect_silent(taxa_map(rare5, unc = 'grouped'))
  expect_silent(taxa_map(rare5, unc = 'drop'))
  expect_silent(taxa_map(rare5, rank = 'Genus', taxa = 4, other = TRUE))
  expect_silent(taxa_map(rare5, rank = 'Genus', taxa = 4, other = FALSE))

  tbl <- taxa_map(rare5)
  expect_s3_class(tbl, 'rbiom_tbl')
  expect_identical(names(tbl)[ncol(tbl)], '.otu')
})


test_that("taxa_map top taxa match taxa_matrix", {

  # The OTU "rank" used to be converted to the wrong column index.
  m <- taxa_map(hmp50, '.otu', taxa = 5, other = TRUE)
  expect_identical(length(m), hmp50$n_otus)
  expect_identical(nlevels(m), 6L)
  expect_setequal(setdiff(levels(m), "Other"), rownames(taxa_matrix(hmp50, '.otu', taxa = 5)))
  expect_true(sum(m == "Other") == hmp50$n_otus - 5)

  m <- taxa_map(hmp50, 'Genus', taxa = 4, other = FALSE)
  expect_setequal(levels(m), rownames(taxa_matrix(hmp50, 'Genus', taxa = 4)))
  expect_true(length(m) < hmp50$n_otus)

  m <- taxa_map(hmp50, 'Genus', taxa = 0.05)
  expect_setequal(levels(m), rownames(taxa_matrix(hmp50, 'Genus', taxa = 0.05)))

  m <- taxa_map(hmp50, 'Genus', taxa = c('Streptococcus', 'Bacteroides'), other = 'Rest')
  expect_setequal(levels(m), c('Streptococcus', 'Bacteroides', 'Rest'))
})


test_that("taxa_map provenance", {
  m   <- taxa_map(rare5, rank = 'Genus', taxa = 4, other = TRUE)
  env <- new.env(parent = globalenv())
  eval(parse(text = m$code), envir = env)
  expect_equal(as.character(env$data), as.character(m))
  expect_equal(names(env$data), names(m))
})
