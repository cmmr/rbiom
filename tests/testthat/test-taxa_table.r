
replay_data <- function (x) {
  env <- new.env(parent = globalenv())
  eval(parse(text = x$code), envir = env)
  env$data
}


test_that("taxa_table", {

  skip_on_cran()

  expect_silent(taxa_table(hmp5, transform = "percent"))
  expect_silent(taxa_table(hmp5, transform = "rank"))
  expect_silent(taxa_table(hmp5, transform = "log1p"))
  expect_silent(taxa_table(hmp5, taxa = 0.01))
  expect_error(taxa_table(hmp5, taxa = FALSE))
  expect_error(taxa_table(hmp5, taxa = 'doesnotexist'), "doesnotexist")

  expect_identical(attr(taxa_table(hmp5), 'resp_label'), "Unrarefied Counts")
  expect_identical(attr(taxa_table(rare5), 'resp_label'), "Rarefied Counts")
  expect_identical(attr(taxa_table(biom_relativize(hmp5)), 'resp_label'), "Relative Abundance")
})


test_that("taxa_matrix sums match a manual rollup", {

  map <- taxa_map(hmp5, 'Phylum')
  man <- rowsum(as.matrix(hmp5$counts)[names(map), ], as.character(map))
  mtx <- taxa_matrix(hmp5, 'Phylum')
  expect_equal(mtx[rownames(man), ], man, ignore_attr = TRUE)

  sums <- taxa_sums(hmp5, 'Phylum')
  expect_equal(unname(sums[rownames(man)]), unname(rowSums(man)))

  # Top taxa are ranked by mean relative abundance.
  rel <- t(t(man) / colSums(man))
  top <- names(sort(rowMeans(rel), decreasing = TRUE))[1:2]
  expect_setequal(rownames(taxa_matrix(hmp5, 'Phylum', taxa = 2)), top)
})


test_that("log transform uses a pseudocount", {
  raw <- taxa_matrix(hmp5, 'Phylum')
  lg  <- taxa_matrix(hmp5, 'Phylum', transform = 'log')
  expect_true(all(is.finite(lg)))
  pc  <- min(raw[raw > 0]) / 2
  expect_equal(lg, log(raw + pc), ignore_attr = TRUE)
  expect_true(all(lg[raw == 0] < min(lg[raw > 0])))
})


test_that("taxa_matrix ranks taxa when a sample has none at the rank", {

  # S3's only OTU has no Genus, so unc = "drop" leaves it an empty column.
  biom <- as_rbiom(list(
    counts = matrix(
      data     = c(5, 3, 0,  2, 4, 0,  1, 1, 9),
      nrow     = 3,
      byrow    = TRUE,
      dimnames = list(c("OTU1", "OTU2", "OTU3"), c("S1", "S2", "S3")) ),
    taxonomy = data.frame(
      .otu   = c("OTU1", "OTU2", "OTU3"),
      Phylum = c("Firmicutes", "Bacteroidetes", "Firmicutes"),
      Genus  = c("Lactobacillus", "Bacteroides", NA) )))

  # Mean shares over S1 and S2: Lactobacillus 8/14, Bacteroides 6/14.
  expect_identical(
    rownames(taxa_matrix(biom, 'Genus', taxa = 1, unc = 'drop')),
    'Lactobacillus' )

  expect_setequal(
    rownames(taxa_matrix(biom, 'Genus', taxa = 0.1, unc = 'drop')),
    c('Lactobacillus', 'Bacteroides') )

  mtx <- taxa_matrix(biom, 'Genus', taxa = 1, other = TRUE, unc = 'drop')
  expect_identical(nrow(mtx), 2L)
  expect_equal(unname(colSums(mtx)), c(7, 7, 0))

})


test_that("subset_taxa", {

  b <- subset_taxa(hmp5, Phylum == 'Bacteroidetes')
  expect_true(all(pull(b$taxonomy, 'Phylum') == 'Bacteroidetes'))
  expect_true(grepl("subset_taxa(biom, Phylum == \"Bacteroidetes\")", b$code, fixed = TRUE))

  # NA in the condition means drop, as in base::subset().
  b2 <- subset_taxa(hmp50, Genus == "Streptococcus")
  expect_true(all(pull(b2$taxonomy, 'Genus') == "Streptococcus"))

  expect_error(subset_taxa(hmp50, Genus == "NoSuchGenus"), "No OTUs")
  expect_error(subset_taxa(hmp50, Genus), "logical")

  # Free variables are inlined so the history replays.
  phy <- "Firmicutes"
  b3  <- subset_taxa(hmp50[1:5], Phylum == phy)
  expect_true(grepl('Phylum == "Firmicutes"', b3$code, fixed = TRUE))
  env <- new.env(parent = globalenv())
  eval(parse(text = b3$code), envir = env)
  expect_equal_rbiom(env$biom, b3)
})


test_that("taxa provenance replays", {

  skip_on_cran()

  b   <- rarefy(hmp5, depth = 1000, seed = 3)
  tbl <- taxa_table(b, rank = "Genus", taxa = 5, md = "Sex")
  expect_equal(as.data.frame(replay_data(tbl)), as.data.frame(tbl), ignore_attr = TRUE)

  mtx <- taxa_matrix(b, rank = "Phylum", taxa = 3)
  expect_equal(unclass(replay_data(mtx)), unclass(mtx), ignore_attr = TRUE)

  rlang::with_options(rbiom.cache_dir = "FALSE", {
    expect_equal(taxa_matrix(b, rank = "Phylum", taxa = 3), mtx)
  })
})
