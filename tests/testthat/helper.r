
# library(testthat)

hmp50  <- as_rbiom(system.file(package = 'rbiom', 'extdata', 'hmp50.bz2'))
rare50 <- rarefy(hmp50)
min50  <- as_rbiom(list(counts = hmp50$counts))

hmp5 <- hmp50[1:5]

# Make sure sample names can handle any characters, including unicode.
hmp5$samples <- c("~4L<  vl%", "\U1F60A]<?^@D'", "x\\T@;`\U26F9/", '{1"[~/q;$8>', "k_s \U2B50[|L:")

rare5 <- rarefy(hmp5)
min5  <- as_rbiom(list(counts = hmp5$counts))


# An rbiom object with awkward OTU names, and metadata/taxonomy columns of
# every type that the file formats must preserve, for I/O round-trip tests.
# OTU names include every newick metacharacter, so the tree must round-trip
# through quoted labels too.
hostile5 <- local({

  b <- hmp5$clone()

  otus <- b$otus
  otus[1:15] <- c(
    "OTU 1", "a,b", "c:d", "(e)", "f;g", "1.10", "#hash", 'q"uote',
    "back\\slash", "ünïcode", "x_y", "007", "it's", "o__[Sapro]", "'quoted'" )
  b$otus <- otus

  md <- b$metadata
  md$flag  <- c(TRUE, FALSE, NA, TRUE, FALSE)
  md$when  <- as.Date(c("2020-01-02", NA, "2021-12-31", "1999-07-04", "2024-02-29"))
  md$score <- c(1.5, NA, -3, 1e-7, 0)
  md$Sex[2] <- NA
  b$metadata <- md

  tx <- b$taxonomy
  tx$Genus[1]  <- NA
  tx$Family[2] <- NA
  tx$Genus[2]  <- NA
  tx[3, -1]    <- NA
  b$taxonomy <- tx

  b$id      <- 'Hostile "dataset" #1'
  b$comment <- "Unicode — 'quotes' and \"doubles\"."

  b
})



expect_equal_rbiom <- function (a, b) {

  expect_identical(class(a), c('rbiom', 'R6'))
  expect_identical(class(b), c('rbiom', 'R6'))

  expect_identical(a$id,             b$id)
  expect_identical(a$comment,        b$comment)
  expect_identical(a$sequences,      b$sequences)
  expect_identical(a$tree$tip.label, b$tree$tip.label)

  expect_identical(lapply(a$metadata, as.character), lapply(b$metadata, as.character))
  expect_identical(lapply(a$taxonomy, as.character), lapply(b$taxonomy, as.character))

  expect_equal(a$depth,              b$depth)
  expect_equal(as.matrix(a$counts),  as.matrix(b$counts))
  expect_equal(a$tree$edge,          b$tree$edge)
  expect_equal(a$tree$edge.length,   b$tree$edge.length)

}

# Column classes must survive a round trip too.
expect_same_types <- function (a, b) {
  expect_identical(lapply(a$metadata, class), lapply(b$metadata, class))
  expect_identical(names(a$metadata),          names(b$metadata))
}

expect_identical_json <- function (a, b) {
  expect_identical(as.character(jsonlite::toJSON(a)), b)
}


# Convert all factors to character before comparing.
expect_equal_tibble <- function (a, b) {
  expect_s3_class(a, 'tbl_df')
  expect_s3_class(b, 'tbl_df')
  a <- mutate(a, across(where(is.factor), as.character))
  b <- mutate(b, across(where(is.factor), as.character))
  expect_equal(a, b)
}


# Write a tab/comma separated file without readr.
write_delim_file <- function (x, file, sep = "\t") {
  write.table(x, file, sep = sep, row.names = FALSE, qmethod = "double", fileEncoding = "UTF-8")
}
