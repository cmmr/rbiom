test_that("s3_methods", {

  skip_on_cran()

  df <- expect_silent(stats_table(adiv_table(rare50), stat.by = 'Sex'))
  expect_silent(df$.n)         # $.rbiom_tbl
  expect_output(print(df))     # tbl_sum.rbiom_tbl
  expect_output(print(df$cmd)) # print.rbiom_code
  expect_output(print(df$code))
  expect_s3_class(df$history, 'rbiom_code')

  expect_silent(as.list(min5))   # as.list.rbiom
  expect_silent(as.matrix(min5)) # as.matrix.rbiom

  expect_silent(pull(min5))             # pull.rbiom
  expect_output(glimpse(min5))          # glimpse.rbiom
  expect_silent(mutate(hmp5, x = 1))    # mutate.rbiom
  expect_silent(rename(hmp5, x = Sex))  # rename.rbiom
  expect_silent(with(hmp5, Age * 2))    # with.rbiom
  expect_silent(within(hmp5, x <- 1))   # within.rbiom
  expect_silent(subset(hmp5, Age > 25)) # subset.rbiom

  expect_silent(hmp5[1:3])      # `[.rbiom`(i)
  expect_silent(hmp5[1:10,1:3]) # `[.rbiom`(i,j)
  expect_silent(hmp5[, 1:3])    # `[.rbiom`(,j)
  expect_silent(hmp5[1:10, ])   # `[.rbiom`(i,)
  expect_silent(na.omit(hmp5))  # na.omit.rbiom

  expect_silent(slice(hmp5, 1:2, 4))         # slice.rbiom
  expect_silent(slice_head(hmp5, n = 2))     # slice_head.rbiom
  expect_silent(slice_tail(hmp5, n = 2))     # slice_tail.rbiom
  expect_silent(slice_min(hmp5, Age, n = 2)) # slice_min.rbiom
  expect_silent(slice_max(hmp5, Age, n = 2)) # slice_max.rbiom
  expect_silent(slice_sample(hmp5, n = 2))   # slice_sample.rbiom

})


test_that("subset() treats NA as FALSE and validates its input", {

  b <- mutate(hmp50, A2 = ifelse(Age > 30, NA, Age))

  s <- subset(b, A2 < 25)
  expect_equal(s$n_samples, sum(b$metadata$A2 < 25, na.rm = TRUE))

  expect_error(subset(b, Age),            "logical")
  expect_error(subset(b, Age > 1000),     "No samples")
  expect_error(subset(b, c(TRUE, FALSE)), "length")
})


test_that("[ validates names and handles NA", {

  expect_error(hmp50["nope"],     "Unknown sample")
  expect_error(hmp50["nope", ],   "Unknown OTU")

  keep <- rep(TRUE, 50); keep[2] <- NA
  expect_equal(hmp50[keep]$n_samples, 49L)

  # clone = FALSE modifies in place
  b <- hmp50$clone()
  hmp50[1:3, clone = FALSE] -> x
  expect_equal(hmp50$n_samples, 3L)
  expect_identical(x, hmp50)
  hmp50 <<- b
})


test_that("transformations are recorded in $history", {

  b <- hmp5$clone()
  n <- length(b$history)

  b <- subset(b, Age > 0)
  b <- mutate(b, x = 1)
  b <- rename(b, y = x)
  b <- within(b, z <- 2)
  b <- slice_head(b, n = 3)
  b <- na.omit(b)
  b <- b[1:2]

  expect_equal(length(b$history), n + 7L)
  expect_true(all(vapply(b$history, is.language, TRUE)))

  # Direct edits are logged as non-reproducible assignments.
  b$id <- "New ID"
  b$metadata <- b$metadata
  code <- b$code
  expect_true(grepl('biom$id <- "New ID"', code, fixed = TRUE))
  expect_true(grepl("not reproducible", code, fixed = TRUE))
})
