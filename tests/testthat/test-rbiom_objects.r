test_that("rbiom_objects", {

  expect_silent(hmp5$id        <- hmp5$id)
  expect_silent(hmp5$comment   <- hmp5$comment)
  expect_silent(hmp5$sequences <- hmp5$sequences)
  expect_silent(hmp5$tree      <- hmp5$tree)
  expect_silent(hmp5$taxonomy  <- hmp5$taxonomy)
  expect_silent(hmp5$metadata  <- hmp5$metadata)
  expect_silent(hmp5$counts    <- hmp5$counts)
  expect_silent(hmp5$date      <- hmp5$date)
  expect_silent(hmp5$samples   <- hmp5$samples)
  expect_silent(hmp5$otus      <- hmp5$otus)

  expect_silent(hmp5$fields)
  expect_silent(hmp5$ranks)
  expect_silent(hmp5$n_otus)
  expect_silent(hmp5$n_samples)
  expect_silent(hmp5$n_fields)
  expect_silent(hmp5$n_ranks)
  expect_silent(hmp5$generated_by)
  expect_silent(hmp5$pkg_version)
  expect_silent(hmp5$depth)
  expect_false(hmp5$relative)
  expect_type(hmp5$history, 'list')
  expect_s3_class(hmp5$code, 'rbiom_code')

  skip_on_cran()


  capture.output(expect_message(print(hmp5)), type = 'message')     # rb_print
  capture.output(expect_message(hmp5$print()), type = 'message')

  expect_error(hmp5$initialize())

  expect_error(  hmp5$id <- 1:3)
  expect_warning(hmp5$id <- paste0(collapse = '', rep(letters, 5)))
  expect_silent( hmp5$id <- NA)

  expect_error(  hmp5$comment <- 1:3)
  expect_warning(hmp5$comment <- paste0(collapse = '', rep(letters, 200)))
  expect_silent( hmp5$comment <- NA)

  expect_error(hmp5$sequences <- 1:3)
  expect_error(hmp5$sequences <- c('a' = 'atcg'))

  expect_error(rare50$tree <- hmp5$tree)

  expect_error(min5$counts <- matrix( integer(25),     nrow = 5))
  expect_error(min5$counts <- matrix(!logical(25),     nrow = 5))
  expect_error(min5$counts <- matrix( integer(25) - 1, nrow = 5))
  expect_error(min5$counts <- matrix( integer(25) + 1, nrow = 5))
  expect_error(min5$counts <- matrix( integer(25) + 1, nrow = 5, dimnames = list(NULL, letters[1:5])))

  m <- min5$counts; rownames(m) %<>% paste0(' ', .); expect_warning(min5$counts <- m)
  m <- min5$counts; colnames(m) %<>% paste0(' ', .); expect_warning(min5$counts <- m)

  expect_error(min5$counts <- matrix( integer(25) + 1, nrow = 5, dimnames = list(rep('a', 5), letters[1:5])))
  expect_error(min5$counts <- matrix( integer(25) + 1, nrow = 5, dimnames = list(letters[1:5], rep('a', 5))))

  m <- min5$counts; rownames(m) <- paste0(seq_len(nrow(m))); expect_error(min5$counts <- m)
  m <- min5$counts; colnames(m) <- paste0(seq_len(ncol(m))); expect_error(min5$counts <- m)

  expect_error(min5$date <- 1:2)
  expect_error(min5$date <- 'notadate')
  expect_error(min5$date <- FALSE)
})


test_that("counts validation gives clear errors", {

  dn <- list(c('a', 'b'), c('x', 'y'))

  expect_error(as_rbiom(matrix(c(1, NA, 2, 3),       2, dimnames = dn)), "NA")
  expect_error(as_rbiom(matrix(c("1", "2", "3", "4"), 2, dimnames = dn)), "numeric")
  expect_error(as_rbiom(matrix(c(1, -1, 2, 3),       2, dimnames = dn)), "negative")
  expect_error(as_rbiom(list(a = 1)))

  # data.frames are accepted
  expect_s3_class(as_rbiom(as.data.frame(matrix(1:4, 2, dimnames = dn))), 'rbiom')
})


test_that("single sample / single OTU objects", {

  b <- hmp50[1]
  expect_equal(b$n_samples, 1L)
  expect_equal(nrow(b$metadata), 1L)

  b <- hmp50[1, ]
  expect_equal(b$n_otus, 1L)

  w <- which(as.matrix(hmp50$counts) > 0, arr.ind = TRUE)[1, ]
  b <- hmp50[w[['row']], w[['col']]]
  expect_equal(dim(b$counts), c(1L, 1L))
})


test_that("relative abundance detection", {

  r <- biom_relativize(hmp50)
  expect_true(r$relative)
  expect_null(r$depth)
  expect_false(rarefy(hmp5)$relative)
  # print() writes with cli, which signals each line as a message. Catch them
  # here, before any outer handler (e.g. suppressMessages() around the whole
  # test run) can muffle them.
  out <- capture_messages(print(r))
  expect_true(any(grepl("relative abundances", out, fixed = TRUE)))
})


test_that("depth after rarefaction", {
  expect_equal(rarefy(hmp50, depth = 1000)$depth, 1000)
  expect_null(hmp50$depth)
})


test_that("sample and OTU names are validated and trimmed", {

  b <- hmp5$clone()
  expect_warning(b$samples <- paste0(b$samples, " "), "Trimming")
  expect_identical(b$samples, trimws(b$samples))
  expect_identical(b$metadata$.sample, b$samples)

  expect_error(b$samples <- rep("a", b$n_samples), "Duplicated")
  expect_error(b$samples <- head(b$samples, 2),   "Expected")
  expect_error(b$otus    <- rep(NA_character_, b$n_otus), "NA")
  expect_error(b$samples <- c(" ", b$samples[-1]),     "empty")
})


test_that("clone independence", {
  a <- hmp5$clone()
  b <- a$clone()
  b$metadata$Age <- 0
  expect_false(all(a$metadata$Age == 0))
  expect_equal(length(a$history), length(b$history) - 1L)
})
