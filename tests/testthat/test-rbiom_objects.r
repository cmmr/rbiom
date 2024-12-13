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
