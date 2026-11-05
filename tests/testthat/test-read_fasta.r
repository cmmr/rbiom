test_that("read_fasta", {

  skip_on_cran()

  fp <- tempfile(fileext = '.fna')
  expect_error(read_fasta(fp))

  expect_silent(write_fasta(c(a = 'atcg', a = 'atcg'), fp))
  expect_error(read_fasta(fp))
  unlink(fp)

  expect_silent(write_fasta(c(a = 'atcg', b = 'atcg'), fp))
  expect_error(read_fasta(fp, ids = c('b', 'c')))
  expect_identical(read_fasta(fp, ids = c('b', 'a')), c(b = 'atcg', a = 'atcg'))
  unlink(fp)

  # A named character vector: names stop at the first whitespace,
  # multi-line sequences are joined, and whitespace is dropped.
  txt <- ">OTU1 first sequence\nATCG\n GGCC \r\n\n>OTU2\tdesc\nTTAA\n>OTU3\n"
  x   <- read_fasta(txt)
  expect_type(x, "character")
  expect_identical(x, c(OTU1 = "ATCGGGCC", OTU2 = "TTAA", OTU3 = ""))

  expect_error(read_fasta("ATCG\n>a\nAT"), "header")
  expect_error(read_fasta("no headers here\n"), "header")

  # Assignable to an rbiom object.
  b <- hmp5$clone()
  f <- write_fasta(hmp5, tempfile(fileext = ".fasta"))
  expect_silent(b$sequences <- read_fasta(f))
  expect_identical(b$sequences, hmp5$sequences)
  unlink(f)

  # Compressed files
  for (ext in c(".fna.gz", ".fna.bz2")) {
    f <- write_fasta(hmp5, tempfile(fileext = ext))
    expect_identical(read_fasta(f), hmp5$sequences)
    unlink(f)
  }

  # O(n): 20,000 sequences in well under a few seconds.
  n    <- 20000
  seqs <- setNames(rep("ACGTACGTACGTACGTACGT", n), paste0("seq", seq_len(n)))
  f    <- write_fasta(seqs, tempfile(fileext = ".fna"))
  t0   <- Sys.time()
  x    <- read_fasta(f)
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 3)
  expect_identical(x, seqs)
  unlink(f)
})
