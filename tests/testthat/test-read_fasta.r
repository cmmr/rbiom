test_that("read_fasta", {
  
  skip_on_cran()
  
  fp <- tempfile(fileext = '.fna')
  expect_error(read_fasta(fp))
  
  expect_silent(write_fasta(c(a = 'atcg', a = 'atcg'), fp))
  expect_error(read_fasta(fp))
  
  expect_silent(write_fasta(c(a = 'atcg', b = 'atcg'), fp))
  expect_error(read_fasta(fp, ids = c('b', 'c')))
  expect_silent(read_fasta(fp, ids = c('b', 'a')))
})
