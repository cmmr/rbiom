test_that("read_biom", {
  
  skip_on_cran()
  
  f <- expect_silent(write_biom(hmp5, tempfile(), 'json'))
  g <- expect_silent(write_biom(min5, tempfile(), 'json'))
  expect_silent(read_biom(src = f))
  expect_silent(read_biom(src = g))
  unlink(c(f, g))
  
  f <- expect_silent(write_biom(hmp5, tempfile(), 'hdf5'))
  g <- expect_silent(write_biom(min5, tempfile(), 'hdf5'))
  expect_silent(read_biom(src = f))
  expect_silent(read_biom(src = g))
  unlink(c(f, g))
  
  f <- expect_silent(write_biom(hmp5, tempfile(fileext = '.gz'), 'tab'))
  g <- expect_silent(write_biom(min5, tempfile(fileext = '.gz'), 'tab'))
  expect_silent(read_biom(src = f))
  expect_silent(read_biom(src = g))
  unlink(c(f, g))
  
  expect_error(read_biom_hdf5(tempfile()))
  expect_silent(parse_hdf5_taxonomy(list(observation = list(ids = c('a', 'b', 'c')))))
  expect_error(parse_json_metadata(list(columns = list(list(metadata = c('a' = 1, 'a' = 1))))))
  expect_silent(parse_hdf5_metadata(list(observation = list(ids = c('a', 'b', 'c')))))
  
  json <- list(matrix_type = 'dense', rows = list(list(id='t1')), columns = list(list(id='s1')), data = list(1))
  expect_silent(parse_json_counts(json))
  expect_error(parse_json_counts(list(matrix_type = 'invalid')))
  expect_error(parse_json_counts(list(matrix_type = 'dense', data = list())))
  
  
  json <- list(generated_by = 'MicrobiomeDB', rows = list(list(id='Archaea;Euryarchaeota', metadata = list())))
  expect_silent(parse_json_taxonomy(json))
  json <- list(rows = list(list(id='Unc01pdq', metadata = list(list("Bacteria", "__Fusobacteriota")))))
  expect_silent(parse_json_taxonomy(json))
  expect_silent(parse_json_taxonomy(list(rows = list(list(id='t1')))))
  json <- list(rows = list(list(metadata=list(taxonomy='t1')), list(metadata=list(taxonomy='t2;t3'))))
  expect_silent(parse_json_taxonomy(json))
  
  expect_null(parse_json_sequences(list(rows = list(list(metadata=list())))))
  expect_null(parse_hdf5_sequences(list(observation = list(list(metadata=list())))))
})
