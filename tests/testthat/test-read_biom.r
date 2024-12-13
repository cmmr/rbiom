test_that("read_biom", {
  
  skip_on_cran()
  
  f <- expect_silent(write_biom(hmp5, tempfile()))
  expect_silent(read_biom_internal(src = f, tree = TRUE))
  expect_silent(read_biom_internal(src = f, tree = FALSE))
  unlink(f)
  
  f <- expect_silent(write_biom(hmp5, tempfile(), 'hdf5'))
  expect_silent(read_biom_internal(src = f))
  unlink(f)
  
  f <- expect_silent(write_biom(hmp5, tempfile(fileext = '.gz'), 'tab'))
  expect_silent(read_biom_internal(src = f))
  expect_error(read_biom_internal(src = f, tree = TRUE))
  unlink(f)
  
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
  
  expect_null(parse_json_tree(list(phylogeny = NULL),         tree_mode = 'auto'))
  expect_null(parse_json_tree(list(phylogeny = NA),           tree_mode = 'auto'))
  expect_null(parse_json_tree(list(phylogeny = FALSE),        tree_mode = 'auto'))
  expect_null(parse_json_tree(list(phylogeny = character(0)), tree_mode = 'auto'))
  expect_null(parse_json_tree(list(phylogeny = c('a', 'b')),  tree_mode = 'auto'))
  expect_null(parse_json_tree(list(phylogeny = ''),           tree_mode = 'auto'))
  
  tree <- expect_silent(write_tree(tree_subset(hmp5$tree, 1:5)))
  hdf5 <- list(observation = list('group-metadata' = list()))
  
  expect_error(parse_json_tree(list(),  tree_mode = 'not a tree'))
  expect_error(parse_json_tree(list(rows = list(list(id='t1'))), tree_mode = tree))
  
  expect_null(parse_hdf5_tree(list(),  tree_mode = FALSE))
  expect_null(parse_hdf5_tree(hdf5,    tree_mode = 'auto'))
  expect_error(parse_hdf5_tree(hdf5,   tree_mode = TRUE))
  expect_error(parse_hdf5_tree(list(), tree_mode = tree))
  expect_error(parse_hdf5_tree(list(), tree_mode = 'not a tree'))
  expect_error(parse_hdf5_tree(list(rows = list(list(id='t1'))), tree_mode = tree))
})
