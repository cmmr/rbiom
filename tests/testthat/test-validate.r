test_that("validate", {
  
  skip_on_cran()
  
  tree <- expect_silent(write_tree(tree_subset(hmp5$tree, 1:5)))
  expect_null(validate_tree())
  expect_s3_class(tree, 'phylo')
  
  tree <- NA; expect_error(validate_tree())
  
  x <- 1; expect_error(validate_bool('x'))
  
  taxa <- NULL; expect_error(validate_taxa())
  taxa <- NA;   expect_error(validate_taxa())
  taxa <- TRUE; expect_error(validate_taxa())
  taxa <- 1:3;  expect_error(validate_taxa())
  taxa <- -1;   expect_error(validate_taxa())
  
  layers <- NA;  expect_error(validate_layers())
  layers <- 10;  expect_error(validate_layers())
  layers <- 'z'; expect_error(validate_layers(choices = 'point'))
  
  x <- character(0); expect_error(validate_string('x'))
  x <- character(2); expect_error(validate_string('x'))
  x <- NA;           expect_error(validate_string('x'))
  x <- NULL;         expect_error(validate_string('x'))
  x <- TRUE;         expect_null(validate_string('x'))
  x <- FALSE;        expect_null(validate_string('x'))
  
  remove('x')
  expect_error(validate_var_choices('x', dne_ok = FALSE))
  expect_null(validate_var_choices('x', dne_ok = TRUE))
  x <- 1;           expect_error(validate_var_choices('x', letters))
  x <- c('a', 'b'); expect_error(validate_var_choices('x', letters))
  x <- '1';         expect_error(validate_var_choices('x', letters))
  x <- NA;          expect_error(validate_var_choices('x', letters))
  x <- '.all'; expect_null(validate_var_choices('x', mtcars, all_option = '.all', max = Inf))
  
  remove('x')
  expect_error(validate_var_range('x', dne_ok = FALSE))
  expect_null(validate_var_range('x', dne_ok = TRUE))
  x <- NA;   expect_error(validate_var_range('x'))
  x <- NULL; expect_error(validate_var_range('x'))
  x <- 'A';  expect_error(validate_var_range('x'))
  x <- 1:3;  expect_error(validate_var_range('x', n = 1))
  x <- 1.4;  expect_error(validate_var_range('x', int = TRUE))
  x <- 8;    expect_error(validate_var_range('x', range = c(1,5)))
  
  x <- NA;   expect_error(validate_var_length('x', x, max = 1, null_ok = FALSE, na_ok = FALSE))
  x <- NULL; expect_error(validate_var_length('x', x, max = 1, null_ok = FALSE, na_ok = FALSE))
  x <- 1:2;  expect_error(validate_var_length('x', x, max = 1, null_ok = FALSE, na_ok = FALSE))
  
  df <- ToothGrowth
  x <- NULL;          expect_error(validate_df_field('x'))
  x <- 1;             expect_error(validate_df_field('x'))
  x <- c('a', 'b');   expect_error(validate_df_field('x'))
  x <- NA_character_; expect_error(validate_df_field('x'))
  x <- '';            expect_error(validate_df_field('x'))
  x <- 'dose';        expect_error(validate_df_field('x', col_type = 'cat', force = FALSE))
  x <- 'supp';        expect_error(validate_df_field('x', col_type = 'num', force = FALSE))
  df <- data.frame('.sample' = letters)
  x <- '.all'; expect_null(validate_df_field('x', null_ok = TRUE)); expect_null(x)
  
})
