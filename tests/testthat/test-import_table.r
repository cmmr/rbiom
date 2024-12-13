test_that("import_table", {
  
  file_csv  <- tempfile(fileext = '.csv')
  file_tsv  <- tempfile(fileext = '.tsv')
  on.exit(unlink(c(file_csv, file_tsv)))
  
  readr::write_delim(hmp5$metadata, file_csv, delim = ",")
  expect_equal_tibble(import_metadata(file_csv, sids = hmp5$samples), hmp5$metadata)
  
  readr::write_delim(hmp5$taxonomy, file_tsv, delim = "\t")
  expect_equal_tibble(import_taxonomy(file_tsv, otus = hmp5$otus), hmp5$taxonomy)
  
  skip_on_cran()
  
  
  file_xlsx <- tempfile(fileext = '.xlsx')
  file_json <- tempfile(fileext = '.json')
  on.exit(unlink(c(file_xlsx, file_json)), add = TRUE)
  
  openxlsx::write.xlsx(hmp5$metadata, file_xlsx, overwrite = TRUE)
  expect_equal_tibble(import_metadata(file_xlsx, sids = hmp5$samples), hmp5$metadata)
  
  jsonlite::write_json(hmp5$metadata, file_json)
  expect_equal_tibble(import_metadata(file_json, sids = hmp5$samples), hmp5$metadata)
  
  # as literal strings of data instead of file names
  file_csv  <- readChar(file_csv,  file.size(file_csv))
  file_tsv  <- readChar(file_tsv,  file.size(file_tsv))
  file_json <- readChar(file_json, file.size(file_json))
  expect_equal_tibble(import_metadata(file_csv,  sids = hmp5$samples), hmp5$metadata)
  expect_equal_tibble(import_taxonomy(file_tsv,  otus = hmp5$otus),    hmp5$taxonomy)
  expect_equal_tibble(import_metadata(file_json, sids = hmp5$samples), hmp5$metadata)
  
  # Duplicated column names.
  expect_error(import_table(sub("BMI", "Age", file_csv, fixed = TRUE)))
  
  # Return tibble with just sample or OTU names.
  expect_equal_tibble(import_metadata(NULL, sids = hmp5$samples), hmp5$metadata[,1])
  expect_equal_tibble(import_taxonomy(NULL, otus = hmp5$otus),    hmp5$taxonomy[,1])
  
  # Sample / OTU names under a non-standard column name / position.
  x <- relocate(rename(hmp5$metadata, sample = .sample), sample, .after = BMI)
  y <- relocate(rename(hmp5$taxonomy, otu    = .otu),    otu,    .after = Kingdom)
  expect_equal_tibble(import_metadata(x, sids = hmp5$samples), hmp5$metadata)
  expect_equal_tibble(import_taxonomy(y, otus = hmp5$otus),    hmp5$taxonomy)
  
  # Sample / OTU names in rownames()
  x <- as.data.frame(hmp5$metadata[,-1]); rownames(x) <- hmp5$samples
  y <- as.data.frame(hmp5$taxonomy)[,-1]; rownames(y) <- hmp5$otus
  expect_equal_tibble(import_metadata(x, sids = hmp5$samples), hmp5$metadata)
  expect_equal_tibble(import_taxonomy(y, otus = hmp5$otus),    hmp5$taxonomy)
  
  # Missing or overabundant row ids.
  x <- as.data.frame(hmp5$metadata); rownames(x) <- rev(x[[1]])
  y <- as.data.frame(hmp5$taxonomy); rownames(y) <- rev(y[[1]])
  expect_error(import_metadata(x, sids = hmp5$samples))
  expect_error(import_taxonomy(y, otus = hmp5$otus))
  expect_error(import_metadata(hmp5$metadata[,-1], sids = hmp5$samples))
  expect_error(import_taxonomy(hmp5$taxonomy[,-1], otus = hmp5$otus))
  
  # Bad column name.
  x <- rename(hmp5$metadata, .BMI     = BMI)
  y <- rename(hmp5$taxonomy, .Kingdom = Kingdom)
  expect_error(import_metadata(x, sids = hmp5$samples))
  expect_error(import_taxonomy(y, otus = hmp5$otus))
  
  # No overlap in row ids.
  x <- mutate(hmp5$metadata, .sample = paste0('x', .sample))
  y <- mutate(hmp5$taxonomy,    .otu = paste0('x', .otu))
  expect_error(import_metadata(x, sids = hmp5$samples))
  expect_error(import_taxonomy(y, otus = hmp5$otus))
  
  # Incomplete overlap in row ids.
  expect_warning(import_metadata(hmp5$metadata, sids = hmp5$samples[-1]))
  expect_warning(import_taxonomy(hmp5$taxonomy, otus = hmp5$otus[-1]))
  expect_warning(import_metadata(hmp5$metadata, sids = c('x', hmp5$samples)))
  expect_warning(import_taxonomy(hmp5$taxonomy, otus = c('x', hmp5$otus)))
  
  # Duplicate row ids with different metadata.
  x <- hmp5$metadata; x$.sample[[1]] <- x$.sample[[2]]
  y <- hmp5$taxonomy; y$.otu[[1]]    <- y$.otu[[2]]
  expect_warning(import_metadata(x, sids = hmp5$samples[-1]))
  expect_warning(import_taxonomy(y, otus = hmp5$otus[-1]))
  
  # Duplicate row ids with identical metadata, quietly merged.
  x <- rbind(hmp5$metadata, hmp5$metadata)
  y <- rbind(hmp5$taxonomy, hmp5$taxonomy[1:5,])
  expect_equal_tibble(import_metadata(x,     sids = hmp5$samples), hmp5$metadata)
  expect_equal_tibble(import_taxonomy(y,     otus = hmp5$otus),    hmp5$taxonomy)
  expect_equal_tibble(import_metadata(x[,1], sids = hmp5$samples), hmp5$metadata[,1])
  expect_equal_tibble(import_taxonomy(y[,1], otus = hmp5$otus),    hmp5$taxonomy[,1])
  
})
