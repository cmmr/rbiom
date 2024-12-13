test_that("s3_methods", {
  
  skip_on_cran()
  
  df <- expect_silent(stats_table(adiv_table(rare50)))
  expect_silent(df$.n)         # $.rbiom_tbl
  expect_output(print(df))     # tbl_sum.rbiom_tbl
  expect_output(print(df$cmd)) # print.rbiom_code
  
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
  expect_silent(na.omit(hmp5))  # na.omit.rbiom
  
  expect_silent(slice(hmp5, 1:2, 4))         # slice.rbiom
  expect_silent(slice_head(hmp5, n = 2))     # slice_head.rbiom
  expect_silent(slice_tail(hmp5, n = 2))     # slice_tail.rbiom
  expect_silent(slice_min(hmp5, Age, n = 2)) # slice_min.rbiom
  expect_silent(slice_max(hmp5, Age, n = 2)) # slice_max.rbiom
  expect_silent(slice_sample(hmp5, n = 2))   # slice_sample.rbiom
  
})
