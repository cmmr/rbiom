test_that("read_tree", {
  
  tree <- "(t9,((t5,t2),(((t10,(t7,t4)),(t6,(t3,t1))),t8)));"
  tree <- expect_s3_class(read_tree(tree), 'phylo')
  expect_silent(tree_subset(hmp5$tree, 1:50))
  
  skip_on_cran()
  
  expect_error(tree_subset(tree, c('a', 'b')))
  expect_error(tree_subset(tree, c(TRUE, TRUE)))
  expect_silent(tree_subset(tree, !logical(10)))
  expect_error(tree_subset(tree, 1:20))
  expect_error(tree_subset(tree, -1:5))
  expect_error(tree_subset(tree, 2:5 / 2))
  expect_error(tree_subset(tree, tree))
  
})
