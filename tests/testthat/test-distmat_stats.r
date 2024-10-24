test_that("distmat_stats", {
  dm <- bdiv_distmat(rare50)
  x <- distmat_stats(dm, pull(rare50, 'Sex'), test = 'none')
  x <- distmat_stats(dm, pull(rare50, 'Sex'))
  x <- distmat_stats(dm, pull(rare50, 'Sex'))
  expect_identical_json(x, '[{".n":49,".stat":1.6041,".z":1.3505,".p.val":0.087}]')
})
