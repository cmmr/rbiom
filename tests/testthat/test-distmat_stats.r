test_that("distmat_stats", {
  dm <- bdiv_distmat(rare50)
  x <- distmat_stats(dm, pull(rare50, 'Sex'), test = 'none')
  x <- distmat_stats(dm, pull(rare50, 'Sex'))
  x <- distmat_stats(dm, pull(rare50, 'Sex'))
  expect_identical_json(round(x, 1), '[{".n":49,".stat":1.6,".z":1.4,".p.val":0.1}]')
})
