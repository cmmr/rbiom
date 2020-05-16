library(rbiom)
library(testthat)

RcppParallel::setThreadOptions(numThreads = 2)

test_check("rbiom")
