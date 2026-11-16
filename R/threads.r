#________________________________________________________
# Threads of the linear algebra library.
#
# `cpus` sets how many threads ecodive uses. The BLAS/LAPACK library that
# R is linked to can run threads of its own: OpenBLAS, MKL, and BLIS
# start up to one per core for each large matrix operation (for example
# the eigendecomposition in stats::cmdscale()), whatever `cpus` says.
#
# When `cpus` is less than the library's thread count, blas_limit() sets
# the library to one thread until the function whose frame is `env`
# exits - by returning, with an error, or when interrupted - and then
# puts the previous count back. One thread rather than `cpus`, because
# more did not pay: measured with OpenBLAS 0.3.32 on six cores, more
# threads made a PCoA of 4,000 samples at most 20% faster on an idle
# machine and twice as slow or worse on a busy one, two threads made it
# five times slower than one, and betadisper and the model fits were
# slower with six threads than with one. The default, `cpus` = all
# cores, leaves the library as it is, and so does a count that is
# already within `cpus` (e.g. from OPENBLAS_NUM_THREADS=1).
#
# validate_cpus() calls it, so every exported function that takes `cpus`
# keeps its linear algebra within `cpus` threads.
#
# The count is read and set with ecodive::blas_threads(), which returns
# NA for a library whose count it can't read (R's reference BLAS, which
# is single-threaded, or Apple's Accelerate) and never raises an error,
# warning, or message. Then blas_limit() does nothing.
#________________________________________________________
blas_limit <- function (cpus, env = parent.frame()) {

  if (!is_scalar_integerish(cpus) || is.na(cpus) || cpus < 1) return (invisible())

  # The previous count is put back by on.exit() in `env`, which therefore
  # has to be the frame of a function that is running.
  if (!any(vapply(sys.frames(), identical, logical(1), env))) return (invisible())

  old <- blas_threads()
  if (is.na(old) || old <= cpus) return (invisible())

  # Register the restore before lowering the count, so that an interrupt
  # in between can't leave the count lowered.
  do.call(base::on.exit, list(as.call(list(blas_threads, old)), add = TRUE), envir = env)
  blas_threads(1L)

  return (invisible())
}
