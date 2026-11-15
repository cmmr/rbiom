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
# The count is read and set with ecodive::blas_threads(). With an ecodive
# that doesn't have it, or a library whose count can't be controlled
# (R's reference BLAS is single-threaded; Apple's Accelerate can't be
# limited at run time), blas_limit() does nothing, silently.
#________________________________________________________
blas_limit <- function (cpus, env = parent.frame()) {

  if (!is_scalar_integerish(cpus) || is.na(cpus) || cpus < 1) return (invisible())
  if (!blas_ctl()) return (invisible())

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


#________________________________________________________
# ecodive::blas_threads(), or NULL with an ecodive version that doesn't
# have it. Looked up when called, so that rbiom works with either.
#________________________________________________________
blas_fun <- function () {
  tryCatch(
    expr  = getExportedValue('ecodive', 'blas_threads'),
    error = function (e) NULL )
}


#________________________________________________________
# Whether the library's thread count can be read and set.
#________________________________________________________
blas_ctl <- function () !is.null(blas_fun())


#________________________________________________________
# Get the library's thread count (NA when it can't be read), or set it
# to `n`, a whole number >= 1; anything else sets nothing. Never an
# error or a message.
#________________________________________________________
blas_threads <- function (n = NULL) {

  fun <- blas_fun()

  if (is.null(n)) {
    n <- if (is.null(fun)) NA else tryCatch(fun(), error = function (e) NA)
    return (if (is_scalar_integerish(n) && !is.na(n)) as.integer(n) else NA_integer_)
  }

  if (!is.null(fun) && is_scalar_integerish(n) && !is.na(n) && n >= 1)
    tryCatch(fun(as.integer(n)), error = function (e) NULL)

  return (invisible())
}
