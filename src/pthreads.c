#include <R.h>
#include <Rinternals.h>

// Detect if pthread is available.
#if defined __has_include
#  if __has_include (<pthread.h>)
#    define HAVE_PTHREAD
#  endif
#endif


//======================================================
// R interface. Returns if pthreads are available.
//======================================================
SEXP C_pthreads(void) {
#ifdef HAVE_PTHREAD
  return ScalarLogical(1);
#else
  return ScalarLogical(0);
#endif
}
