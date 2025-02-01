/* Generated with tools::package_native_routine_registration_skeleton('.',,,FALSE) */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP C_alpha_div(SEXP, SEXP, SEXP);
extern SEXP C_beta_div(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pthreads(void);
extern SEXP C_rarefy(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_read_tree(SEXP, SEXP);
extern SEXP C_unifrac(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_alpha_div", (DL_FUNC) &C_alpha_div, 3},
  {"C_beta_div",  (DL_FUNC) &C_beta_div,  4},
  {"C_pthreads",  (DL_FUNC) &C_pthreads,  0},
  {"C_rarefy",    (DL_FUNC) &C_rarefy,    5},
  {"C_read_tree", (DL_FUNC) &C_read_tree, 2},
  {"C_unifrac",   (DL_FUNC) &C_unifrac,   6},
  {NULL, NULL, 0}
};

void R_init_rbiom(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
