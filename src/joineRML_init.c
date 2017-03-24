#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP joineRML_bSim(SEXP, SEXP, SEXP);
extern SEXP joineRML_expWArma(SEXP, SEXP, SEXP, SEXP);
extern SEXP joineRML_gammaUpdate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP joineRML_gammaUpdate_approx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP joineRML_hazHat(SEXP, SEXP, SEXP);
extern SEXP joineRML_lambdaUpdate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP joineRML_mvrnormArma(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"joineRML_bSim",               (DL_FUNC) &joineRML_bSim,                3},
  {"joineRML_expWArma",           (DL_FUNC) &joineRML_expWArma,            4},
  {"joineRML_gammaUpdate",        (DL_FUNC) &joineRML_gammaUpdate,        11},
  {"joineRML_gammaUpdate_approx", (DL_FUNC) &joineRML_gammaUpdate_approx, 10},
  {"joineRML_hazHat",             (DL_FUNC) &joineRML_hazHat,              3},
  {"joineRML_lambdaUpdate",       (DL_FUNC) &joineRML_lambdaUpdate,       10},
  {"joineRML_mvrnormArma",        (DL_FUNC) &joineRML_mvrnormArma,         3},
  {NULL, NULL, 0}
};

void R_init_joineRML(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
