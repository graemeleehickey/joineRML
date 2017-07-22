#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _joineRML_bSim(SEXP, SEXP, SEXP);
extern SEXP _joineRML_expWArma(SEXP, SEXP, SEXP, SEXP);
extern SEXP _joineRML_gammaUpdate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _joineRML_gammaUpdate_approx(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _joineRML_hazHat(SEXP, SEXP, SEXP);
extern SEXP _joineRML_lambdaUpdate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _joineRML_mvrnormArma(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_joineRML_bSim",               (DL_FUNC) &_joineRML_bSim,                3},
  {"_joineRML_expWArma",           (DL_FUNC) &_joineRML_expWArma,            4},
  {"_joineRML_gammaUpdate",        (DL_FUNC) &_joineRML_gammaUpdate,        11},
  {"_joineRML_gammaUpdate_approx", (DL_FUNC) &_joineRML_gammaUpdate_approx, 10},
  {"_joineRML_hazHat",             (DL_FUNC) &_joineRML_hazHat,              3},
  {"_joineRML_lambdaUpdate",       (DL_FUNC) &_joineRML_lambdaUpdate,       10},
  {"_joineRML_mvrnormArma",        (DL_FUNC) &_joineRML_mvrnormArma,         3},
  {NULL, NULL, 0}
};

void R_init_joineRML(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
