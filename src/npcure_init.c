#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP berannp0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP berannp0confband(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP berannp0cv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP latencynp0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP latencynp0confband(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP latencynp0hboot(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP probcurenp0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP probcurenp0confband(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP probcurenp0hboot(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"berannp0",            (DL_FUNC) &berannp0,            11},
    {"berannp0confband",    (DL_FUNC) &berannp0confband,    14},
    {"berannp0cv",          (DL_FUNC) &berannp0cv,           8},
    {"latencynp0",          (DL_FUNC) &latencynp0,          11},
    {"latencynp0confband",  (DL_FUNC) &latencynp0confband,  14},
    {"latencynp0hboot",     (DL_FUNC) &latencynp0hboot,     13},
    {"probcurenp0",         (DL_FUNC) &probcurenp0,          9},
    {"probcurenp0confband", (DL_FUNC) &probcurenp0confband, 13},
    {"probcurenp0hboot",    (DL_FUNC) &probcurenp0hboot,    11},
    {NULL, NULL, 0}
};

void R_init_npcure(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
