#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "etasFLP.h"

static const R_FortranMethodDef FortEntries[] = {
    {"density2serial", (DL_FUNC) &F77_SUB(density2serial), 11},
    {"etasfull8newserial", (DL_FUNC) &F77_SUB(etasfull8newserial), 15},
    {"etasfull8tintegratednew", (DL_FUNC) &F77_SUB(etasfull8tintegratednew), 18},
    {"etasfull8tfixednew", (DL_FUNC) &F77_SUB(etasfull8tfixednew), 18},
    {"deltafl1kspacevar", (DL_FUNC) &F77_SUB(deltafl1kspacevar), 17},
    {NULL, NULL, 0}
};


void attribute_visible R_init_etasFLP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

