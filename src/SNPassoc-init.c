#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(permutation)(int *, int *, int *, int *, int *,
                     int *, int *, double *);
extern void F77_NAME(wgassociation)(int *, int *, int *, int *, int *,
                     int *, double *);

static const R_FortranMethodDef FortranEntries[] = {
  {"permutation",   (DL_FUNC) &F77_NAME(permutation),   8},
  {"wgassociation", (DL_FUNC) &F77_NAME(wgassociation), 7},
  {NULL, NULL, 0}
};

void R_init_SNPassoc(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}