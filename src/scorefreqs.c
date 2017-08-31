#include<R.h>
#include<Rinternals.h>
#include<Rmath.h>

SEXP scorefreqs(SEXP FMS){
  int nsnps, s = 2, i, j;
  double *u, *v, aux, *fms = REAL(FMS);
  SEXP RES;

  nsnps = LENGTH(FMS);
  PROTECT(RES = allocVector(REALSXP,2*nsnps+1));
  u = REAL(RES);
  v = (double *)malloc((2*nsnps+1)*sizeof(double));     

  for(i=0;i<=s*2;i++){
    aux = 0;
    for(j=0;j<=i;j++)
      aux += dbinom(j,2,fms[0],0)*dbinom(i-j,2,fms[1],0);
    u[i] = aux;
  }

  while(s<nsnps){
    s++;
    for(i=0;i<=s*2;i++){
      aux = 0;
      for(j=0;(j<=i)&(j<=2*(s-1));j++)
	aux += u[j]*dbinom(i-j,2,fms[s-1],0);
      v[i] = aux;
      if(i==s*2)
	for(j=0;j<=i;j++)
	  u[j] = v[j];
    }
  }
  
  free(v);
  UNPROTECT(1);
  return RES;
}
