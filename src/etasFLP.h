#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/RS.h>

void F77_SUB(density2parallel)(double *x,double *y,int *m,double  *xkern,double *ykern,int *nkern,double *h,double *w,double  *hvarx,double *hvary,double *dens);
void F77_SUB(etasfull8newparallel)(int *tflag,int *n,double  *mu,double *k,double *c,double *p,double *g,double *d,double  *q,double *x,double *y,double *t,double *m,double *predictor,double  *l);
void F77_SUB(density2serial)(double *x,double *y,int *m,double  *xkern,double *ykern,int *nkern,double *h,double *w,double  *hvarx,double *hvary,double *dens);
void F77_SUB(etasfull8fast)(int *tflag,int *n,double *mu,double  *k,double *c,double *p,double *g,double *d,double *q,double  *x,double *y,double *t,double *m,double *predictor,int *ind,int *nindex,int *index,double *l);
void F77_SUB(etasfull8newserial)(int *tflag,int *n,double *mu,double  *k,double *c,double *p,double *g,double *d,double *q,double  *x,double *y,double *t,double *m,double *predictor,double *l);
void F77_SUB(etasfull8tintegratednew)(int *n,double *mu,double  *k,double *c,double *p,double *g,double *d,double *q,double  *x,double *y,double *t,double *m,double *predictor,double *l,int  *ngridtot,double *xgrid,double *ygrid,double *tmax);
void F77_SUB(etasfull8tfixednew)(int *n,double *mu,double *k,double  *c,double *p,double *g,double *d,double *q,double *x,double  *y,double *t,double *m,double *predictor,double *l,int  *ngridtot,double *xgrid,double *ygrid,double *tmax);
void F77_SUB(deltafl1kspacevar)(double *x,double *t,double *w,int  *n,int *k,int *m1,int *m2,int *nh,double *rangex,double *h,double *hdef,double *dens,double *integr,double *delta,double *expweight,int  *indweight,int *allocationerr);



