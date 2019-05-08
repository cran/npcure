#include "npcure.h"

// Cross-validation bandwidth for Beran's estimator
SEXP berannp0cv(SEXP Datat,
		SEXP Datax,
		SEXP Datadelta,
		SEXP Nrow,
		SEXP X,
		SEXP Lx,
		SEXP H,
		SEXP Lh) {
  int i, j, contx, conth, nrow = asInteger(Nrow), nrowm1, lx = asInteger(Lx), lh = asInteger(Lh), indexmincvh, dummy = 1, *pdatadelta, *pdatadeltam1;
  double cv = 0.0, auxmincvh, *pdatat, *pdatax, *px, *ph, *pdatatm1, *pdataxm1, *pcvh, *presult, cvterm;
  bool found = false, start = true;
  pdatat = REAL(Datat);
  pdatax = REAL(Datax);
  pdatadelta = INTEGER(Datadelta);
  px = REAL(X);
  ph = REAL(H);
  SEXP local = PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(local)[0] = TRUE;
  SEXP Nrowm1 = PROTECT(allocVector(INTSXP, 1));
  nrowm1 = nrow - 1;
  INTEGER(Nrowm1)[0] = nrowm1;
  SEXP Datatm1 = PROTECT(allocVector(REALSXP, nrowm1));
  pdatatm1 = REAL(Datatm1);
  SEXP Dataxm1 = PROTECT(allocVector(REALSXP, nrowm1));
  pdataxm1 = REAL(Dataxm1);
  SEXP Datadeltam1 = PROTECT(allocVector(INTSXP, nrowm1));
  pdatadeltam1 = INTEGER(Datadeltam1);
  SEXP x0 = PROTECT(allocVector(REALSXP, 1));
  SEXP h0 = PROTECT(allocVector(REALSXP, 1));
  SEXP lx0 = PROTECT(ScalarInteger(dummy));
  SEXP sberan = PROTECT(allocVector(REALSXP, nrow));
  SEXP cvh = PROTECT(allocVector(REALSXP, lh));
  pcvh = REAL(cvh);
  SEXP result = PROTECT(allocVector(REALSXP, lx));
  presult = REAL(result);
  for (contx = 0; contx < lx; contx++) {
    REAL(x0)[0] = px[contx];
    for (conth = 0; conth < lh; conth++) {
      REAL(h0)[0] = ph[conth];
      for (i = 0; i < nrow; i++) {
	for (j = 0; j < i; j++) {
	  pdatatm1[j] = pdatat[j];
	  pdataxm1[j] = pdatax[j];
	  pdatadeltam1[j] = pdatadelta[j];
	}
	for (j = i + 1; j < nrow; j++) {
	  pdatatm1[j - 1] = pdatat[j];
	  pdataxm1[j - 1] = pdatax[j];
	  pdatadeltam1[j - 1] = pdatadelta[j]; 
	}
	sberan = berannp0(Datatm1, Dataxm1, Datadeltam1, Nrowm1, x0, lx0, h0, lx0, local, Datat, Nrow);
	for (j = 0; j < nrow; j++) {
	  if ((pdatadelta[i] == 1 && pdatadelta[j] == 1) || (i == j) || (pdatadelta[i] == 1 && pdatadelta[j] == 0 && pdatat[i] <= pdatat[j]) || (pdatadelta[i] == 0 && pdatadelta[j] == 1 && pdatat[i] >= pdatat[j])) {// 'useful pair'
	    cvterm = REAL(VECTOR_ELT(sberan, 0))[j];
	    if (pdatat[i] <= pdatat[j])
	      cv += 1 - cvterm*(2 - cvterm);
	    else
	      cv += cvterm*cvterm;
	  }
	}
      }
      pcvh[conth] = cv;
      cv = 0.0;
    }
    // search of the largest local minimum
    auxmincvh = pcvh[lh - 1];
    indexmincvh = lh - 1;
    conth = lh - 2;
    while (!found && conth > 0) {
      if (pcvh[conth] < auxmincvh) {
	auxmincvh = pcvh[conth];
	indexmincvh = conth;
	conth--;
	if (start)
	  start = false;
      }
      else {
	if (start) {
	  auxmincvh = pcvh[conth];
	  indexmincvh = conth;
	  conth--;
	}
	else
	  found = true;
      }
    }
    if (!found && start) {
      presult[contx] = ph[lh - 1];      
    }
    else {
      presult[contx] = ph[indexmincvh];
      found = false;
      start = true;
    }
  }
  UNPROTECT(11);
  return result;
}
