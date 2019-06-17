#include "npcure.h"

SEXP probcurenp0hboot(SEXP Datat,
		      SEXP Datax,
		      SEXP Datadelta,
		      SEXP Nrow,
		      SEXP X,
		      SEXP Lx,
		      SEXP H,
		      SEXP Lh,
		      SEXP Pilot,
		      SEXP Probcurepilot,
		      SEXP B) {
  int conth, contb, contx, i, j, nrow = asInteger(Nrow), lx = asInteger(Lx), lh = asInteger(Lh), b = asInteger(B), *pdatadelta, deltaboot[nrow], torder[nrow], *pdeltaboot, indexmsemin, dummy = 1, index, *ponemdatadelta;
  double pdataxi, pilotcontx, wijp, *pdatat, *pdatax, *px, *ph, *ppilot, *pprobcurepilot, *pxboot, *ptboot, mse[lh], auxmsemin, *phboot, *ptbootpre, probi[nrow], cumprobi, sumprobi, unif, w[nrow][nrow], mseterm, pilotterm;
  SEXP probcureboot = PROTECT(allocVector(REALSXP, lh));
  SEXP deltabootorder = PROTECT(allocVector(INTSXP, nrow));
  SEXP xbootorder = PROTECT(allocVector(REALSXP, nrow));
  SEXP tbootpre = PROTECT(allocVector(REALSXP, nrow));
  SEXP tbootorder = PROTECT(allocVector(REALSXP, nrow));
  SEXP x0 = PROTECT(allocVector(REALSXP, 1));
  SEXP lx0 = PROTECT(ScalarInteger(dummy));
  px = REAL(X);
  pdatat = REAL(Datat);
  pdatax = REAL(Datax);
  pdatadelta = INTEGER(Datadelta);
  ph = REAL(H);
  ppilot = REAL(Pilot);
  pxboot = REAL(xbootorder);
  ptbootpre = REAL(tbootpre);
  ptboot = REAL(tbootorder);
  pdeltaboot = INTEGER(deltabootorder);
  SEXP local = PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(local)[0] = FALSE; 
  pprobcurepilot = REAL(Probcurepilot);			     
  SEXP hboot = PROTECT(allocVector(REALSXP, lx));
  phboot = REAL(hboot);
  SEXP onemdatadelta = PROTECT(allocVector(INTSXP, nrow));
  ponemdatadelta = INTEGER(onemdatadelta);
  for (i = 0; i < nrow; i++)
    ponemdatadelta[i] = 1 - pdatadelta[i];
  GetRNGstate();
  for (contx = 0; contx < lx; contx++) { // 1. for
    pilotterm = pprobcurepilot[contx];
    memset(mse, 0.0, lh*sizeof(double));
    REAL(x0)[0] = px[contx];
    pilotcontx = ppilot[contx];
    for (i = 0; i < nrow; i++) {
      pdataxi = pdatax[i];
      for (j = 0; j < nrow; j++) {
	wijp = fabs(pdataxi - pdatax[j])/pilotcontx;
	w[i][j] = (wijp < 1.0) ? 1.0 - wijp*wijp : 0.0;
      }
    }
    for (contb = 0; contb < b; contb++) { // 2. for
      for (i = 0; i < nrow; i++) { // 3. for
	sumprobi = 0.0;
	for (j = 0; j < nrow; j++) {
	  probi[j] = w[i][j];
	  sumprobi += probi[j];
	}
	unif = runif(0, sumprobi);
	cumprobi = 0.0;
	index = 0;       
	while (cumprobi < unif && index < nrow) {
	  cumprobi += probi[index];  
	  index++;
	}
	ptbootpre[i] = pdatat[index - 1];
	deltaboot[i] = pdatadelta[index - 1];
      } // end of 3. for (in i)
      R_orderVector(torder, nrow, PROTECT(Rf_lang2(tbootpre, onemdatadelta)), TRUE, FALSE);
      for (i = 0; i < nrow; i++) {
	ptboot[i] = ptbootpre[torder[i]];
	pxboot[i] = pdatax[torder[i]];
	pdeltaboot[i] = deltaboot[torder[i]];
      }
      probcureboot = probcurenp0(tbootorder, xbootorder, deltabootorder, Nrow, x0, lx0, H, Lh, local);
      for (conth = 0; conth < lh; conth++) {
	mseterm = REAL(probcureboot)[conth] - pilotterm;
	mse[conth] += mseterm*mseterm;
      }
      UNPROTECT(1);
    } // end of 2. for (in contb)
    auxmsemin = mse[0];
    indexmsemin = 0;
    for (conth = 1; conth < lh; conth++) {
      if (mse[conth] < auxmsemin) {
	auxmsemin = mse[conth];
	indexmsemin = conth;
      }
    }
    phboot[contx] = ph[indexmsemin];
  } // end of 1. for (in contx)
  PutRNGstate();
  UNPROTECT(10);
  return hboot;
}

inline static double ise(double *latencyboot,
			 double *latencypilot,
			 double *t,
			 int lt,
			 double tmax) {
  int i = 1;
  double result, aux;
  bool notdone = true;
  aux = latencyboot[0] - latencypilot[0];
  result = t[0]*aux*aux;
  while (i < lt && notdone) {
    aux = latencyboot[i] - latencypilot[i];
    if (t[i] > tmax) {
      result += (tmax - t[i-1])*aux*aux;
      notdone = false;
    }
    else {
      result += (t[i] - t[i-1])*aux*aux;
    }
    i++;
  }
  return result;
}

SEXP latencynp0hboot(SEXP Datat,
		     SEXP Datax,
		     SEXP Datadelta,
		     SEXP Nrow,
		     SEXP X,
		     SEXP Lx,
		     SEXP H,
		     SEXP Lh,
		     SEXP Pilot,
		     SEXP Probcurepilot,
		     SEXP Latencypilot,
		     SEXP B,
		     SEXP Tmax) {
  int conth, contb, contx, i, j, nrow = asInteger(Nrow), lx = asInteger(Lx), lh = asInteger(Lh), b = asInteger(B), *pdatadelta, deltaboot[nrow], torder[nrow], *pdeltaboot, indexmisemin, dummy = 1, cprobi1, index0, index1, *ponemdatadelta;
  double pdataxi, pilotxi, *pdatat, *pdatax, *px, *ph, *ppilot, *pprobcurepilot, *pxboot, *ptboot, auxmisemin, *phboot, temp0, temp1, tempcens0, tempcens1 = 1.0, *cboot, unif, xgij, Bx[nrow], sumBx = 0.0, *yboot, *ptbootpre, probi0[nrow], probi1[nrow], *platencypilot, *platencyboot, *miseh, sumprobi0, cumprobi0, sumprobi1[nrow], cumprobi1, probi1array[nrow][nrow], tmax = asReal(Tmax);
  bool probi1zero;
  SEXP latencyboot = PROTECT(allocVector(VECSXP, lh));
  SEXP deltabootorder = PROTECT(allocVector(INTSXP, nrow));
  SEXP xbootorder = PROTECT(allocVector(REALSXP, nrow));
  SEXP x0 = PROTECT(allocVector(REALSXP, 1));
  SEXP lx0 = PROTECT(ScalarInteger(dummy));
  px = REAL(X); 
  pdatat = REAL(Datat);
  pdatax = REAL(Datax);
  pdatadelta = INTEGER(Datadelta);
  ph = REAL(H);
  ppilot = REAL(Pilot);
  pxboot = REAL(xbootorder);
  cboot = malloc(nrow*sizeof(double));
  yboot = malloc(nrow*sizeof(double));
  miseh = malloc(lh*sizeof(double));
  SEXP tbootpre = PROTECT(allocVector(REALSXP, nrow));
  ptbootpre = REAL(tbootpre);
  SEXP tbootorder = PROTECT(allocVector(REALSXP, nrow));
  ptboot = REAL(tbootorder);
  pdeltaboot = INTEGER(deltabootorder);
  SEXP local = PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(local)[0] = FALSE;			     
  SEXP hboot = PROTECT(allocVector(REALSXP, lx));
  phboot = REAL(hboot);
  pprobcurepilot = REAL(Probcurepilot);
  for (i = 0; i < nrow; i++) {
    pilotxi = ppilot[i];
    pdataxi = pdatax[i];
    tempcens0 = tempcens1;
    if (pdatadelta[i] == 0)
      tempcens1 *= 1.0 - 1.0/(nrow - i + 1);
    probi0[i] = tempcens0 - tempcens1;
    sumBx = 0.0;
    for (j = 0; j < nrow; j++) {
      xgij = fabs(pdataxi - pdatax[j])/pilotxi;
      Bx[j] = (xgij < 1.0) ? 1.0 - xgij*xgij : 0.0;
      sumBx += Bx[j];
    }
    temp1 = 1;
    for (j = 0; j < nrow; j++) {
      temp0 = temp1;
      if (sumBx > 0 && pdatadelta[j] == 1)
	temp1 *= 1.0 - Bx[j]/sumBx;
      sumBx -= Bx[j];
      probi1array[i][j] = temp0 - temp1;
    }
    sumprobi1[i] = 1 - temp1;
  }
  sumprobi0 = 1 - tempcens1;
  memset(deltaboot, 0, nrow*sizeof(int));
  SEXP onemdatadelta = PROTECT(allocVector(INTSXP, nrow));
  ponemdatadelta = INTEGER(onemdatadelta);
  for (i = 0; i < nrow; i++)
    ponemdatadelta[i] = 1 - pdatadelta[i];
  GetRNGstate();
  for (contx = 0; contx < lx; contx++) { // 1. for
    memset(miseh, 0.0, lh*sizeof(double));
    REAL(x0)[0] = px[contx];
    platencypilot = REAL(VECTOR_ELT(Latencypilot, contx));
    for (contb = 0; contb < b; contb++) { // 2. for
      for (i = 0; i < nrow; i++) { // 3. for
	unif = runif(0, sumprobi0);
	cumprobi0 = 0; 
	index0 = 0; 
	while (cumprobi0 < unif && index0 < nrow) {
	  cumprobi0 += probi0[index0]; 
	  index0++;
	}
	cboot[i] = pdatat[index0 - 1];
	unif = runif(0, 1);
	if (unif < pprobcurepilot[i]) {
	  ptbootpre[i] = cboot[i];
	  deltaboot[i] = 0;
	}
	else {
	  for (j = 0; j < nrow; j++)
	    probi1[j] = probi1array[i][j];
	  probi1zero = true;
	  cprobi1 = 0;
	  while (probi1zero && cprobi1 < nrow) {
	    if (probi1[cprobi1] > 0)
	      probi1zero = false;
	    cprobi1++;
	  }
	  if (probi1zero) {
	    yboot[i] = pdatat[i];
	  }
	  else {
	    unif = runif(0, sumprobi1[i]);
	    cumprobi1 = 0;
	    index1 = 0;       
	    while (cumprobi1 < unif && index1 < nrow) {
	      cumprobi1 += probi1[index1];  
	      index1++;
	    }
	    yboot[i] = pdatat[index1 - 1];
	  }
	  ptbootpre[i] = fmin(yboot[i], cboot[i]);
	  deltaboot[i] = (yboot[i] <= cboot[i]) ? 1 : 0;
	}
      } // end of 3. for (in i)
      R_orderVector(torder, nrow, PROTECT(Rf_lang2(tbootpre, onemdatadelta)), TRUE, FALSE);
      for (i = 0; i < nrow; i++) {
	pxboot[i] = pdatax[torder[i]];
	ptboot[i] = ptbootpre[torder[i]];
	pdeltaboot[i] = deltaboot[torder[i]];
      }
      latencyboot = latencynp0(tbootorder, xbootorder, deltabootorder, Nrow, x0, lx0, H, Lh, local, Datat, Nrow);
      for (conth = 0; conth < lh; conth++) {
	platencyboot = REAL(VECTOR_ELT(latencyboot, conth));
	miseh[conth] += ise(platencyboot, platencypilot, pdatat, nrow, tmax);
      }
      UNPROTECT(1);
    } // end of 2. for (in contb)
    auxmisemin = miseh[0];
    indexmisemin = 0;
    for (conth = 1; conth < lh; conth++) {
      if (miseh[conth] < auxmisemin) {
	auxmisemin = miseh[conth];
	indexmisemin = conth;
      }
    }
    phboot[contx] = ph[indexmisemin];
  } // end of 1. for (in contx)
  PutRNGstate(); 
  free(cboot);
  free(yboot);
  free(miseh);
  UNPROTECT(10);
  return hboot;
}
