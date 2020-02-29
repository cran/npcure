#include "npcure.h"

// Inline function called by latencynp0
inline static void latencynp0core(double *pdatax,
				  double *pdatat,
				  int *pdatadelta,
				  int nrow,
				  double x,
				  double h,
				  int lt,
				  bool testim,
				  double *ppt,
				  double *presult) {
  int i, lasti, j;
  double xhi, Bx[nrow], sumBx = 0.0, temp = 1.0, pberan[nrow], ps[nrow];
  bool found;
  for (i = 0; i < nrow; i++) {
    xhi = fabs(x - pdatax[i])/h;
    Bx[i] = (xhi < 1.0) ? 1.0 - xhi*xhi : 0.0;
    sumBx += Bx[i];
  }
  for (i = 0; i < nrow; i++) {
    if (sumBx > 0 && pdatadelta[i] == 1)
      temp *= 1.0 - Bx[i]/sumBx;
    // very small values are set to 0
    if (temp < 1e-10)
      temp = 0.0;
    pberan[i] = temp;
    sumBx -= Bx[i];
  }
  for (i = 0; i < nrow; i++)
    ps[i] = (1.0 - temp < 1e-10) ? 0.0 : (pberan[i] - temp)/(1.0 - temp);
  if (testim) {
    lasti = 0;
    for (j = 0; j < lt; j++) {
      found = false;
      i = lasti;
      while (!found && i < nrow) {
	if (ppt[j] < pdatat[i]) {
	  presult[j] = (i == 0) ? 1.0 : ps[i - 1];
	  found = true;
	  lasti = i;
	}
	i++;
      }
      if (!found)
	presult[j] = ps[nrow - 1];
    }
  }
  else {
    for (i = 0; i < nrow; i++)
      presult[i] = ps[i];
  }
}

// NP estimator of the latency function
SEXP latencynp0(SEXP Datat,
		SEXP Datax,
		SEXP Datadelta,
		SEXP Nrow,
		SEXP X,
		SEXP Lx,
		SEXP H,
		SEXP Lh,
		SEXP Local,
		SEXP T,
		SEXP Lt) {
  int conth, i, j, k, nrow = asInteger(Nrow), lx = asInteger(Lx), lh = asInteger(Lh), *pdatadelta, lt = asInteger(Lt), lastj, nroworlt;
  double x, xhj, h, Bx[nrow], sumBx = 0.0, temp = 1.0, *px, *pdatax, *ph, *presult, *pberan, *pdatat, *ppt = NULL, *ps;
  bool found, testim = false, local = LOGICAL(Local)[0];
  px = REAL(X);
  pdatax = REAL(Datax);
  pdatadelta = INTEGER(Datadelta);
  ph = REAL(H);
  pdatat = REAL(Datat);
  if (!isNull(T)) {
    ppt = REAL(T);
    testim = true;
  }
  nroworlt = testim ? lt : nrow;
  if (local) {
    SEXP resultlist = PROTECT(allocVector(VECSXP, lx));
    for (i = 0; i < lx; i++) {
      SEXP result = PROTECT(allocVector(REALSXP, nroworlt));
      presult = REAL(result);
      h = ph[i];
      x = px[i];
      latencynp0core(pdatax, pdatat, pdatadelta, nrow, x, h, lt, testim, ppt, presult);
      SET_VECTOR_ELT(resultlist, i, result);
      UNPROTECT(1);
    }
    UNPROTECT(1);
    return resultlist;
  }
  else { // local == 0
    if (lx == 1) { // in this case, for efficiency, the inline function is not called
      ps = malloc(nrow*sizeof(double));
      pberan = malloc(nrow*sizeof(double));
      SEXP resultlist = PROTECT(allocVector(VECSXP, lh));
      x = px[0];
      for (conth = 0; conth < lh; conth++) {
	SEXP result = PROTECT(allocVector(REALSXP, nroworlt));
	presult = REAL(result);
	h = ph[conth];
	for (j = 0; j < nrow; j++) {
	  xhj = fabs(x - pdatax[j])/h;
	  Bx[j] = (xhj < 1.0) ? 1.0 - xhj*xhj : 0.0;
	  sumBx += Bx[j];
	}
	for (j = 0; j < nrow; j++) {
	  if (sumBx > 0 && pdatadelta[j] == 1)
	    temp *= 1.0 - Bx[j]/sumBx;
	  // very small values are set to 0
	  if (temp < 1e-10)
	    temp = 0.0;
	  pberan[j] = temp;
	  sumBx -= Bx[j];
	}
	for (j = 0; j < nrow; j++)
	  ps[j] = (1.0 - temp < 1e-10) ? 0.0 : (pberan[j] - temp)/(1.0 - temp);
	if (testim) {
	  lastj = 0;
	  for (k = 0; k < lt; k++) {
	    found = false;
	    j = lastj;
	    while (!found && j < nrow) {
	      if (ppt[k] < pdatat[j]) {
		presult[k] = (j == 0) ? 1.0 : ps[j - 1];
		found = true;
		lastj = j;
	      }
	      j++;
	    }
	    if (!found)
	      presult[k] = ps[nrow - 1];
	  }
	}
	else {
	  for (j = 0; j < nrow; j++)
	    presult[j] = ps[j];
	}
	sumBx = 0.0;
	temp = 1.0;
	SET_VECTOR_ELT(resultlist, conth, result);
	UNPROTECT(1);
      }
      free(pberan);
      free(ps);
      UNPROTECT(1);
      return resultlist;
    }
    else { // lx > 1
      SEXP resultlistlist = PROTECT(allocVector(VECSXP, lh));
      for (conth = 0; conth < lh; conth++) {
	SEXP resultlist = PROTECT(allocVector(VECSXP, lx));
	for (i = 0; i < lx; i++) {	  
	  SEXP result = PROTECT(allocVector(REALSXP, nroworlt));
	  presult = REAL(result);
	  h = ph[conth];
	  x = px[i];
	  latencynp0core(pdatax, pdatat, pdatadelta, nrow, x, h, lt, testim, ppt, presult);
	  SET_VECTOR_ELT(resultlist, i, result);
	  UNPROTECT(1);
	}
	SET_VECTOR_ELT(resultlistlist, conth, resultlist);
	UNPROTECT(1);
      }
      UNPROTECT(1);
      return resultlistlist;
    }
  }
}

// Function for computing confidence bands for the NP estimator of latency
SEXP latencynp0confband(SEXP Datat,
			SEXP Datax,
			SEXP Datadelta,
			SEXP Nrow,
			SEXP X,
			SEXP Lx,
			SEXP H,
			SEXP Lh,
			SEXP Pilot,
			SEXP Probcurepilot,
			SEXP Conf,
			SEXP B,
			SEXP S,
			SEXP Local) {
  int contb, conth, contx, i, j, nrow = asInteger(Nrow), lx = asInteger(Lx), b = asInteger(B), *pdatadelta, deltaboot[nrow], torder[nrow], *pdeltaboot, dummy = 1, cprobi1, index0, index1, *ponemdatadelta, lh = asInteger(Lh);
  double pdataxi, pilotxi, *pdatat, *pdatax, *px, *ph, *ppilot, *pprobcurepilot, *pxboot, *ptboot, temp0, temp1, tempcens0, tempcens1 = 1.0, *cboot, unif, xgij, Bx[nrow], sumBx = 0.0, *yboot, *ptbootpre, probi0[nrow], probi1[nrow], *platencyboot, *plow, *pup, sumprobi0, cumprobi0, sumprobi1[nrow], cumprobi1, probi1array[nrow][nrow], conf = asReal(Conf), Si, Sim1, mterm, m[nrow], m2[nrow], sd[nrow], quant = qnorm(conf, 0, 1, TRUE, FALSE);
  bool probi1zero;
  SEXP deltabootorder = PROTECT(allocVector(INTSXP, nrow));
  SEXP xbootorder = PROTECT(allocVector(REALSXP, nrow));
  SEXP x0 = PROTECT(allocVector(REALSXP, 1));
  SEXP lx0 = PROTECT(ScalarInteger(dummy));
  SEXP h0 = PROTECT(allocVector(REALSXP, 1));
  px = REAL(X);
  pdatat = REAL(Datat);
  pdatax = REAL(Datax);
  pdatadelta = INTEGER(Datadelta);
  ph = REAL(H);
  ppilot = REAL(Pilot);
  pprobcurepilot = REAL(Probcurepilot);
  pxboot = REAL(xbootorder);
  cboot = malloc(nrow*sizeof(double));
  yboot = malloc(nrow*sizeof(double));
  SEXP tbootpre = PROTECT(allocVector(REALSXP, nrow));
  ptbootpre = REAL(tbootpre);
  SEXP tbootorder = PROTECT(allocVector(REALSXP, nrow));
  ptboot = REAL(tbootorder);
  pdeltaboot = INTEGER(deltabootorder);
  SEXP localdummy = PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(localdummy)[0] = TRUE;
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
  SEXP onemdatadelta = PROTECT(allocVector(INTSXP, nrow));
  ponemdatadelta = INTEGER(onemdatadelta);
  for (i = 0; i < nrow; i++)
    ponemdatadelta[i] = 1 - pdatadelta[i];
  GetRNGstate();
  if (LOGICAL(Local)[0]) {
    SEXP confbandlist = PROTECT(allocVector(VECSXP, lx));
    for (contx = 0; contx < lx; contx++) { // 1. for (in contx)
      SEXP confband = PROTECT(allocVector(VECSXP, 2));
      SEXP low = PROTECT(allocVector(REALSXP, nrow));
      plow = REAL(low);
      SEXP up = PROTECT(allocVector(REALSXP, nrow));
      pup = REAL(up);
      SEXP latencyboot = PROTECT(allocVector(REALSXP, nrow));
      REAL(x0)[0] = px[contx];
      REAL(h0)[0] = ph[contx];
      memset(m, 0, nrow*sizeof(double));
      memset(m2, 0, nrow*sizeof(double));
      for (contb = 0; contb < b; contb++) { // 2. for (in contb)
	for (i = 0; i < nrow; i++) { // 3. for (in i)
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
	} // of 3. for (in i)
	R_orderVector(torder, nrow, PROTECT(Rf_lang2(tbootpre, onemdatadelta)), TRUE, FALSE);
	for (i = 0; i < nrow; i++) {
	  pxboot[i] = pdatax[torder[i]];
	  ptboot[i] = ptbootpre[torder[i]];
	  pdeltaboot[i] = deltaboot[torder[i]];
	}
	latencyboot = latencynp0(tbootorder, xbootorder, deltabootorder, Nrow, x0, lx0, h0, lx0, localdummy, Datat, Nrow);
	platencyboot = REAL(VECTOR_ELT(latencyboot, 0));
	for (i = 0; i < nrow; i++) {
	  mterm = platencyboot[i];
	  m[i] += mterm;
	  m2[i] += mterm*mterm;
	}
	UNPROTECT(1);
      } // end of 2. for (in contb)
      for (i = 0; i < nrow; i++) {
	sd[i] = sqrt((m2[i] - m[i]*m[i]/b)/(b - 1));
	Si = REAL(VECTOR_ELT(S, contx))[i];
	if (Si < 1e-10) {
	  plow[i] = nan("");
	  pup[i] = nan("");
	}
	else {
	  plow[i] = fmax(fmin(Si - quant*sd[i], 1), 0);
	  pup[i] = fmax(fmin(Si + quant*sd[i], 1),0);
	  if (plow[i] > Si)
	    plow[i] = Si;
	  if (pup[i] < Si)
	    pup[i] = Si;
	  if (i > 0 && plow[i] > plow[i - 1])
	    plow[i] = plow[i - 1];
	  if (i > 0 && pup[i] > pup[i - 1])
	    pup[i] = pup[i - 1];
	}
      }
      for (i = nrow - 1; i > 1; i--) {
	Si = REAL(VECTOR_ELT(S, contx))[i];
	Sim1 = REAL(VECTOR_ELT(S, contx))[i - 1];
	if (Si == Sim1) {
	  plow[i - 1] = plow[i];
	  pup[i - 1] = pup[i]; 
	}
      }
      SET_VECTOR_ELT(confband, 0, low);
      SET_VECTOR_ELT(confband, 1, up);
      SET_VECTOR_ELT(confbandlist, contx, confband);
      UNPROTECT(4);
    } // end of 1. for (in contx)
    PutRNGstate();
    free(cboot);
    free(yboot);
    UNPROTECT(10);
    return confbandlist;
  }
  else { // Local == FALSE
    SEXP confbandlistlist = PROTECT(allocVector(VECSXP, lh));
    for (conth = 0; conth < lh; conth++) { // 1. for (in conth)
      SEXP confbandlist = PROTECT(allocVector(VECSXP, lx));
      REAL(h0)[0] = ph[conth];
      for (contx = 0; contx < lx; contx++) { // 2. for (in contx)
	SEXP confband = PROTECT(allocVector(VECSXP, 2));
	SEXP low = PROTECT(allocVector(REALSXP, nrow));
	plow = REAL(low);
	SEXP up = PROTECT(allocVector(REALSXP, nrow));
	pup = REAL(up);
	SEXP latencyboot = PROTECT(allocVector(REALSXP, nrow));
	REAL(x0)[0] = px[contx];
	memset(m, 0, nrow*sizeof(double));
	memset(m2, 0, nrow*sizeof(double));
	for (contb = 0; contb < b; contb++) { // 3. for (in contb)
	  for (i = 0; i < nrow; i++) { // 4. for (in i)
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
	  } // of 4. for (in i)
	  R_orderVector(torder, nrow, PROTECT(Rf_lang2(tbootpre, onemdatadelta)), TRUE, FALSE);
	  for (i = 0; i < nrow; i++) {
	    pxboot[i] = pdatax[torder[i]];
	    ptboot[i] = ptbootpre[torder[i]];
	    pdeltaboot[i] = deltaboot[torder[i]];
	  }
	  latencyboot = latencynp0(tbootorder, xbootorder, deltabootorder, Nrow, x0, lx0, h0, lx0, localdummy, Datat, Nrow);
	  platencyboot = REAL(VECTOR_ELT(latencyboot, 0));
	  for (i = 0; i < nrow; i++) {
	    mterm = platencyboot[i];
	    m[i] += mterm;
	    m2[i] += mterm*mterm;
	  }
	  UNPROTECT(1);
	} // end of 3. for (in contb)
	for (i = 0; i < nrow; i++) {
	  sd[i] = sqrt((m2[i] - m[i]*m[i]/b)/(b - 1));
	  Si = REAL(VECTOR_ELT(VECTOR_ELT(S, conth), contx))[i]; 
	  if (Si < 1e-10) {
	    plow[i] = nan("");
	    pup[i] = nan("");
	  }
	  else {
	    plow[i] = fmax(fmin(Si - quant*sd[i], 1), 0);
	    pup[i] = fmax(fmin(Si + quant*sd[i], 1), 0);
	    if (plow[i] > Si)
	      plow[i] = Si;
	    if (pup[i] < Si)
	      pup[i] = Si;
	    if (i > 0 && plow[i] > plow[i - 1])
	      plow[i] = plow[i - 1];
	    if (i > 0 && pup[i] > pup[i - 1])
	      pup[i] = pup[i - 1];
	  }
	}
	for (i = nrow - 1; i > 1; i--) {
	  Si = REAL(VECTOR_ELT(VECTOR_ELT(S, conth), contx))[i];
	  Sim1 = REAL(VECTOR_ELT(VECTOR_ELT(S, conth), contx))[i - 1];
	  if (Si == Sim1) {
	    plow[i - 1] = plow[i];
	    pup[i - 1] = pup[i]; 
	  }
	}
	SET_VECTOR_ELT(confband, 0, low);
	SET_VECTOR_ELT(confband, 1, up);
	SET_VECTOR_ELT(confbandlist, contx, confband);
	UNPROTECT(4);
      } // end of 2. for (in contx)
      SET_VECTOR_ELT(confbandlistlist, conth, confbandlist);
      UNPROTECT(1);
    } // end of 1. for (in conth)
    PutRNGstate();
    free(cboot);
    free(yboot);
    UNPROTECT(10);
    return confbandlistlist;
  }
}
