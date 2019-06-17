#include "npcure.h"

// Inline function called by probcurenp0
inline static double probcurenp0core(double *pdatax,
				     int *pdatadelta,
				     int nrow,
				     double x,
				     double h) {
  int i;
  double xhi, Bx[nrow], sumBx = 0.0, temp = 1.0;
  for (i = 0; i < nrow; i++) {
    xhi = fabs(x - pdatax[i])/h;
    Bx[i] = (xhi < 1.0) ? 1.0 - xhi*xhi : 0.0;
    sumBx += Bx[i];
  }
  for (i = 0; i < nrow; i++) {
    if (sumBx > 0 && pdatadelta[i] == 1)
      temp *= 1.0 - Bx[i]/sumBx;
    sumBx -= Bx[i];
  }
  return fabs(temp);
}

// NP estimator of the conditional probability of cure
SEXP probcurenp0(SEXP Datat, 
		 SEXP Datax,
		 SEXP Datadelta,
		 SEXP Nrow,
		 SEXP X,
		 SEXP Lx,
		 SEXP H,
		 SEXP Lh,
		 SEXP Local) {
  int conth, i, j, nrow = asInteger(Nrow), lx = asInteger(Lx), lh = asInteger(Lh), *pdatadelta;
  double x, h, xhj, sumBx = 0.0, temp = 1.0, Bx[nrow], *px, *pdatax, *ph, *presult;
  bool local = LOGICAL(Local)[0];
  px = REAL(X);
  pdatax = REAL(Datax);
  pdatadelta = INTEGER(Datadelta);
  ph = REAL(H);
  if (local) { // 1st if-else: local bandwidths
    SEXP result = PROTECT(allocVector(REALSXP, lx));
    presult = REAL(result);
    for (i = 0; i < lx; i++) {
      h = ph[i];
      x = px[i];
      presult[i] = probcurenp0core(pdatax, pdatadelta, nrow, x, h);
    }
    UNPROTECT(1);
    return result;
  }
  else { // else of 1st if-else: global bandwidths
    if (lx == 1) { // 2nd if-else
      SEXP result = PROTECT(allocVector(REALSXP, lh));
      presult = REAL(result);
      x = px[0];
      for (conth = 0; conth < lh; conth++) {
	h = ph[conth];
	for (j = 0; j < nrow; j++) {
	  xhj = fabs(x - pdatax[j])/h;
	  Bx[j] = (xhj < 1.0) ? 1.0 - xhj*xhj : 0.0;
	  sumBx += Bx[j];
	}
	for (j = 0; j < nrow; j++) {
	  if (sumBx > 0 && pdatadelta[j] == 1)
	    temp *= 1.0 - Bx[j]/sumBx;
	  sumBx -= Bx[j];
	}
	presult[conth] = fabs(temp);
	sumBx = 0.0;
	temp = 1.0;
      }
      UNPROTECT(1);
      return result;
    }
    else { // lx != 1
      SEXP resultlist = PROTECT(allocVector(VECSXP, lh));
      for (conth = 0; conth < lh; conth++) {
	SEXP result = PROTECT(allocVector(REALSXP, lx));
	presult = REAL(result);
	h = ph[conth];
	for (i = 0; i < lx; i++) {
	  x = px[i];
	  presult[i] = probcurenp0core(pdatax, pdatadelta, nrow, x, h);
	}
	SET_VECTOR_ELT(resultlist, conth, result);
	UNPROTECT(1);
      }
      UNPROTECT(1);
      return resultlist;
    } // end of else of 2nd if-else
  }  // end of else of 1st if-else
}

// Confidence bands for the conditional probability of cure
SEXP probcurenp0confband(SEXP Datat,
			 SEXP Datax,
			 SEXP Datadelta,
			 SEXP Nrow,
			 SEXP X,
			 SEXP Lx,
			 SEXP H,
			 SEXP Lh,
			 SEXP Conf,
			 SEXP B,
			 SEXP Pilot,
			 SEXP Q,
			 SEXP Local) {
  int contb, conth, contx, i, j, nrow = asInteger(Nrow), lx = asInteger(Lx), b = asInteger(B), *pdatadelta, deltaboot[nrow], torder[nrow], *pdeltaboot, dummy = 1, index, *ponemdatadelta, lh = asInteger(Lh);
  double pdataxi, wijp, *pdatat, *pdatax, *px, *ph, *pxboot, *ptboot, *plow, *pup, *ptbootpre, *ppilot, probi[nrow], cumprobi, sumprobi, unif, w[nrow][nrow], conf = asReal(Conf), hcontx, mterm, m, m2, sd, qterm, quant = qnorm(conf, 0, 1, TRUE, FALSE);
  SEXP probcureboot = PROTECT(allocVector(REALSXP, dummy));
  SEXP deltabootorder = PROTECT(allocVector(INTSXP, nrow));
  SEXP xbootorder = PROTECT(allocVector(REALSXP, nrow));
  SEXP tbootpre = PROTECT(allocVector(REALSXP, nrow));
  SEXP tbootorder = PROTECT(allocVector(REALSXP, nrow));
  SEXP x0 = PROTECT(allocVector(REALSXP, 1));
  SEXP lx0 = PROTECT(ScalarInteger(dummy));
  SEXP h0 = PROTECT(allocVector(REALSXP, 1));
  px = REAL(X);
  pdatat = REAL(Datat);
  pdatax = REAL(Datax);
  pdatadelta = INTEGER(Datadelta);
  ph = REAL(H);
  pxboot = REAL(xbootorder);
  ptbootpre = REAL(tbootpre);
  ptboot = REAL(tbootorder);
  pdeltaboot = INTEGER(deltabootorder);
  ppilot = REAL(Pilot);
  SEXP localdummy = PROTECT(allocVector(LGLSXP, 1));
  LOGICAL(localdummy)[0] = TRUE;
  SEXP onemdatadelta = PROTECT(allocVector(INTSXP, nrow));
  ponemdatadelta = INTEGER(onemdatadelta);
  for (i = 0; i < nrow; i++)
    ponemdatadelta[i] = 1 - pdatadelta[i];
  GetRNGstate();
  if (LOGICAL(Local)[0]) {
    SEXP low = PROTECT(allocVector(REALSXP, lx));
    plow = REAL(low);
    SEXP up = PROTECT(allocVector(REALSXP, lx));
    pup = REAL(up);
    SEXP confband = PROTECT(allocVector(VECSXP, 2));
    for (contx = 0; contx < lx; contx++) { // 1. for 
      REAL(x0)[0] = px[contx];
      REAL(h0)[0] = ph[contx];
      qterm = REAL(Q)[contx];
      m = 0;
      m2 = 0;
      hcontx = ppilot[contx];
      for (i = 0; i < nrow; i++) {
	pdataxi = pdatax[i];
	for (j = 0; j < nrow; j++) {
	  wijp = fabs(pdataxi - pdatax[j])/hcontx;
	w[i][j] = (wijp < 1.0) ?  1.0 - wijp*wijp : 0.0;
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
	probcureboot = probcurenp0(tbootorder, xbootorder, deltabootorder, Nrow, x0, lx0, h0, lx0, localdummy);
	mterm = REAL(probcureboot)[0];
	m += mterm;
	m2 += mterm*mterm;
	UNPROTECT(1);
      } // end of 2. for (in contb)
      sd = sqrt((m2 - m*m/b)/(b - 1));
      plow[contx] = fmax(fmin(qterm - quant*sd, 1), 0);
      pup[contx] = fmax(fmin(qterm + quant*sd, 1), 0);
    } // end of 1. for (in contx)
    PutRNGstate();
    SET_VECTOR_ELT(confband, 0, low);
    SET_VECTOR_ELT(confband, 1, up);
    UNPROTECT(13);
    return confband;
  }
  else { // Local == FALSE
    SEXP confbandlist = PROTECT(allocVector(VECSXP, lh));
    for (conth = 0; conth < lh; conth++) { // 1. for
      SEXP low = PROTECT(allocVector(REALSXP, lx));
      plow = REAL(low);
      SEXP up = PROTECT(allocVector(REALSXP, lx));
      pup = REAL(up);
      SEXP confband = PROTECT(allocVector(VECSXP, 2));
      REAL(h0)[0] = ph[conth];
      for (contx = 0; contx < lx; contx++) { // 2. for
	REAL(x0)[0] = px[contx];
	qterm = REAL(VECTOR_ELT(Q, conth))[contx];
	m = 0;
	m2 = 0;
	hcontx = ppilot[contx];
	for (i = 0; i < nrow; i++) {
	  pdataxi = pdatax[i];
	  for (j = 0; j < nrow; j++) {
	    wijp = fabs(pdataxi - pdatax[j])/hcontx;
	    w[i][j] = (wijp < 1.0) ?  1.0 - wijp*wijp : 0.0;
	  }
	}
	for (contb = 0; contb < b; contb++) { // 3. for
	  for (i = 0; i < nrow; i++) { // 4. for
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
	  } // end of 4. for (in i)
	  R_orderVector(torder, nrow, PROTECT(Rf_lang2(tbootpre, onemdatadelta)), TRUE, FALSE);
	  for (i = 0; i < nrow; i++) {
	    ptboot[i] = ptbootpre[torder[i]];
	    pxboot[i] = pdatax[torder[i]];
	    pdeltaboot[i] = deltaboot[torder[i]];
	  }
	  probcureboot = probcurenp0(tbootorder, xbootorder, deltabootorder, Nrow, x0, lx0, h0, lx0, localdummy);
	  mterm = REAL(probcureboot)[0];
	  m += mterm;
	  m2 += mterm*mterm;
	  UNPROTECT(1);
	} // end of 3. for (in contb)
	sd = sqrt((m2 - m*m/b)/(b - 1));
	plow[contx] = fmax(fmin(qterm - quant*sd, 1), 0);
	pup[contx] = fmax(fmin(qterm + quant*sd, 1),0);
      } // end of 2. for (in contx)
      SET_VECTOR_ELT(confband, 0, low);
      SET_VECTOR_ELT(confband, 1, up);
      SET_VECTOR_ELT(confbandlist, conth, confband);
      UNPROTECT(3);
    } // end of 1. for (in conth)
    PutRNGstate();
    UNPROTECT(11);
    return confbandlist;
  }
}
