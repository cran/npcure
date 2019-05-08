#ifndef NPCURE_H
#define NPCURE_H

#include <math.h>
#include <stdlib.h>
#include <string.h> // for memset
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdbool.h> // for bool type

SEXP berannp0(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh, SEXP Local, SEXP T, SEXP Lt);

SEXP berannp0confband(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh, SEXP Pilot, SEXP Probcurepilot, SEXP Conf, SEXP B, SEXP S, SEXP Local);

SEXP berannp0cv(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh);

SEXP latencynp0(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh, SEXP Local, SEXP T, SEXP Lt);

SEXP latencynp0hboot(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh, SEXP Pilot, SEXP Pilotprobcure, SEXP Latencypilot, SEXP B, SEXP Tmax);

SEXP latencynp0confband(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh, SEXP Pilot, SEXP Probcurepilot, SEXP Conf, SEXP B, SEXP S, SEXP Local);
 
SEXP probcurenp0(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh, SEXP Local);

SEXP probcurenp0hboot(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh, SEXP Pilot, SEXP Probcurepilot, SEXP B);

SEXP probcurenp0confband(SEXP Datat, SEXP Datax, SEXP Datadelta, SEXP Nrow, SEXP X, SEXP Lx, SEXP H, SEXP Lh, SEXP Conf, SEXP B, SEXP Pilot, SEXP Q, SEXP Local);

#endif


