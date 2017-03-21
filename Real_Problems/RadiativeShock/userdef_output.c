#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k, nv;
  double ***te, ***lmd;
  double v[NVAR], rhs[NVAR];
  double *x1, *x2, *x3;

  /* New variables - names must exist under uservar */
  te = GetUserVar("te");
  lmd = GetUserVar("lmd");

  /* These are the geometrical central points */
  //x1 = grid[IDIR].x;
  //x2 = grid[JDIR].x;
  //x3 = grid[KDIR].x;

  /* These are the volumetric central points */
  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  DOM_LOOP(k, j, i) {

    for (nv = 0; nv < NVAR; nv++) v[nv] = d->Vc[nv][k][j][i];
    double mu = MeanMolecularWeight(v);

    /* Temperature */
    te[k][j][i] = v[PRS] / v[RHO] * mu;

    /* Cooling rate */
    v[RHOE] /= g_gamma - 1.;
    Radiat(v, rhs);
    lmd[k][j][i] = rhs[RHOE];

  }

}

/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{
  Image *image;

  /* HDF5 output cannot be controlled yet. Everything is output.*/

  /* VTK output */
  SetDumpVar("prs",  VTK_OUTPUT, NO);
  SetDumpVar("te",   VTK_OUTPUT, YES);
  SetDumpVar("lmd",  VTK_OUTPUT, YES);


  /* FLT output */
  SetDumpVar("prs",  FLT_OUTPUT, NO);
  SetDumpVar("te",   FLT_OUTPUT, YES);
  SetDumpVar("lmd",  FLT_OUTPUT, YES);

}

