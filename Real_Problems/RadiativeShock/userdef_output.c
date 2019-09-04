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
  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  DOM_LOOP(k, j, i) {

    NVAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];
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
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{
  Image *image;

  /* HDF5 output cannot be controlled yet. Everything is output.*/

  /* VTK output */
  SetOutputVar("prs",  VTK_OUTPUT, NO);
  SetOutputVar("te",   VTK_OUTPUT, YES);
  SetOutputVar("lmd",  VTK_OUTPUT, YES);


  /* FLT output */
  SetOutputVar("prs",  FLT_OUTPUT, NO);
  SetOutputVar("te",   FLT_OUTPUT, YES);
  SetOutputVar("lmd",  FLT_OUTPUT, YES);


#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}

