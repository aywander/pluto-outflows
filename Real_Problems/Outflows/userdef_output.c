#include "pluto.h"
#include "pluto_usr.h"
#include "nrEOS.h"
#include "abundances.h"

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
  double ***te, ***spd;
  double ***prs, ***rho, ***vx1, ***vx2, ***vx3, dummy[NVAR];
  double mu, sp1, sp2, sp3;
  double *x1, *x2, *x3;
#if USE_FOUR_VELOCITY == YES
  double ***v1, ***v2, ***v3;
  double vel, speed, lorentz;
#endif

  /* New variables - names must exist under uservar */
  te   = GetUserVar("te");
  spd  = GetUserVar("spd");

/* Change to v instead of u = lorentz v */
#if USE_FOUR_VELOCITY == YES
  v1  = GetUserVar("v1");
  v2  = GetUserVar("v2");
  v3  = GetUserVar("v3");
#endif

  /* State variables */
  rho = d->Vc[RHO];
  prs = d->Vc[PRS];
  vx1 = d->Vc[VX1];
  vx2 = d->Vc[VX2];
  vx3 = d->Vc[VX3];

  /* These are the geometrical central points */
  //x1 = grid[IDIR].x;
  //x2 = grid[JDIR].x;
  //x3 = grid[KDIR].x;

  /* These are the volumetric central points */
  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;
  
  DOM_LOOP(k,j,i){

    /* Temperature */
    for (nv = 0; nv < NVAR; nv++) dummy[nv] = d->Vc[nv][k][j][i];
    mu = MeanMolecularWeight(dummy);
    te[k][j][i] = TempNrEOS(rho[k][j][i], prs[k][j][i], mu);

    /* Speed */
    sp1 = vx1[k][j][i]; 
    sp2 = vx2[k][j][i];
    sp3 = vx3[k][j][i];
    spd[k][j][i] = VMAG(x1[i], x2[j], x3[k], sp1, sp2, sp3);

#if USE_FOUR_VELOCITY == YES
    /* spd at this point is gamma * v. Solve for v. Then get gamma. 
     * c=1 if USE_FOUR_VELOCITY = YES. 
     * speed = Gamma * vel */
    speed = spd[k][j][i];
    vel = 1./(1./(speed*speed) + 1.);
    lorentz = speed/vel;

    spd[k][j][i] = vel;
    EXPAND(v1[k][j][i] = sp1/lorentz;,
           v2[k][j][i] = sp2/lorentz;,
           v3[k][j][i] = sp3/lorentz;);
#endif

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
#if USE_FOUR_VELOCITY
  EXPAND(SetDumpVar("v1",  VTK_OUTPUT, YES);,
         SetDumpVar("v2",  VTK_OUTPUT, YES);,
         SetDumpVar("v3",  VTK_OUTPUT, YES););
  EXPAND(SetDumpVar("vx1",  VTK_OUTPUT, YES);,
         SetDumpVar("vx2",  VTK_OUTPUT, YES);,
         SetDumpVar("vx3",  VTK_OUTPUT, YES););
#else
  EXPAND(SetDumpVar("v1",  VTK_OUTPUT, NO);,
         SetDumpVar("v2",  VTK_OUTPUT, NO);,
         SetDumpVar("v3",  VTK_OUTPUT, NO););
  EXPAND(SetDumpVar("vx1",  VTK_OUTPUT, YES);,
         SetDumpVar("vx2",  VTK_OUTPUT, YES);,
         SetDumpVar("vx3",  VTK_OUTPUT, YES););
#endif
  SetDumpVar("prs",  VTK_OUTPUT, YES);
  SetDumpVar("tr1",  VTK_OUTPUT, YES);
  SetDumpVar("tr2",  VTK_OUTPUT, YES);
  SetDumpVar("te",  VTK_OUTPUT, YES);
  SetDumpVar("spd",  VTK_OUTPUT, YES);


  /* FLT output */
#if USE_FOUR_VELOCITY
  EXPAND(SetDumpVar("v1",  FLT_OUTPUT, YES);,
         SetDumpVar("v2",  FLT_OUTPUT, YES);,
         SetDumpVar("v3",  FLT_OUTPUT, YES););
  EXPAND(SetDumpVar("vx1",  FLT_OUTPUT, YES);,
         SetDumpVar("vx2",  FLT_OUTPUT, YES);,
         SetDumpVar("vx3",  FLT_OUTPUT, YES););
#else
  EXPAND(SetDumpVar("v1",  FLT_OUTPUT, NO);,
         SetDumpVar("v2",  FLT_OUTPUT, NO);,
         SetDumpVar("v3",  FLT_OUTPUT, NO););
  EXPAND(SetDumpVar("vx1",  FLT_OUTPUT, YES);,
         SetDumpVar("vx2",  FLT_OUTPUT, YES);,
         SetDumpVar("vx3",  FLT_OUTPUT, YES););
#endif
  SetDumpVar("prs",  FLT_OUTPUT, YES);
  SetDumpVar("tr1",  FLT_OUTPUT, YES);
  SetDumpVar("tr2",  FLT_OUTPUT, YES);
  SetDumpVar("te",  FLT_OUTPUT, YES);
  SetDumpVar("spd",  FLT_OUTPUT, YES);


  /* PNG output */
#if USE_FOUR_VELOCITY
  EXPAND(SetDumpVar("v1",  PNG_OUTPUT, YES);,
         SetDumpVar("v2",  PNG_OUTPUT, YES);,
         SetDumpVar("v3",  PNG_OUTPUT, YES););
  EXPAND(SetDumpVar("vx1",  PNG_OUTPUT, YES);,
         SetDumpVar("vx2",  PNG_OUTPUT, YES);,
         SetDumpVar("vx3",  PNG_OUTPUT, YES););
#else
  EXPAND(SetDumpVar("v1",  PNG_OUTPUT, NO);,
         SetDumpVar("v2",  PNG_OUTPUT, NO);,
         SetDumpVar("v3",  PNG_OUTPUT, NO););
  EXPAND(SetDumpVar("vx1",  PNG_OUTPUT, YES);,
         SetDumpVar("vx2",  PNG_OUTPUT, YES);,
         SetDumpVar("vx3",  PNG_OUTPUT, YES););
#endif
  SetDumpVar("prs",  PNG_OUTPUT, YES);
  SetDumpVar("tr1",  PNG_OUTPUT, YES);
  SetDumpVar("tr2",  PNG_OUTPUT, YES);
  SetDumpVar("te",  PNG_OUTPUT, YES);
  SetDumpVar("spd",  PNG_OUTPUT, YES);



  /* density slice */
  image = GetImage("rho");
#if COMPONENTS > 2
  image->slice_plane = X13_PLANE;
  image->slice_coord = 0.0;
#endif
  image->max = image->min = 0.;
  image->logscale = 1;

  /* density slice */
  image = GetImage("prs");
#if COMPONENTS > 2
  image->slice_plane = X13_PLANE;
  image->slice_coord = 0.0;
#endif
  image->max = image->min = 0.;
  image->logscale = 1;


  /* first tracer slice */
#if NTRACER > 0
  image = GetImage("tr1");
#if COMPONENTS > 2
  image->slice_plane = X13_PLANE;
  image->slice_coord = 0.0;
#endif
  image->max = 1.;
  image->min = 0.;
  image->logscale = 1;
#endif

  /* second tracer slice */
#if NTRACER > 1
  image = GetImage("tr2");
#if COMPONENTS > 2
  image->slice_plane = X13_PLANE;
  image->slice_coord = 0.0;
#endif
  image->max = 1.;
  image->min = 0.;
  image->logscale = 1;
#endif


}

