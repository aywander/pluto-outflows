#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pluto.h"
#include "pluto_usr.h"
#include "interpolation.h"

//TODO: Use PLUTO interpolation routines

int hunter(const double arr[], const int narr, const double val){

  /* 
   * This function returns index, il, of an array, arr, such that
   * arr[il] < val < arr[il + 1]
   */

  int il, ir, shift;
  double arrb, arre;
  double vl, vr;
  double ldelta, rdelta, delta, dprod;


  /* Beginning and end values */
  arrb = arr[0];
  arre = arr[narr-1];

  /* Bounds check */
  if (val < arrb){
    print1("Error: interpolation.c: hunter: interpolation out of (lower) bounds.\n");
    print1("val :  %e\n", val);
    print1("arrb:  %e\n", arrb);
    exit(1);
  }
  else if (val > arre){
    print1("Error: interpolation.c: hunter: interpolation out of (upper) bounds.\n");
    print1("val  %e\n", val);
    print1("arre %e\n", arre);
    exit(1);
  }

  /* Initial linear guess, left and right indices */
  il = (int) (narr - 1)*(val - arrb)/(arre - arrb);
  ir = il + 1;

	/* Left and right values from arrays */
  vl = arr[il];
  vr = arr[ir];

	/* Partial and total differences in r */
	ldelta = val - vl;
	rdelta = vr - val;
  delta = vr - vl;

  /* If positive, then we're in the right interval */
  dprod = ldelta*rdelta;

  while (dprod < 0){
    /* Calculate single cell shift direction */
    shift = SGN(ldelta);

    /* New indices */
    il += shift; 
    ir = il + 1;

    /* Left and right values from arrays */
    vl = arr[il];
    vr = arr[ir];

    /* Partial and total differences in r */
    ldelta = val - vl;
    rdelta = vr - val;
    delta = vr - vl;

    /* If positive, then we're in the right interval */
    dprod = ldelta*rdelta;
  }

  return il;

}

double LinearInterpolate(double y1, double y2, double fc){
  /* Linear interpolator */

  return ( y1*(1 - fc) + y2*fc );
}


double CosineInterpolate(double y1, double y2, double fc){
  /* Cosine interpolator */

  double fc2;
  double pi = 3.14159265358979323846;
  fc2 = (1. - cos(fc*pi))/2.;
  return ( y1*(1. - fc2) + y2*fc2 );
}


double CubicInterpolate(double y0, double y1, 
                        double y2, double y3, 
                        double fc){
  /* Cubic interpolator */

  double a0, a1, a2, a3, fc2;

  fc2 = fc*fc;
  a0 = y3 - y2 - y0 + y1;
  a1 = y0 - y1 - a0;
  a2 = y2 - y0;
  a3 = y1;

  return ( a0*fc*fc2 + a1*fc2 + a2*fc + a3 );
}


double CubicCatmullRomInterpolate(double y0, double y1, 
                                  double y2, double y3, 
                                  double fc){
  /* Cubic interpolator Catmull-Rom Spline */

  double a0, a1, a2, a3, fc2;

  fc2 = fc*fc;
  a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3;
  a1 = y0 - 2.5*y1 + 2*y2 - 0.5*y3;
  a2 = -0.5*y0 + 0.5*y2;
  a3 = y1;

  return ( a0*fc*fc2 + a1*fc2 + a2*fc + a3 );
}


double HermiteInterpolate(double y0, double y1,
                          double y2, double y3,
                          double fc, double tension, 
                          double bias){
/* Hermite interpolator
 * Tension: 1 is high, 0 normal, -1 is low
 * Bias: 0 is even,
 *       positive is towards first segment,
 *       negative towards the other
 */


  double m0, m1, fc2, fc3;
  double a0, a1, a2, a3;

  fc2 = fc*fc;
  fc3 = fc2*fc;
  m0  = (y1 - y0)*(1. + bias)*(1. - tension)/2.;
  m0 += (y2 - y1)*(1. - bias)*(1. - tension)/2.;
  m1  = (y2 - y1)*(1. + bias)*(1. - tension)/2.;
  m1 += (y3 - y2)*(1. - bias)*(1. - tension)/2.;
  a0 =  2.*fc3 - 3.*fc2 + 1.;
  a1 =     fc3 - 2.*fc2 + fc;
  a2 =     fc3 -    fc2;
  a3 = -2.*fc3 + 3.*fc2;

  return ( a0*y1 + a1*m0 + a2*m1 + a3*y2 );
}


/* ************************************************************* */
void x3u_3d_extrapol(double ***a, int kb, int i, int j, int k, Grid *grid)
/*
 *
 * Quadratic extrapolation in the upper boundary of X3
 *
 *************************************************************** */
{
  double *x1, *x2, *x3;
  double y0, y1, y2, z0, z1, z2;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  y0 = a[kb][j][i];
  y1 = a[kb-1][j][i];
  y2 = a[kb-2][j][i];

  z0 = x3[kb];
  z1 = x3[kb-1];
  z2 = x3[kb-2];

  a[k][j][i] = y0 + (y1 - y0) / (z1 - z0) * (x3[k] - z0) + 
               ((y2 - y1) / (z2 - z1) / (z2 - z0) - 
                (y1 - y0) / (z1 - z0) / (z2 - z0)) * 
               (x3[k] - z0) * (x3[k] - z1);

  }



/* ************************************************************* */
void x2l_3d_extrapol (double ***a, int jb, int i, int j, int k, Grid *grid)
/*
 *
 * Quadratic extrapolation in the lower boundary of X2
 *
 *************************************************************** */
{
  double *x1, *x2, *x3;
  double y0, y1, y2, z0, z1, z2;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  y0 = a[k][jb][i];
  y1 = a[k][jb+1][i];
  y2 = a[k][jb+2][i];

  z0 = x2[jb];
  z1 = x2[jb+1];
  z2 = x2[jb+2];

  a[k][j][i] = y0 + (y1 - y0) / (z1 - z0) * (x3[k] - z0) + 
               ((y2 - y1) / (z2 - z1) / (z2 - z0) - 
                (y1 - y0) / (z1 - z0) / (z2 - z0)) * 
               (x3[k] - z0) * (x3[k] - z1);

}

