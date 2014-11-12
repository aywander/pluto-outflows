#include "pluto.h"
/* AYW -- 2013-01-08 23:05 JST 
 * To get MeanMolecularWeight function */
#include "abundances.h"
/* -- AYW */
/* AYW -- 2014-05-24 10:36 JST
 * To use tabular MU */
#if MU_CALC == MU_TABLE
#include "rMuTable.h"
#endif
/* -- AYW */
#define frac_Z   1.e-3   /*   = N(Z) / N(H), fractional number density of metals (Z)
                                with respect to hydrogen (H) */ 
#define frac_He  0.082   /*   = N(Z) / N(H), fractional number density of helium (He)
                                with respect to hydrogen (H) */ 
#define A_Z      30.0    /*   mean atomic weight of heavy elements  */
#define A_He     4.004   /*   atomic weight of Helium  */
#define A_H      1.008   /*   atomic weight of Hydrogen  */

/* ***************************************************************** */
void Radiat (real *v, real *rhs)
/*
 *
 * NAME
 *
 *   Radiat
 *
 *
 * PURPOSE
 *
 *   Provide r.h.s. for tabulated cooling.
 * 
 *
 ******************************************************************* */
{
  int    klo, khi, kmid;
  real   mu, T, Tmid, scrh, dT;
  static int ntab;
  static real *L_tab, *T_tab, E_cost;
  
  FILE *fcool;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (T_tab == NULL){
    print1 (" > Reading table from disk...\n");
    fcool = fopen("cooltable.dat","r");
    if (fcool == NULL){
      print1 ("! cooltable.dat does not exists\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, 
                                       L_tab + ntab)!=EOF) {
      ntab++;
    }
    E_cost    = g_unitLength/g_unitDensity/pow(g_unitVelocity, 3.0);
  }

/* ---------------------------------------------
            Get temperature 
   --------------------------------------------- */

  if (v[PRS] < 0.0) v[PRS] = g_smallPressure;
  mu  = MeanMolecularWeight(v);
  T   = v[PRS]/v[RHO]*KELVIN*mu;

  if (T != T){
    printf (" ! Nan found in radiat \n");
    printf (" ! rho = %12.6e, pr = %12.6e\n",v[RHO], v[PRS]);
    QUIT_PLUTO(1);
  }

  if (T < g_minCoolingTemp) { 
    rhs[PRS] = 0.0;
    return;
  }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (T > T_tab[khi] || T < T_tab[klo]){
    print (" ! T out of range   %12.6e\n",T);
    QUIT_PLUTO(1);
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = T_tab[kmid];
    if (T <= Tmid){
      khi = kmid;
    }else if (T > Tmid){
      klo = kmid;
    }
  }

  dT      = T_tab[khi] - T_tab[klo];
  scrh    = L_tab[klo]*(T_tab[khi] - T)/dT + L_tab[khi]*(T - T_tab[klo])/dT;
  rhs[PRS] = -(g_gamma - 1.0)*scrh*v[RHO]*v[RHO];
  rhs[PRS] *= E_cost*g_unitDensity*g_unitDensity/(CONST_mp*CONST_mp);
  
}
#undef T_MIN

/* ******************************************************************* */
double MeanMolecularWeight (real *V)
/*
 *
 *
 *
 ********************************************************************* */
{
  /* AYW -- 2012-10-02 15:28 JST
   * Note this file is just for tabulated cooling!
   * When no cooling is available there is a MeanMolecularWeight
   * function in abundances.c.
   * */

#if MU_CALC == MU_TABLE

  int il;
  double por, r1, r2, frac;
  double y0, y1, y2, y3;

  if (mu_por == NULL){
    rMuTable()
  }

  /* Value of T/mu in cgs */
  por = V[PRS]/V[RHO]*KELVIN

  /* Find cell left index to interpolate at */
  il = hunter(mu_rad, mu_ndata, por);

  /* Linear fractional location of por in cell */
  r1 = mu_por[il];
  r2 = mu_por[il+1];
  frac = (r - r1)/(r2 - r1);

  /* Mu interpolation */
  y0 = mu_mu[il-1];
  y1 = mu_mu[il];
  y2 = mu_mu[il+1];
  y3 = mu_mu[il+2];
  return CubicCatmullRomInterpolate(y0, y1, y2, y3, frac);

#elif MU_CALC == MU_FRACTIONS

  return  ( (A_H + frac_He*A_He + frac_Z*A_Z) /
            (2.0 + frac_He + 2.0*frac_Z - 0.0));

#else

  return (MU_NORM);

#endif

  /* --AYW */

}


