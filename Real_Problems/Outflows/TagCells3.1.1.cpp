#include <cstdio>
#include <string>
#include <functional>
#include <numeric>
/* AYW -- 2012-06-14 10:42 JST
 * to get max function */
#include <algorithm>
/* -- AYW */
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"
#include "pluto_usr.h"


/* A typesafe sgn function */
template <typename T> int sgn(T val) { return (val > T(0)) - (val < T(0)); }


/* A typesafe heaviside function */
template <typename T> int heaviside(T val) { return !(val > T(0)); }


/* ************************************************************************* */
void PatchPluto::GETRELGRAD(FArrayBox& gFab, FArrayBox& UFab, const Box& b, int level, 
    int maxLevel, std::vector<int> refRatios)
/*
 *
 * PURPOSE
 *
 *   Tag cells for refinement by computing grad[k][j][i]. 
 *   By default a convex combination of the first and second
 *   derivative of the total energy density is used.
 *   alpha = 0 --> triggers refinement towards the 2nd derivative
 *   alpha = 1 --> triggers refinement towards the 1st derivative
 *
 * 
 *
 *************************************************************************** */
{
  CH_assert(m_isDefined);

  int nv, i, j, k;
  int Uib, Uie, Ujb=0, Uje=0, Ukb=0, Uke=0;
  int Gib, Gie, Gjb=0, Gje=0, Gkb=0, Gke=0;
  int i1, j1, k1; 

  double rp, rm, r;
  double x, dqx_p, dqx_m, d2qx, den_x;
  double y, dqy_p, dqy_m, d2qy, den_y;
  double z, dqz_p, dqz_m, d2qz, den_z;

  double alpha, qref, gr1, gr2;

  double eps = 0.01;
  double ***UU[NVAR], ***q, ***grad;

  //double us[NVAR], vs[NVAR], mu;
  //static double **T;

  /* AYW --
   * 2011-04-21 19:19 JST 
   * My new local variables */
  double cyl_radius2, radius2, xlo, ylo, zlo, xhi, yhi, zhi;
  double boxdist2, xcomp, ycomp, zcomp;
  int iref;
  float refRatioProd;
  float cxc, cyc, czc, csz;
  float cxl, cxu, cyl, cyu, czl, czu;
  float xel, yel, zel;
  float large = 1.e37;
  double tr_small = 1.e-12;
  double ***qj, ***qc, ***qn;
  double pj, pc;
  /* -- AYW */


  rp = rm = r = 1.0;

/* -----------------------------------------------
    The solution array U is defined on the box 
    [Uib, Uie] x [Ujb, Uje] x [Ukb, Uke], which 
    differs from that of gFab ([Gib,...Gke]), 
    typically one point larger in each direction. 
   ----------------------------------------------- */

  D_EXPAND(Uib = UFab.loVect()[IDIR]; Uie = UFab.hiVect()[IDIR]; ,
           Ujb = UFab.loVect()[JDIR]; Uje = UFab.hiVect()[JDIR]; ,
           Ukb = UFab.loVect()[KDIR]; Uke = UFab.hiVect()[KDIR]; );

  D_EXPAND(Gib = gFab.loVect()[IDIR]; Gie = gFab.hiVect()[IDIR]; ,
           Gjb = gFab.loVect()[JDIR]; Gje = gFab.hiVect()[JDIR]; ,
           Gkb = gFab.loVect()[KDIR]; Gke = gFab.hiVect()[KDIR]; );

  for (nv=0 ; nv<NVAR ; nv ++){
    UU[nv] = chmatrix(Ukb, Uke, Ujb, Uje, Uib, Uie, UFab.dataPtr(nv));
  }
  grad = chmatrix(Gkb, Gke, Gjb, Gje, Gib, Gie,gFab.dataPtr(0));

/* -- build temperature array -- */
/*  
  if (T == NULL) T = Array_2D(NMAX_POINT, NMAX_POINT, double);

  for (k = Ukb; k <= Uke; k++) {
  for (j = Ujb; j <= Uje; j++) {
  for (i = Uib; i <= Uie; i++) {

    #if GEOMETRY == CYLINDRICAL 
     r = m_dx*(i+0.5);
    #endif
    for (nv = 0; nv < NVAR; nv++) us[nv] = UU[nv][k][j][i]/r;

    alpha = 1.0/us[DN];
    for (nv = NFLX; nv < NVAR; nv++) vs[nv] = us[nv]*alpha;
    #if INCLUDE_COOLING == NO
     mu = 1.237190; 
    #else
     mu = MEAN_MOLECULAR_WEIGHT(vs);
    #endif

    #if AMR_EN_SWITCH == YES
     d2q = 0.0;

//    T[j-Ujb][i-Uib] = us[EN]/pow(alpha,gmm-2.0)*KELVIN*mu;
    T[j-Ujb][i-Uib] = us[EN]/pow(alpha,gmm-2.0);    this is T = p/rho  
    #else
     d2q  = EXPAND(us[MX]*us[MX], + us[MY]*us[MY], + us[MZ]*us[MZ]);
     d2q /= us[DN];
     #if PHYSICS == MHD
      d2q += EXPAND(us[BX]*us[BX], + us[BY]*us[BY], + us[BZ]*us[BZ]);
     #endif
//    T[j-Ujb][i-Uib] = (gmm - 1.0)*(us[EN] - 0.5*d2q)/us[DN]*KELVIN*mu;
    T[j-Ujb][i-Uib] = us[EN] - 0.5*d2q;  // internal energy 
    #endif
  }}}
*/
/* -----------------------------------------------
    the parameter alpha controls the bias towards
    1st derivative criterion (alpha = 1) or 2nd
    derivative (alpha = 0)
   ----------------------------------------------- */


  alpha = 0.0;
  #if (EOS != ISOTHERMAL) && (AMR_EN_SWITCH == NO)
   q  = UU[EN];   
   /* AYW -- 2012-06-14 09:32 JST 
    * Also get conserved quantity of continuity equation */
   qn = UU[DN];   
  /* -- AYW */

  #else
   q  = UU[DN];
   /* AYW -- 2012-06-14 09:32 JST 
    * Also get conserved quantity of continuity equation */
   qn = q;
  /* -- AYW */
  #endif

  /* AYW -- 2012-06-13 19:51 JST
   * Wind tracer */
  qj  = UU[TR];

  /* Clouds tracer */
  #if CLOUDS != NONE
  qc  = UU[TR+1];
  #endif
  /* -- AYW */

  /* AYW -- 2012-06-13 20:07 JST 
   * Level check 
   * The first time I implemented stuff here (PLUTO 3.1.0) 
   * LEVEL didn't exist as a global variable yet.
   * LEVEL does not give the same as level. Dunno why yet. 
   * Note, LEVEL is set in LevelPluto.cpp. */
  //printf("level = %d, LEVEL = %d\n", level, LEVEL);
  //CH_assert(level == LEVEL);
  /* AYW */

  for (k = Gkb; k <= Gke; k++) { z = k*m_dx + DOM_XBEG[KDIR];
  for (j = Gjb; j <= Gje; j++) { y = j*m_dx + DOM_XBEG[JDIR];
  for (i = Gib; i <= Gie; i++) { x = i*m_dx + DOM_XBEG[IDIR];

   /* AYW -- 2011-04-21 19:19 JST
    * Always refine near wind base
    * */

    //radius2 =  z*z + y*y + x*x;
    cyl_radius2 =  z*z + y*y;

    // Primitive advected scalar
    pj = qj[k][j][i]/max(qn[k][j][i], tr_small);


    if ( ( cyl_radius2 <= aux[URAD]+aux[UTHK] ) && 
        ( x <= aux[UTHK] ) && 
        ( pj > 0.8) ) {
      grad[k][j][i] = large;
    }

    /* While wind is refined at the highest level given,
     * other regions are refined up to aux[LEV1] */
    else if (level < aux[LEV1]){

      #if CLOUDS == SINGLE_CUBE

      pc = qc[k][j][i]/max(qn[k][j][i], tr_small);

      if ( pc > 0.95 ) {
        grad[k][j][i] = large;
      }

      else

      #endif

      {

      /* Normal refinement criterion from here. 
       * NOTE: normal refinement criterium only up to LEV1. */
      /* -- AYW */


      //i1 = i-Uib; j1 = j-Ujb; k1 = k-Ukb; /* -- use these index for a newly 
      //                                         generated array like T      -- */

      #if GEOMETRY == CYLINDRICAL
       rp = (i+0.5)/(i+1.5);
       rm = (i+0.5)/(i-0.5);
      #endif
  /*
      D_EXPAND(dpx = T[j1][i1+1]*rp - T[j1][i1];
               dmx = T[j1][i1-1]*rm - T[j1][i1];  ,
               dpy = T[j1+1][i1] - T[j1][i1];
               dmy = T[j1-1][i1] - T[j1][i1];     ,
               dpz = q[k+1][j][i] - q[k][j][i];
               dmz = q[k-1][j][i] - q[k][j][i];)

      qref = fabs(T[j1][i1]);
  */

      D_EXPAND(dqx_p = q[k][j][i+1]*rp - q[k][j][i];
               dqx_m = q[k][j][i-1]*rm - q[k][j][i];  ,
               dqy_p = q[k][j+1][i] - q[k][j][i];
               dqy_m = q[k][j-1][i] - q[k][j][i];     ,
               dqz_p = q[k+1][j][i] - q[k][j][i];
               dqz_m = q[k-1][j][i] - q[k][j][i];)

    /* ---------------------------------------------------
                      New version
       --------------------------------------------------- */

      /*
      D_EXPAND(d2qx = dqx_p + dqx_m;  ,
               d2qy = dqy_p + dqy_m;  ,
               d2qz = dqz_p + dqz_m;)

      D_EXPAND(
        den_x = 2.0*fabs(q[k][j][i]) + fabs(q[k][j][i+1]) + fabs(q[k][j][i-1]);
        den_x = fabs(dqx_p) + fabs(dqx_m) + eps*den_x;    ,

        den_y = 2.0*fabs(q[k][j][i]) + fabs(q[k][j+1][i]) + fabs(q[k][j-1][i]);
        den_y = fabs(dqy_p) + fabs(dqy_m) + eps*den_y;    ,

        den_z = 2.0*fabs(q[k][j][i]) + fabs(q[k+1][j][i]) + fabs(q[k-1][j][i]);
        den_z = fabs(dqz_p) + fabs(dqz_m) + eps*den_z;
      )

      gr2  = D_EXPAND(d2qx*d2qx, + d2qy*d2qy, + d2qz*d2qz);
      gr2 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);
      grad[k][j][i] = sqrt(gr2);

      
      */

    /* ---------------------------------------------------
                        Old version
       --------------------------------------------------- */


   // -- first derivative -- 

      eps   = 0.05;
      qref = fabs(q[k][j][i]); 
      gr1 = D_EXPAND( fabs(dqx_p - dqx_m)/qref ,
                    + fabs(dqy_p - dqy_m)/qref ,
                    + fabs(dqz_p - dqz_m)/qref);

   // -- second derivative --

      gr2 = D_EXPAND( fabs(dqx_p + dqx_m)/(fabs(dqx_p) + fabs(dqx_m) + eps*qref) ,
                    + fabs(dqy_p + dqy_m)/(fabs(dqy_p) + fabs(dqy_m) + eps*qref) ,
                    + fabs(dqz_p + dqz_m)/(fabs(dqz_p) + fabs(dqz_m) + eps*qref));


      /* AYW -- 
       * 2011-11-08 14:02 JST 
       * Different threshold for RHD. Not sure why, but 
       * U[EN] must be quite different. Future chekcs are needed */

#if PHYSICS == RHD || PHYSICS == RMHD

      grad[k][j][i] = 0.125*(alpha*gr1 + (1.0 - alpha)*gr2);
#else
      /* -- AYW */

      grad[k][j][i] = alpha*gr1 + (1.0 - alpha)*gr2;

      /* AYW -- 
       * 2011-11-08 14:02 JST  */
#endif
      /* -- AYW */



      
    /* --------------------------------------------------- */

    /* AYW -- 
     * 2011-04-21 19:20 JST */
      }
    }

    /* No refinement if we're at or above level aux[LEV1] */
    else {
      grad[k][j][i] = 0.;
    }
    /* -- AYW */

  }}}

  for (nv=0 ; nv<NVAR ; nv ++){
    free_chmatrix(UU[nv], Ukb, Uke, Ujb, Uje, Uib, Uie);
  }
  free_chmatrix(grad, Gkb, Gke, Gjb, Gje, Gib, Gie);
}


