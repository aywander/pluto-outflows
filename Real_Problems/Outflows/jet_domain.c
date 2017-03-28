/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Adjust the effective domain size when the -xNjet option is given.

  The SetJetDomain() function shortens the domain integration index in the 
  direction specified by either '\c -x1jet', '\c -x2jet' or '\c -x3jet' to save
  computational time.
  This is done by setting the final integration index in the direction
  of jet propagation to that of the first zone (counting from the top) with 
  non-zero pressure gradient (GetRightmostIndex()).
  The function UnsetJetDomain() restores the initial computational domain
  size.
  Useful for problems involving jet propagation.

  \note In parallel, the domain is \e not decomposed along the
        propagation direction (see ParseCmdLineArgs()).
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Dec 24, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"
/* AYW -- 2012-12-27 17:01 JST 
 * to use JD_CONST and JD_MODE */
#include "pluto_usr.h"
/*-- AYW */

static int jd_nbeg, jd_nend, jd_npt, jd_ntot, rbound;

/* AYW -- 2012-11-15 16:05 JST */
//static int GetRightmostIndex(int, double ***);
static int GetRightmostIndex(int dir, double ***q, double ***rho, Grid *grid);
/* -- AYW */

/* ********************************************************************* */
void SetJetDomain (const Data *d, int dir, int log_freq, Grid *grid)
/*!
 *
 * Adjust the size of computational domain by reducing the final
 * computationa index in the direction of propagation.
 *
 * \param [in]         d   pointer to Data structure
 * \param [in]       dir   the direction of propagation
 * \param [in]  log_freq   the output log frequency
 * \param [in,out]  grid   pointer to array of Grid structures
 * 
 *********************************************************************** */
{
  int i, j, k, ngh;
  int n, n_glob;
  static int first_call = 1;
  double ***pr, ***dn, dp;

  dn = d->Vc[RHO];
  #if HAVE_ENERGY
   pr = d->Vc[PRS];
  #else
   pr = d->Vc[RHO];
  #endif

  ngh = grid[dir].nghost;

/* -- save original domain offsets 
      and return for the first time -- */

  if (first_call){
    if (dir == IDIR) {
      jd_nbeg = IBEG; jd_nend = IEND;
      jd_npt  = NX1;  jd_ntot = NX1_TOT;
    }else if (dir == JDIR){
      jd_nbeg = JBEG; jd_nend = JEND;
      jd_npt  = NX2;  jd_ntot = NX2_TOT;
    }else if (dir == KDIR){
      jd_nbeg = KBEG; jd_nend = KEND;
      jd_npt  = NX3;  jd_ntot = NX3_TOT;
    }

    rbound = grid[dir].rbound;
    first_call = 0;
    return;
  }

/* --------------------------------------
    find where grad(p) is not zero and 
    add safety guard cells 
   -------------------------------------- */

  /* AYW -- 2012-11-15 16:07 JST */
  //n = GetRightmostIndex(dir, pr) + 2*ngh;
  n = GetRightmostIndex(dir, pr, dn, grid) + 2*ngh;
  /* -- AYW */

  if (n < jd_nbeg + ngh) n = jd_nbeg + ngh;
  if (n > jd_nend)       n = jd_nend;

  #ifdef PARALLEL
   MPI_Allreduce (&n, &n_glob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   n = n_glob;
  #endif

  if (g_stepNumber%log_freq==0){
/*    print1 ("- SetJetDomain: index %d / %d\n",n,jd_nend); */
  }

/* ------------------------------------------------------------------
    Change global integer variables giving the rightmost integration 
    indexes. Also, suppress the rightmost physics boundary (rbound) 
    if solution has not yet reached it. 
   ------------------------------------------------------------------ */

  if (dir == IDIR) {
    IEND = grid[IDIR].lend = n;
    NX1_TOT = IEND + ngh;
    NX1     = IEND - ngh;

    grid[IDIR].rbound = rbound;
    if (jd_nend != IEND) grid[IDIR].rbound = 0;
    
  } else if (dir == JDIR){

    JEND = grid[JDIR].lend = n;
    NX2_TOT = JEND + ngh;
    NX2     = JEND - ngh;
    grid[JDIR].rbound = rbound;
    if (jd_nend != JEND) grid[JDIR].rbound = 0;

  }else if (dir == KDIR){

    KEND = grid[KDIR].lend = n;
    NX3_TOT = KEND + ngh;
    NX3     = KEND - ngh;

    grid[KDIR].rbound = rbound;
    if (jd_nend != KEND) grid[KDIR].rbound = 0;
  }
  
/* ------------------------------------------------------
    Recompute RBox(es) after indices have been changed
   ------------------------------------------------------ */

  SetRBox();

/*
  print1 ("Box = %d/%d\n", n, jd_nend);
  {
    double t0, t1, i0, i1;
    FILE *fp;
    t0 = g_time;
    t1 = t0 + g_dt;
    i0 = (int)(t0/0.5);
    i1 = (int)(t1/0.5);

    if (i0 != i1){
      fp = fopen("box.out","w");
      fprintf (fp,"%f  %f  %f  %f\n",0.0, grid[IDIR].x[IEND], 
                                     0.0, grid[JDIR].x[JEND]);
      fclose(fp);
    } 
  }
*/
}

/* ********************************************************************* */
void UnsetJetDomain (const Data *d, int dir, Grid *grid)
/*!
 *  Restore original (full domain) indexes.
 *
 *********************************************************************** */
{

  if (dir == IDIR){
    IBEG = grid[IDIR].lbeg = jd_nbeg;
    IEND = grid[IDIR].lend = jd_nend;
    NX1_TOT = jd_ntot;
    NX1     = jd_npt;
    grid[IDIR].rbound = rbound;
  }

  if (dir == JDIR){
    JBEG = grid[JDIR].lbeg = jd_nbeg;
    JEND = grid[JDIR].lend = jd_nend;
    NX2_TOT = jd_ntot;
    NX2     = jd_npt;
    grid[JDIR].rbound = rbound;
  }

  if (dir == KDIR){
    KBEG = grid[KDIR].lbeg = jd_nbeg;
    KEND = grid[KDIR].lend = jd_nend;
    NX3_TOT = jd_ntot;
    NX3     = jd_npt;
    grid[KDIR].rbound = rbound;
  }

/* ------------------------------------------------------
    Recompute RBox(es) after indices have been changed
   ------------------------------------------------------ */

  SetRBox();

}

/* ********************************************************************* */
/* AYW -- 2012-11-15 16:09 JST */
//int GetRightmostIndex (int dir, double ***q)
int GetRightmostIndex (int dir, double ***q, double ***rho, Grid *grid)
/* -- AYW  */
/*!
 *
 * Find the local index j where grad(p) exceeds
 * a ceratin threshold
 *
 *********************************************************************** */
{
  int    i, j, k;
  double dp;

  /* AYW -- 2012-11-15 16:11 JST */
#if JD_MODE == JD_PRES
  double *x1, *x2, *x3;
  double halo_primitives[NVAR];
  double pres, pres_hot;
#endif
  /* -- AYW */


  /* AYW -- 2012-11-15 16:11 JST */
#if JD_MODE == JD_PRES
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;
#endif
  /* -- AYW */

/* -- backward sweep on x-axis -- */

  if (dir == IDIR){
    for (i = jd_nend; i >= jd_nbeg; i--){
    KDOM_LOOP(k){
    JDOM_LOOP(j){

      /* AYW -- 2012-11-15 16:13 JST */
#if JD_MODE == JD_PRES
      pres = q[k][j][i];
      HotHaloPrimitives(halo_primitives, x1[i], x2[j], x3[k]);
      pres_hot = halo_primitives[PRS];
      if (pres > JD_CONST*pres_hot) return(i);
#else
      dp  = q[k][j][i+1] - q[k][j][i-1];
      dp /= q[k][j][i+1] + q[k][j][i-1];

      if (fabs(dp) > JD_CONST) return(i);
      //if (fabs(dp) > 1.e-5) return(i);
#endif

      /* -- AYW */
    }}}
  }

/* -- backward sweep on y-axis -- */

  if (dir == JDIR){
    for (j = jd_nend; j >= jd_nbeg; j--){
    KDOM_LOOP(k){
    IDOM_LOOP(i){

      /* AYW -- 2012-11-15 16:13 JST */
#if JD_MODE == JD_PRES
      pres = q[k][j][i];
      HotHaloPrimitives(halo_primitives, x1[i], x2[j], x3[k]);
      pres_hot = halo_primitives[PRS];
      if (pres > JD_CONST*pres_hot) return(j);
#else

      dp  = q[k][j+1][i] - q[k][j-1][i];
      dp /= q[k][j+1][i] + q[k][j-1][i];

      //if (fabs(dp) > 1.e-5) return(j);
      if (fabs(dp) > JD_CONST) return(j);
#endif
    }}}
  }

/* -- backward sweep on z-axis -- */

  if (dir == KDIR){
    for (k = jd_nend; k >= jd_nbeg; k--){
    JDOM_LOOP(j){
    IDOM_LOOP(i){

      /* AYW -- 2012-11-15 16:13 JST */
#if JD_MODE == JD_PRES
      pres = q[k][j][i];
      HotHaloPrimitives(halo_primitives, x1[i], x2[j], x3[k]);
      pres_hot = halo_primitives[PRS];
      if (pres > JD_CONST*pres_hot) return(k);
#else

      dp  = q[k+1][j][i] - q[k-1][j][i];
      dp /= q[k+1][j][i] + q[k-1][j][i];

      if (fabs(dp) > JD_CONST) return(k);
      //if (fabs(dp) > 1.e-5) return(k);
#endif
    }}}
  }

  return(-1);
}
