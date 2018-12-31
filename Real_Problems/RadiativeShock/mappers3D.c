/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief 3D wrapper for conservative/primitive conversion.

  Provide 3D wrappers to the standard 1D conversion functions
  ConsToPrim() and PrimToCons().

  \authors A. Mignone (mignone@ph.unito.it)
  \date    June 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "macros_usr.h"

/* ********************************************************************* */
void ConsToPrim3D (Data_Arr U, Data_Arr V, unsigned char ***flag, RBox *box)
/*!
 *  Convert a 3D array of conservative variables \c U to
 *  an array of primitive variables \c V.
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 * \param [in]   U      pointer to 3D array of conserved variables,
 *                      with array indexing <tt>[k][j][i][nv]</tt>
 * \param [out]  V      pointer to 3D array of primitive variables,
 *                      with array indexing <tt>[nv][k][j][i]</tt>
 * \param [in,out] flag   pointer to 3D array of flags.
 * \param [in]     box    pointer to RBox structure containing the domain
 *                        portion over which conversion must be performed.
 *
 *********************************************************************** */
{
  int   i, j, k, nv, err;
  int   ibeg, iend, jbeg, jend, kbeg, kend;
    int current_dir;
    static double **v, **u;
    double prsfix, g, scrh, rhofix, engfix, m2;

    if (v == NULL) {
      v = ARRAY_2D(NMAX_POINT, NVAR, double);
      u = ARRAY_2D(NMAX_POINT, NVAR, double);
    }

/* ----------------------------------------------
    Save current sweep direction and by default,
    perform the conversion along X1 stripes
   ---------------------------------------------- */

  current_dir = g_dir;
    g_dir = IDIR;

/* -----------------------------------------------
    Set (beg,end) indices in ascending order for
    proper call to ConsToPrim()
   ----------------------------------------------- */

  ibeg = (box->ibeg <= box->iend) ? (iend=box->iend, box->ibeg):(iend=box->ibeg, box->iend);
  jbeg = (box->jbeg <= box->jend) ? (jend=box->jend, box->jbeg):(jend=box->jbeg, box->jend);
  kbeg = (box->kbeg <= box->kend) ? (kend=box->kend, box->kbeg):(kend=box->kbeg, box->kend);

  for (k = kbeg; k <= kend; k++){ g_k = k;
  for (j = jbeg; j <= jend; j++) {
      g_j = j;
#if (defined CHOMBO) && (COOLING == MINEq || COOLING == H2_COOL)
      if (g_intStage == 1) {
        for (i = ibeg; i <= iend; i++)  NormalizeIons(U[k][j][i]);
      }
#endif
      err = ConsToPrim(U[k][j], v, ibeg, iend, flag[k][j]);
      for (i = ibeg; i <= iend; i++) NVAR_LOOP(nv) V[nv][k][j][i] = v[i][nv];
  }}


  /* AYW --
   * Based on modifications by DM in 2015.
   * This routine has changed a lot from PLUTO 4.1 to 4.2.
   * Simplified the smoothing here - macros RHO_FAIL, PRS_FAIL etc, don't exist anymore.
   * Instead, there is a FLAG_CONS2PRIM_FAIL, for all cases.
   * So we smooth density and pressure for this cell. */

  for (k = kbeg; k <= kend; k++) { g_k = k;
  for (j = jbeg; j <= jend; j++) { g_j = j;
  for (i = ibeg; i <= iend; i++) {

            // TODO: AYW -- This is probalby not the best way of doing this.

            if ((flag[k][j][i] & FLAG_CONS2PRIM_FAIL) == FLAG_CONS2PRIM_FAIL) {

                // Density fix

                rhofix = BURY(V[RHO], k, j, i);

                if (rhofix <= 0.0) {
                    print("RHO still negative [%d,%d,%d] \n", i, j, k);
                    rhofix = g_smallDensity;
                }

                if (rhofix != rhofix) {
                    print("RHO is NAN [%d,%d,%d] \n", i, j, k);
                    rhofix = g_smallDensity;
                }

                v[i][RHO] = rhofix;

                // Pressure fix

                prsfix = BURY(V[PRS], k, j, i);

                if (prsfix <= 0.0) {
                    print("PRS still negative [%d,%d,%d] \n", i, j, k);
                    prsfix = g_smallPressure;
                }

                if (prsfix != prsfix) {
                    print("PRS is NAN [%d,%d,%d] \n", i, j, k);
                    prsfix = g_smallPressure;
                }

                v[i][PRS] = prsfix;


            } // CONS2PRIM_FAIL
  }}} // k, j, i-loops
  g_dir = current_dir;

}
/* ********************************************************************* */
void PrimToCons3D (Data_Arr V, Data_Arr U, RBox *box)
/*!
 *  Convert a 3D array of primitive variables \c V  to
 *  an array of conservative variables \c U.
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 * \param [in]    V     pointer to 3D array of primitive variables,
 *                      with array indexing <tt>[nv][k][j][i]</tt>
 * \param [out]   U     pointer to 3D array of conserved variables,
 *                      with array indexing <tt>[k][j][i][nv]</tt>
 * \param [in]    box   pointer to RBox structure containing the domain
 *                      portion over which conversion must be performed.
 *
 *********************************************************************** */
{
    int i, j, k, nv;
    int   ibeg, iend, jbeg, jend, kbeg, kend;
    int current_dir;
    static double **v, **u;

    if (v == NULL) {
      v = ARRAY_2D(NMAX_POINT, NVAR, double);
      u = ARRAY_2D(NMAX_POINT, NVAR, double);
    }

    current_dir = g_dir; /* save current direction */
    g_dir = IDIR;

/* -----------------------------------------------
    Set (beg,end) indices in ascending order for
    proper call to ConsToPrim()
   ----------------------------------------------- */

  ibeg = (box->ibeg <= box->iend) ? (iend=box->iend, box->ibeg):(iend=box->ibeg, box->iend);
  jbeg = (box->jbeg <= box->jend) ? (jend=box->jend, box->jbeg):(jend=box->jbeg, box->jend);
  kbeg = (box->kbeg <= box->kend) ? (kend=box->kend, box->kbeg):(kend=box->kbeg, box->kend);

  for (k = kbeg; k <= kend; k++){ g_k = k;
  for (j = jbeg; j <= jend; j++){ g_j = j;
    for (i = ibeg; i <= iend; i++) NVAR_LOOP(nv) v[i][nv] = V[nv][k][j][i];
    PrimToCons (v, U[k][j], ibeg, iend);
  }}
    g_dir = current_dir; /* restore current direction */
}
