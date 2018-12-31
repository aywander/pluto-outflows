/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sepy 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include <assert.h>
#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"
#include "shock.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rd dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{

  static int once0 = 0;

  /* Some initializations */
  if (!once0) {
    SetBaseNormalization();
    SetIniNormalization();
  }


  // TODO: Allow variable mu upstream and downstream
  double mu = MU_NORM;

  /* Input parameters */
  double mach = g_inputParam[PAR_MACH] * ini_code[PAR_MACH];
  double spd = g_inputParam[PAR_SSPD] * ini_code[PAR_SSPD];
  double rho = g_inputParam[PAR_DENS] * ini_code[PAR_DENS];
  double te = g_inputParam[PAR_TEMP] * ini_code[PAR_TEMP];
  double bprp = g_inputParam[PAR_BPRP] * ini_code[PAR_BPRP];
  double bpar = g_inputParam[PAR_BPAR] * ini_code[PAR_BPAR];
  double shock_loc = g_inputParam[PAR_SLOC] * ini_code[PAR_SLOC];
  double minCoolingTemp = g_inputParam[PAR_TMIN];

  /* Some derived parameters. Fill shock structure struct.
   * mach is always the distant upstream mach number. */
  double csound;

  // TODO: Do MHD jumps (May not be possible analytically for this mode)


#if SHOCK_COND_MODE == SC_UPSTREAM

  csound = sqrt(g_gamma * te / mu);
  if (mach <= 0) mach = spd / csound;


  sh.s0.rho = rho;
  sh.s0.u = mach * csound;
  sh.s0.te = te;
  sh.s0.mu = mu;
  sh.s0.bprp = bprp;
  sh.s0.bpar = bpar;

  sh.sp.rho = sh.s0.rho;
  sh.sp.u = sh.s0.u;
  sh.sp.te = sh.s0.te;
  sh.sp.mu = sh.s0.mu;
  sh.sp.bprp = sh.s0.bprp;
  sh.sp.bpar = sh.s0.bpar;

  // TODO: Do MHD jumps
#if PHYSICS == MHD || PHYSICS == RMHD

   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = 0.0;

   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = 0.0;
#endif

  sh.s1.rho = sh.s0.rho * rho_jump(mach);
  sh.s1.u = sh.s0.u / rho_jump(mach);
  sh.s1.te = sh.s0.te * te_jump(mach);
  sh.s1.mu = sh.s0.mu;

  sh.sd.te = minCoolingTemp < 0 ? sh.s0.te : minCoolingTemp / KELVIN;
  sh.sd.rho = sh.s0.rho * rho_jump_cdl(mach, sh.s0.te, sh.sd.te);
  sh.sd.u = sh.s0.u / rho_jump_cdl(mach, sh.s0.te, sh.sd.te);
  sh.sd.mu = sh.s0.mu;

#elif SHOCK_COND_MODE == SC_SHOCK_TE

  if (mach <= 0) {
    mach = mach_from_te1(spd, te, mu);
  }

  if (mach <= 1) {
    print("Solution unphysical: Mach number <= 1.");
    QUIT_PLUTO(1);
  }

  sh.s0.rho = rho;
  sh.s0.te = te / te_jump(mach);
  csound = sqrt(g_gamma * sh.s0.te / mu);
  sh.s0.u = mach * csound;
  sh.s0.mu = mu;
  sh.s0.bprp = bprp;
  sh.s0.bpar = bpar;

  sh.sp.rho = sh.s0.rho;
  sh.sp.u = sh.s0.u;
  sh.sp.te = sh.s0.te;
  sh.sp.mu = sh.s0.mu;
  sh.sp.bprp = sh.s0.bprp;
  sh.sp.bpar = sh.s0.bpar;

  // TODO: Do MHD jumps
#if PHYSICS == MHD || PHYSICS == RMHD

   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = 0.0;

   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = 0.0;
#endif

  sh.s1.rho = sh.s0.rho * rho_jump(mach);
  sh.s1.u = sh.s0.u / rho_jump(mach);
  sh.s1.te = te;
  sh.s1.mu = sh.s0.mu;

  sh.sd.te = minCoolingTemp < 0 ? sh.s0.te : minCoolingTemp / KELVIN;
  sh.sd.rho = sh.s0.rho * rho_jump_cdl(mach, sh.s0.te, sh.sd.te);
  sh.sd.u = sh.s0.u / rho_jump_cdl(mach, sh.s0.te, sh.sd.te);
  sh.sd.mu = sh.s0.mu;

#endif

#if SHOCK_INIT_MODE == SI_IMPULSIVE

  // TODO: Generalize to 3D
  /* Location of shock center */
  double xc = g_domBeg[IDIR] + shock_loc * (g_domEnd[IDIR] - g_domBeg[IDIR]);

  /* Impulsive shock generation setup. Upstream is right in order to use -x1jet. */
  if (x1 < xc){
    v[RHO] = sh.s0.rho;
    v[VX1] = sh.s0.u;
    v[PRS] = sh.s0.rho * sh.s0.te / sh.s0.mu;
  }

  else {
    v[RHO] = sh.sd.rho;
    v[VX1] = sh.sd.u;
    v[PRS] = sh.sd.rho * sh.sd.te / sh.sd.mu;
  }


#elif SHOCK_INIT_MODE == SI_WALL

  v[RHO] = sh.s0.rho;
  v[VX1] = sh.s0.u - sh.s1.u;
  v[PRS] = sh.s0.rho * sh.s0.te / sh.s0.mu;

#endif


  /* Cooling integration constraints */
  /* Remember, there is a hack in cooling_source.c that prevents upstream gas from cooling */
  g_minCoolingTemp = minCoolingTemp < 0 ? sh.s0.te * KELVIN : minCoolingTemp;
  g_maxCoolingRate = 0.1;

  /* Print the initial conditions */
  if (!once0) {
    PrintInitData01();
    once0 = 1;
  }

}


/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){};
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  return 0.0;
}
#endif
