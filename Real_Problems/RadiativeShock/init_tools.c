#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"
#include "abundances.h"

/* ********************************************************************* */
void RankineHugoniotJump(double *v1, double *v2, double *temp2, int isoth){
/*! 
 * Rankine Hugoniot Jump conditions for stationary shock.
 *
 * \param [out] v1      a pointer to a vector of primitive downstream variables
 * \param [out] temp2   temperature behind shock
 * \param [in]  v0      a pointer to a vector of primitive upstream variables
 * \param [in]  isoth   a switch as to whether shock jumps are 
 *                      adiabatic (isoth == 0) or isothermal (isoth != 0).
 * 
 *********************************************************************** */

  /* Shorthand for v0 components */
  double rho = v1[RHO];
  EXPAND(double vx1 = v1[VX1];,
         double vx2 = v1[VX2];,
         double vx3 = v1[VX3];);
  double prs = v1[PRS];

  /* Calculate temperature too */
  real mu = MeanMolecularWeight(v1);
  double temp1 = mu*prs/rho;

  /* Calculate mach number */
  double vel = sqrt(EXPAND(vx1*vx1, + vx2*vx2, + vx3*vx3));
  double csound = sqrt(g_gamma*v1[PRS]/v1[RHO]);
  double mach = vel/csound;

  /* Differentiate between isothermal and adiabatic jump conditions */
  double gamma;
  if (isoth == 0) gamma = g_gamma;
  else            gamma = 1.;

  /* Caluclate jumps 
   * The factor g_gamma/gamma generalizes to both isothermal and adiabatic jumps */
  double r = (g_gamma/gamma)*(gamma+1.)*mach*mach/((gamma-1.)*mach*mach + 2.);
  v2[RHO] = rho*r;
  EXPAND(v2[VX1] = vx1/r;, v2[VX2] = vx2/r;, v2[VX3] = vx3/r;);
  v2[PRS] = prs*(g_gamma/gamma)*(2.*gamma*mach*mach - (gamma-1.))/(gamma+1.);
  *temp2  = temp1*(((gamma-1.)*mach*mach + 2.)*
          (2.*gamma*mach*mach - (gamma-1.)))/
          ((gamma+1.)*(gamma+1.)*mach*mach);

}

