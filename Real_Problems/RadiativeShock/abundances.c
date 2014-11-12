#include "pluto.h"
#include "pluto_usr.h"
#include "abundances.h"

#if COOLING == NO

#define frac_Z   1.e-3   /*   = N(Z) / N(H), fractional number density of metals (Z)
                                with respect to hydrogen (H) */ 
#define frac_He  0.0574  /*    This would bring mmw in line with what we used to use 
                                but is inconsistent with neqstd2009.T4 abundance pattern see file */
/*#define frac_He  0.082      = N(Z) / N(H), fractional number density of helium (He) 
                                with respect to hydrogen (H) */ 
#define A_Z      30.0    /*   mean atomic weight of heavy elements  */
#define A_He     4.004   /*   atomic weight of Helium  */
#define A_H      1.008   /*   atomic weight of Hydrogen  */

/* ******************************************************************* */
real MeanMolecularWeight (real *V)
/* This function is reproduced here because without the cooling function
 * the mean molecular weight function is not compiled in. 
 *   Changing the abundances every time a new problem is set up is 
 * tedious, though.
 *
 ********************************************************************* */
{
  //return (0.5)
  return (0.6156);
  //return  ( (A_H + frac_He*A_He + frac_Z*A_Z) /
  //           (2.0 + frac_He + 2.0*frac_Z - 0.0)); 
  
}

#endif
