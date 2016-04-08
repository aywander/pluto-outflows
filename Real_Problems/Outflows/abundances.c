#include "abundances.h"
#include "init_tools.h"

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
double MeanMolecularWeight (real *V)
/* This function is reproduced here because without the cooling function
 * the mean molecular weight function is not compiled in. 
 *   Changing the abundances every time a new problem is set up is 
 * tedious, though.
 *
 ********************************************************************* */
{
#if MU_CALC == MU_TABLE

    int il;
    double por;

    if (mu_por == NULL) {
        ReadMuTable();
    }

    /* Value of p/rho in code units */
    por = V[PRS] / V[RHO];

    /* Interpolate */
    return InterpolationWrapper(mu_por, mu_mu, mu_ndata, por);


#elif MU_CALC == MU_FRACTIONS

    return  ( (A_H + frac_He*A_He + frac_Z*A_Z) /
                  (2.0 + frac_He + 2.0*frac_Z - 0.0));

#elif MU_CALC == MU_ANALYTIC

    /* Value of log10(p/rho) in cgs units */
    double por;
    por = log10(V[PRS] / V[RHO] * vn.pres_norm / vn.dens_norm);

    static double a = 11.48648535;
    static double w = 0.62276904;
    static double m = 1.25;
    double tanh_factor;

    tanh_factor = tanh((por - a) / w);
    return 0.5 * (MU_NORM + m) + 0.5 * tanh_factor * (MU_NORM - m);

#elif MU_CALC == MU_CONST

    return (MU_NORM);

#else

    return (MU_NORM);

#endif

}

#endif
