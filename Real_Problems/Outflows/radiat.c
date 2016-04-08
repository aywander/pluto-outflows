#include "pluto.h"
/* AYW -- 2013-01-08 23:05 JST */
#include "pluto_usr.h"
#include "read_mu_table.h"
#include "interpolation.h"
#include "init_tools.h"

/* -- AYW */

#define frac_Z   1.e-3   /* = N(Z) / N(H), fractional number density of 
                              metals (Z) with respect to hydrogen (H) */ 
#define frac_He  0.082   /* = N(Z) / N(H), fractional number density of 
                              helium (He) with respect to hydrogen (H) */ 
/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */
{
    int klo, khi, kmid;
    static int ntab;
    double mu, T, Tmid, scrh, dT, prs;
    static double *L_tab, *T_tab, E_cost;

    FILE *fcool;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

    if (T_tab == NULL) {
        print1(" > Reading table from disk...\n");
        /*  Read in input table with 1st column as P/rho in cgs and
            second column being Lambda/mu^2 */
        fcool = fopen("cooltable.dat", "r");
        if (fcool == NULL) {
            print1("! Radiat: cooltable.dat could not be found.\n");
            QUIT_PLUTO(1);
        }
        L_tab = ARRAY_1D(20000, double);
        T_tab = ARRAY_1D(20000, double);

        ntab = 0;
        while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab,
                      L_tab + ntab) != EOF) {
            ntab++;
        }
        /* Normalization for Lambda * n^2 [erg cm-3 s-1] into code units when multiplied */
        E_cost = UNIT_LENGTH / UNIT_DENSITY / pow(UNIT_VELOCITY, 3.0);
    }

/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */

    // TODO: Check if this is correct for RHD
    prs = v[RHOE] * (g_gamma - 1.0);
    if (prs < 0.0) {
        prs = g_smallPressure;
        v[RHOE] = prs / (g_gamma - 1.0);
    }

    mu = MeanMolecularWeight(v);

    /* DM 11 Jul 2015: Now T corresponds to P / rho and not temperature in Kelvin */
    // T   = prs / v[RHO] * KELVIN * mu;
    T = prs / v[RHO] * UNIT_VELOCITY * UNIT_VELOCITY;

    if (T != T) {
        print1(" ! Nan found in radiat \n");
        print1(" ! rho = %12.6e, prs = %12.6e\n", v[RHO], prs);
        QUIT_PLUTO(1);
    }

    // NOTE: This is not consistent with the value of mu from the cooling table, if MU_CALC = MU_CONST
//  if (T < g_minCoolingTemp) { 
    if (prs / v[RHO] * KELVIN * mu < g_minCoolingTemp) {
        rhs[RHOE] = 0.0;
        return;
    }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

    klo = 0;
    khi = ntab - 1;

    if (T > T_tab[khi] || T < T_tab[klo]) {
        print(" ! T out of range   %12.6e %12.6e %12.6e  %12.6e \n", T, prs, v[RHO], prs / v[RHO] * KELVIN * mu);
        QUIT_PLUTO(1);
    }

    while (klo != (khi - 1)) {
        kmid = (klo + khi) / 2;
        Tmid = T_tab[kmid];
        if (T <= Tmid) {
            khi = kmid;
        } else if (T > Tmid) {
            klo = kmid;
        }
    }

    /* Cooling rate over mu squared. Mu is the mean mass per particle in units of atomic units (amu).
     * Mu itself does not contain amu, we divide by amu^2 below. Note, Ecost is normalization for
     * Lambda * n^2 [erg cm-3 s-1] into code units when multiplied. */
    dT = T_tab[khi] - T_tab[klo];
    scrh = L_tab[klo] * (T_tab[khi] - T) / dT + L_tab[khi] * (T - T_tab[klo]) / dT;
    rhs[RHOE] = -scrh * v[RHO] * v[RHO];
    // rhs[RHOE] *= E_cost * UNIT_DENSITY * UNIT_DENSITY / (CONST_amu * CONST_amu);
    rhs[RHOE] *= E_cost * UNIT_DENSITY * UNIT_DENSITY / (CONST_amu * CONST_amu * mu * mu);
}
#undef T_MIN

/* ******************************************************************* */
double MeanMolecularWeight (double *V)
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

