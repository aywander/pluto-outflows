#include "pluto.h"
/* AYW -- 2013-01-08 23:05 JST */
#include "pluto_usr.h"
#include "read_mu_table.h"
#include "interpolation.h"
#include "init_tools.h"
/* -- AYW */

int HuntCooltableIndex(double T, const double *T_tab, const int ntab);

double InterpolatedCoolingRate(double T, const double *T_tab, const double *L_tab, int klo);

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
    static double por_min, L_tmin;
    double Te;

    FILE *fcool;

/* ------------------------------------7XhPmOk%-------
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
        E_cost = vn.t_norm / vn.pres_norm;

        /* Get cooling rate at minimum temperature */
        por_min = g_minCoolingTemp / KELVIN / MU_NORM * vn.pres_norm / vn.dens_norm;
        klo = HuntCooltableIndex(por_min, T_tab, ntab);
        L_tmin = InterpolatedCoolingRate(por_min, T_tab, L_tab, klo);

    }

/* ---------------------------------------------
            Get pressure and temperature 
   --------------------------------------------- */

    /* Get primitive pressure and mean molecular weight. */

    // TODO: Check if this is correct for RHD
    prs = v[RHOE] * (g_gamma - 1.0);
    if (prs < 0.0) {
        prs = g_smallPressure;
        v[RHOE] = prs / (g_gamma - 1.0);
    }

    double prims[NVAR] = {v[RHO], ARG_EXPAND(0, 0, 0), prs, 0, TRC};
#if CLOUDS
    prims[TRC+1] = v[TRC+1] / v[RHO];
#endif
    mu = MeanMolecularWeight(prims);



    /* DM 11 Jul 2015: Now T corresponds to P / rho and not temperature in Kelvin */
    // T   = prs / v[RHO] * KELVIN * mu;
    T = prs / v[RHO] * UNIT_VELOCITY * UNIT_VELOCITY;
    Te = prs / v[RHO] * KELVIN * mu;

    if (T != T) {
        print1(" ! Nan found in radiat \n");
        print1(" ! rho = %12.6e, prs = %12.6e\n", v[RHO], prs);
        QUIT_PLUTO(1);
    }

    // NOTE: This is not consistent with the value of mu from the cooling table, if MU_CALC = MU_CONST
    // AYW -- Add gentle heating below instead of rhs[RHOE] = 0
//  if (T < g_minCoolingTemp) {
//    if (Te < g_minCoolingTemp) {
//        rhs[RHOE] = 0.0;
//        return;
//    }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

    klo = 0;
    khi = ntab - 1;

    if (T > T_tab[khi] || T < T_tab[klo]) {
        print(" ! T out of range   %12.6e %12.6e %12.6e  %12.6e \n", T, prs, v[RHO], prs / v[RHO] * KELVIN * mu);
        QUIT_PLUTO(1);
    }

    klo = HuntCooltableIndex(T, T_tab, ntab);

/* -----------------------------------------------
    Compute r.h.s
   ----------------------------------------------- */

    /* Cooling rate over mu squared. Mu is the mean mass per particle in units of atomic units (amu).
     * Mu itself does not contain amu, we divide by amu^2 below. Note, Ecost is normalization for
     * Lambda * n^2 [erg cm-3 s-1] into code units when multiplied. */

    scrh = InterpolatedCoolingRate(T, T_tab, L_tab, klo);

    rhs[RHOE] = -(scrh - L_tmin);                     // Cooling and Heating
    rhs[RHOE] *=  v[RHO] * v[RHO];

    scrh       = UNIT_DENSITY / (CONST_amu * mu);
    rhs[RHOE] *= E_cost * scrh * scrh;
}

double InterpolatedCoolingRate(double T, const double *T_tab, const double *L_tab, int klo) {
    double scrh;
    double dT;
    int khi = klo + 1;
    dT = T_tab[khi] - T_tab[klo];
    scrh = L_tab[klo] * (T_tab[khi] - T) / dT + L_tab[khi] * (T - T_tab[klo]) / dT;
    return scrh;
}

int HuntCooltableIndex(double T, const double *T_tab, const int ntab) {

    int klo, khi, kmid;
    double Tmid;

    klo = 0;
    khi = ntab - 1;

    while (klo != (khi - 1)) {
        kmid = (klo + khi) / 2;
        Tmid = T_tab[kmid];
        if (T <= Tmid) {
            khi = kmid;
        } else if (T > Tmid) {
            klo = kmid;
        }
    }
    return klo;
}
