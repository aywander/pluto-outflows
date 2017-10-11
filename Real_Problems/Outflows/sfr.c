//
// Created by Alexander Y. Wagner on 2017/06/29.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "clouds.h"
#include "sfr.h"
#include "accretion.h"
#include "idealEOS.h"

/*********************************************************************** */
void FederrathAccretion(Data *d, const Grid *grid) {
/*!
 * Driver routine for removing mass in gas around central BH
 * according to Federrath accretion.
 *
 *********************************************************************** */
    // TODO: Reprogram this for star-formation

    /* Accretion - General variables */
    int i, j, k;
    double *x1, *x2, *x3;
    double *dV1, *dV2, *dV3;
    double accr, accr_rate = 0;
    double result[NVAR];


    /* These are the geometrical central points */
    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    /* These are cell volumes */
    dV1 = grid[IDIR].dV;
    dV2 = grid[JDIR].dV;
    dV3 = grid[KDIR].dV;


    DOM_LOOP(k, j, i) {
                /* Remove mass according to Federrath's sink particle method */
                accr = FederrathSinkInternalBoundary(d->Vc, i, j, k, x1, x2, x3, dV1, dV2, dV3, result);
                accr /= g_dt;
                accr_rate += accr;
            }


    /* MPI reductions and analysis */

#ifdef PARALLEL
        MPI_Allreduce(&accr_rate, &ac.accr_rate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        ac.accr_rate = accr_rate;
#endif

    /* Increase BH mass by measured accretion rate * dt */
    ac.mbh += ac.accr_rate * g_dt * (1. - ac.eff);

    /* Update Eddington luminosity - calculate after increasing BH mass */
    ac.edd = EddingtonLuminosity(ac.mbh);

}

/* ************************************************ */
double JeansResolvedDensity(const double *prim){
/*!
 * Return the difference between the cell density and maximum density
 * obtained from the criterion that the Jeans length should be resolved
 * within the sink region is satisfied.
 *
 ************************************************** */

    // TODO: Reprogram this for star-formation
    double c2 = SoundSpeed2IdealEOS(prim[RHO], prim[PRS]);

    return CONST_PI * c2 * vn.newton_norm / (ac.snk * ac.snk * CONST_G);

}

/* ************************************************ */
double FederrathSinkInternalBoundary(const double ****Vc, int i, int j, int k, const double *x1, const double *x2,
                                     const double *x3, const double *dV1, const double *dV2, const double *dV3,
                                     double *result) {
/*!
 * Remove delta_mass in cells so as to satisfy the Jeans criterion
 * that the Jeans length should be
 * resolved by at least 5 cells.
 *
 ************************************************** */

    // TODO: Reprogram this for star-formation

    /* Accretion - General variables */
    int iv;
    double xi, xj, xk;
    double vx1, vx2, vx3;
    double vsph;
    double cell_vol, delta_rho, delta_mass;

    xi = x1[i];
    xj = x2[j];
    xk = x3[k];

    /* Default is that nothing is changed */
    for (iv = 0; iv < NVAR; iv++) result[iv] = Vc[iv][k][j][i];

    /* Check whether velocity is pointing inward toward sink
     * currently fixed at (0, 0, 0) */
    EXPAND(vx1 = result[VX1];,
           vx2 = result[VX2];,
           vx3 = result[VX3];);
    vsph = VSPH1(xi, xj, xk, vx1, vx2, vx3);
    if (vsph > 0) return 0;

    /* Jeans density criterion */
    delta_rho = result[RHO] - JeansResolvedDensity(result);
    if (delta_rho < 0) return 0;

    /* Accreted Mass */
    cell_vol = dV1[i] * dV2[j] * dV3[k];
    delta_mass = delta_rho * cell_vol;

    /* Check whether gas mass is gravitationally bound */
    if (!(GravitationallyBound(result, delta_mass, cell_vol, xi, xj, xk))) return 0;

    /* Remove mass */
    result[RHO] -= delta_mass / cell_vol;

    return delta_mass;

}

/* ************************************************ */
double VirialParameter(const double * prim, const double mass,
                       const double x1, const double x2, const double x3){
/*!
 * Compute the virial parameter 2 * Ekin / |Epot| (Federrath & Klessen 2012, eq. 16)
 * of gas in a cell relative to a stationary point (sink particle motion = 0).
 *
 ************************************************** */
    // TODO: Reprogram this for star-formation

    EXPAND(double vx1 = prim[VX1];,
           double vx2 = prim[VX2];,
           double vx3 = prim[VX3];);

    double vmag = VMAG(x1, x2, x3, vx1, vx2, vx3);
    double sph1 = SPH1(x1, x2, x3);
    double rho = prim[RHO];

    double ekin = 0.5 * rho * vmag * vmag;
    double epot = CONST_G * mass * rho / (sph1 * vn.newton_norm);

    return 2 * ekin / epot;

}

/* ************************************************ */
double GravitationallyBound(const double *prim, const double mass, const double vol,
                            const double x1, const double x2, const double x3){
/*!
 * Check whether mass element mass is gravitationally bound.
 * Velocity is relative to a stationary sink at (0, 0, 0).
 * E_pot + Ekin + Eth + Emag > 0;
 *
 ************************************************** */

    // TODO: Reprogram this for star-formation

    EXPAND(double vx1 = prim[VX1];,
           double vx2 = prim[VX2];,
           double vx3 = prim[VX3];);

    double vmag = VMAG(x1, x2, x3, vx1, vx2, vx3);
    double rho = prim[RHO];
    double prs = prim[PRS];

    double ekin = 0.5 * mass * vmag * vmag;
    double epot = 0;  // TODO: Update this when self-gravity becomes available
    double eth = 1. / (g_gamma - 1.) * prs * vol;

#if PHYSICS == MHD || PHYSICS == RMHD
    EXPAND(double bx1 = prim[BX1];,
           double bx2 = prim[BX2];,
           double bx3 = prim[BX3];);
    double bmag = VMAG(x1, x2, x3, bx1, bx2, bx3);
    double emag = 1. / (8. * CONST_PI) * bmag * bmag * vol;
#else
    double emag = 0;
#endif

    return (epot + ekin + eth + emag < 0);

}


double StarFormationRateDensity(Data *d, Grid *grid, int i, int j, int k) {

    double sfr, sfe;
    double rho, ff_time;

    /* Total mass in domain */

    sfe = StarFormationEfficiency(d, grid, i, j, k);

    /* Primitives */
    rho = d->Vc[RHO][k][j][i];

    sfr = 0;
    if ( StarFormationCriteria(d, grid, i, j, k) ) {
        ff_time = 0.5427 / sqrt(CONST_G / vn.newton_norm * rho);
        sfr = sfe * rho / ff_time;
    }

    return sfr;

}


double StarFormationEfficiency(Data *d, Grid *grid, int i, int j, int k) {

    return 0.02;

}

double StarFormationCriteria(Data *d, Grid *grid, int i, int j, int k) {

    double rho, tr2, te;
    double vx1, vx2, vx3;

    /* Region in phase to consider as warm and where stars can form */
    double rho_c, te_c, tr2_c;
    WarmPhaseConditions(&rho_c, &te_c, &tr2_c);

    /* Primitives */
    rho = d->Vc[RHO][k][j][i];
    te = d->Vc[PRS][k][j][i] / rho * MU_NORM * KELVIN;
    tr2 = d->Vc[TRC+1][k][j][i];

    // TODO: add more criteria (div v < 0, virial param, grav. bound)
    if (te < te_c && rho < rho_c
//      && tr2 > tr2_c
            )
    {
        return 1;
    }
    else {
        return 0;
    }



}