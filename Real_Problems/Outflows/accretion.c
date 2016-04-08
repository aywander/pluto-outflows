//
// Created by Alexander Y. Wagner on 4/6/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "accretion.h"
#include "grid_geometry.h"
#include "idealEOS.h"
#include "outflow.h"


/* Global struct for accretion */
Accretion ac;


/* ************************************************ */
void SetAccretionPhysics() {
/*
 * Set values for the accretion structure
 *
 ************************************************** */


    /* Quantities from input parameters */
    ac.rad = g_inputParam[PAR_ARAD] * ini_code[PAR_ARAD];
    ac.mbh = g_inputParam[PAR_AMBH] * ini_code[PAR_AMBH];
    ac.eff = g_inputParam[PAR_AEFF] * ini_code[PAR_AEFF];
    ac.snk = g_inputParam[PAR_ASNK] * ini_code[PAR_ASNK];

    /* Bondi Accretion parameters */
    ac.rho_acc = g_smallDensity;
    ac.prs_acc = g_smallPressure;
    ac.snd_acc = sqrt(g_gamma * g_smallPressure / g_smallDensity);
    ac.accr_rate_bondi = BondiAccretionRate(ac.mbh, ac.rho_acc, ac.snd_acc);

    // TODO: Some of these quantities need to be uptdated after restart.

    /* Eddington rate */
    ac.edd = EddingtonLuminosity(ac.mbh);

    /* Global accretion rate */
    ac.accr_rate = 0;

    /* Area of accretion surface. */
    ac.area = 4 * CONST_PI * ac.rad * ac.rad;

    /* Number of cells in the spherical surface of radius ac.rad is
     * calculated in Analysis analysis */
}

/* ************************************************ */
int InSinkRegion(const double x1, const double x2, const double x3) {
/*
 * Returns 1 if sphere intersects a cartesian cell local domain.
 * Currently only for Cartesian and cylindrical case,
 * but other geometries are not difficult.
 *
 * This is almost identical to InFlankRegion, except for the radii.
 *
 * Sphere is assumed to be at (0,0,0) and has a radius
 * of rad.
 ************************************************** */

    /* Nozzle is be deactivated if OSPH is zero */
    if (ac.snk == 0.) return 0;

    /* Do radius check. Get spherical r coordinate */
    double sr = SPH1(x1, x2, x3);
    if (sr > ac.snk) return 0;

    /* Grid point in cartesian coordinates */
    double cx1, cx2, cx3;
    cx1 = CART1(x1, x2, x3);
    cx2 = CART2(x1, x2, x3);
    cx3 = CART3(x1, x2, x3);

    /* Rotate cartesian coordinates */
    double cx1p, cx2p, cx3p;
    RotateGrid2Nozzle(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);

    /* For internal boundary, we include a mirrored component with fabs.
     * This must occur after rotation but before translation. */
#if INTERNAL_BOUNDARY == YES
    D_SELECT(cx1p, cx2p, cx3p) = fabs(D_SELECT(cx1p, cx2p, cx3p));
#endif

    if (nz.isfan) {

        /* Shift so that cone apex is at (0,0,0) */
        D_SELECT(cx1p, cx2p, cx3p) -= nz.cone_apex;

        /* Do angle checks. Get spherical theta coordinate. */
        double st = CART2SPH2(cx1p, cx2p, cx3p);
        return (st > nz.ang);

    }
    else {

        /* Do radius check. Get cylindrical coordinate. */
        double cr = CART2CYL1(cx1p, cx2p, cx3p);
        return (cr > nz.rad);

    }

}



/* ************************************************ */
void SphericalAccretion(const Data *d, Grid *grid) {
/*!
 * Calculate the spherical accretion rate through the surface of a sphere.
 *
 ************************************************** */

    /* Accretion */

    int i, j, k;
    double rho, vs1;
    double vx1, vx2, vx3;
    double *x1, *x2, *x3;
    double accr, accr_rate = 0;

#if SINK_METHOD == SINK_BONDI
    int gcount = 0, lcount = 0;
    double rho_acc = 0, snd_acc = 0, prs, prs_acc = 0;
    double tmp_far, snd_far;
#endif


    /* Measure accretion rate through a spherical surface defined by ac->rad (ARAD) */

    if (SphereSurfaceIntersectsDomain(grid, ac.rad)) {

        /* Get number of cells lying on spherical surface ac.rad. */

        static int once = 0;
        static int ncells;
        static double area_per_cell;

        if (!once) {
            ncells = SphereSurfaceIntersectsNCells(grid[IDIR].dx[0], grid[JDIR].dx[0], grid[KDIR].dx[0], ac.rad);
            area_per_cell = ac.area / ncells;
            once = 1;
        }


        /* These are the geometrical central points */
        x1 = grid[IDIR].x;
        x2 = grid[JDIR].x;
        x3 = grid[KDIR].x;

        /* Intersection "booleans" */
        int sicr = 0;
        int sicc = 0;

        DOM_LOOP(k, j, i) {

#if SIC_METHOD == SIC_RADIUS || SIC_METHOD == SIC_HYBRID
                    sicr = SphereSurfaceIntersectsCellByRadius(x1[i], x2[j], x3[k],
                                                               grid[IDIR].dV[i], grid[JDIR].dV[j], grid[KDIR].dV[k],
                                                               ac.rad);
#endif

#if SIC_METHOD == SIC_CORNERS || SIC_METHOD == SIC_HYBRID
                    sicc = SphereSurfaceIntersectsCellByCorners(x1[i], x2[j], x3[k],
                                                                grid[IDIR].dx[i], grid[JDIR].dx[j], grid[KDIR].dx[k],
                                                                ac.rad);
#endif
                    /* Only for cells "on" the surface, where we measure the accretion rate. */
                    if (sicr || sicc) {

                        /* Calculate and sum accretion rate */
                        rho = d->Vc[RHO][k][j][i];
                        vx1 = vx2 = vx3 = 0;
                        EXPAND(vx1 = d->Vc[VX1][k][j][i];,
                               vx2 = d->Vc[VX2][k][j][i];,
                               vx3 = d->Vc[VX3][k][j][i];);
                        vs1 = VSPH1(x1[i], x2[j], x3[k], vx1, vx2, vx3);
                        vs1 = fabs(MIN(vs1, 0));

                        accr = rho * vs1 * area_per_cell;
#if SIC_METHOD == SIC_HYBRID
                        if (sicr && sicc) accr *= 2;
#endif
                        accr_rate += accr;


#if SINK_METHOD == SINK_BONDI
                        /* Bondi accretion parameters */
                        rho_acc += rho;
                        prs = d->Vc[PRS][k][j][i];
                        prs_acc += prs;
                        snd_acc += sqrt(g_gamma * prs / rho);
                        lcount++;
#endif

                    } // if sicr || sicc

                }  // DOM_LOOP

        /* Averaging in HYBRID mode */
#if SIC_METHOD == SIC_HYBRID
        accr_rate /= 2;
#endif

    }



    /* Reductions, including MPI reductions,
     * and calculation of other quantities. */

#ifdef PARALLEL
    MPI_Allreduce(&accr_rate, &ac.accr_rate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#if SINK_METHOD == SINK_BONDI

    /* MPI reduce Bondi accretion parameters */
    MPI_Allreduce(&lcount, &gcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&rho_acc, &ac.rho_acc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&prs_acc, &ac.prs_acc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&snd_acc, &ac.snd_acc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    ac.rho_acc /= gcount;
    ac.prs_acc /= gcount;
    ac.snd_acc /= gcount;
#endif

    MPI_Barrier(MPI_COMM_WORLD);

#else // If not parallel

    ac.accr_rate = accr_rate;

#if SINK_METHOD == SINK_BONDI

    // TODO: Do we need all these in the structure?
    ac.rho_acc = rho_acc / lcount;
    ac.prs_acc = prs_acc / lcount;
    ac.snd_acc = snd_acc / lcount;

#endif

#endif // ifdef PARALLEL

    /* Bondi accretion rate - calculate before increasing BH mass */
#if SINK_METHOD == SINK_BONDI

    /* Bondi rate rewritten relating far-away density to density and sound speed in accretion region. */
    tmp_far = g_inputParam[PAR_HTMP] * ini_cgs[PAR_HTMP];
    snd_far = sqrt(g_gamma * CONST_kB * tmp_far / (MU_NORM * CONST_amu)) / vn.v_norm;
    ac.accr_rate_bondi = BondiAccretionRateLocal(ac.mbh, ac.rho_acc, ac.snd_acc, snd_far);

#endif

    /* Increase BH mass by measured accretion rate * dt */
    ac.mbh += ac.accr_rate * g_dt;

    /* Update Eddington luminosity - calculate after increasing BH mass */
    ac.edd = EddingtonLuminosity(ac.mbh);


}


/*********************************************************************** */
void SphericalAccretionOutput() {
/*!
 * Perform output of spherical accretion calculations
 *
 *********************************************************************** */

    if (prank == 0) {

        char *dir, fname[512];
        FILE *fp_acc;
        static double next_output = -1;
        double accr_rate_msun_yr;
        double accr_rate_bondi_msun_yr;

        dir = GetOutputDir();
        sprintf (fname, "%s/accretion_rate.dat", dir);

        /* Open file if first timestep (but not restart).
         * We always write out the first timestep. */
        if (g_stepNumber == 0) {
            fp_acc = fopen(fname, "w");
        }

            /* Prepare for restart or append if this is not step 0  */
        else {

            /* In case of restart, get last timestamp
             * and determine next timestamp */
            if (next_output < 0.0) {
                char sline[512];
                fp_acc = fopen(fname, "r");
                while (fgets(sline, 512, fp_acc)) { }
                sscanf(sline, "%lf\n", &next_output);
                next_output += ACCRETION_OUTPUT_RATE + 1;
                fclose(fp_acc);
            }

            /* Append if next output step has been reacehd */
            if (g_time > next_output) fp_acc = fopen(fname, "a");

        }

        /* Write data */
        if (g_time > next_output) {

            /* Accretion rate in cgs */
            accr_rate_msun_yr = ac.accr_rate * vn.mdot_norm / (CONST_Msun / (CONST_ly / CONST_c));

#if SINK_METHOD == SINK_BONDI
            accr_rate_bondi_msun_yr = ac.accr_rate_bondi * vn.mdot_norm / (CONST_Msun / (CONST_ly / CONST_c));
            fprintf(fp_acc, "%12.6e  %12.6e  %12.6e %12.6e %12.6e \n", g_time * vn.t_norm / (CONST_ly / CONST_c),
                    g_dt * vn.t_norm / (CONST_ly / CONST_c), accr_rate_msun_yr, accr_rate_bondi_msun_yr,
                    ac.mbh * vn.m_norm / CONST_Msun);

#else
            fprintf(fp_acc, "%12.6e  %12.6e  %12.6e %12.6e \n", g_time * vn.t_norm / (CONST_ly / CONST_c),
                    g_dt * vn.t_norm / (CONST_ly / CONST_c), accr_rate_msun_yr, ac.mbh * vn.m_norm / CONST_Msun);
#endif

            next_output += ACCRETION_OUTPUT_RATE;
            fclose(fp_acc);

        }

    }

}



/* ************************************************ */
void SphericalFreeflow(double *prims, double ****Vc, const double *x1, const double *x2, const double *x3,
                         const int k, const int j, const int i) {
/*!
 * Applies the spherical equivalent of free flowing boundary condition for cell at x1[i], x2[j], x3[k].
 * Assumes sphere center is at (0, 0, 0).
 *
 ************************************************** */

    /* Cental cell radius */
    double radc = SPH1(x1[i], x2[j], x3[k]);

    /* Number of cells in a direction to look at - hardcoded here */
    int e = 1;

    /* Loop over all surrounding cells and see if they're outside the radius of the central cell */
    double rad;
    double weight, weights;
    int kk, jj, ii, nv;

    weights = 0;

    /* Zero primitives array */
    for (nv = 0; nv < NVAR; nv++) prims[nv] = 0;

    for (kk = MAX(k - e, 0); kk < MIN(k + e + 1, NX3_TOT); kk++) {
        for (jj = MAX(j - e, 0); jj < MIN(j + e + 1, NX2_TOT); jj++) {
            for (ii = MAX(i - e, 0); ii < MIN(i + e + 1, NX1_TOT); ii++) {

                if (!((ii == i) && (jj == j) && (kk == k))) {

                    rad = SPH1(x1[ii], x2[jj], x3[kk]);
                    if (rad > radc) {
                        weight = 1. / rad;
                        weights += weight;

                        for (nv = 0; nv < NVAR; nv++) {
                            prims[nv] += weight * Vc[nv][kk][jj][ii];
                        }

                    }

                }
            }
        }
    }

    for (nv = 0; nv < NVAR; nv++) {
        if (weights > 0) {
            prims[nv] /= weights;
        }
        else {
            prims[nv] = Vc[nv][k][j][i];
        }
    }

}

/* ************************************************ */
double BondiAccretionRate(const double mbh, const double rho_far, const double snd_far){
/*!
 * Return Bondi accretion rate
 *
 ************************************************** */

    double lambda = BondiLambda();

    return  4 * CONST_PI * CONST_G * CONST_G * mbh * mbh * lambda * rho_far /
            (snd_far * snd_far * snd_far * vn.newton_norm * vn.newton_norm);

}

/* ************************************************ */
double BondiAccretionRateLocal(const double mbh, const double rho_acc, const double snd_acc, const double snd_far) {
/*!
 * Return Bondi accretion rate based only on the
 * sound speed far away and values of the sound speed
 * and the density at some point in the inflow.
 * It is the Bondi rate rewritten relating far-away density
 * to density and sound speed in accretion region.
 *
 ************************************************** */

    double lambda = BondiLambda();

    double rho_far;

    rho_far = pow(snd_far / snd_acc, 2 / (g_gamma - 1)) * rho_acc;

    return  4 * CONST_PI * CONST_G * CONST_G * mbh * mbh * lambda * rho_far /
            (snd_far * snd_far * snd_far * vn.newton_norm * vn.newton_norm);


}

/* ************************************************ */
double BondiLambda(){
/*!
 * Return the value for the constant lambda needed in the
 * Bondi accretion rate
 *
 ************************************************** */

    return pow(2., (9. - 7. * g_gamma) / (2. * (g_gamma - 1.))) *
           pow((5. - 3. * g_gamma), (5. - 3. * g_gamma) / (2. - 2. * g_gamma));

}

/* ************************************************ */
double EddingtonLuminosity(const double mbh) {
/*!
 * Return Eddington Luminosity for a given black hole mass.
 *
 ************************************************** */

    return 4 * CONST_PI * CONST_G * mbh * vn.m_norm * CONST_mH * CONST_c /
           (CONST_sigmaT * vn.power_norm);

}



/* ************************************************ */
void BondiFlowInternalBoundary(const double x1, const double x2, const double x3, double *result) {
/*!
 * Apply Bondi solution with accretion rate measured
 * We merely apply the solution for radii much smaller than the
 * critical radius, which reduces to the freefall solution.
 *
 * This routine is called from UserDefBoundary
 *
 ************************************************** */

    int nv;
    double sx1, sx2, sx3;
    double vs1, rho, prs, tmp, snd, snd_far, tmp_far;

    /* Radial velocity = Freefall velocity */
    sx1 = SPH1(x1, x2, x3);
    sx2 = SPH2(x1, x2, x3);
    sx3 = SPH3(x1, x2, x3);

    /* Decide if we are to take the measured rate or the theoretical Bondi accretion rate */
    if (ac.accr_rate < ac.accr_rate_bondi) rho = ac.accr_rate_bondi;
    else rho = ac.accr_rate;

    /* Solution for g_gamma < 5/3 */
    if (g_gamma < 5 / 3) {

        /* Flow speed */
        vs1 = -sqrt(2 * CONST_G * ac.mbh / (vn.newton_norm * sx1));

        /* "far-away" values */
        tmp_far = g_inputParam[PAR_HTMP] * ini_cgs[PAR_HTMP];
        snd_far = sqrt(g_gamma * CONST_kB * tmp_far / (MU_NORM * CONST_amu)) / vn.v_norm;
        tmp_far /= vn.temp_norm;

        /* Density, from conservation of mass */
        rho /= -(4 * CONST_PI * sx1 * sx1 * vs1);

        /* Pressure from temperature */
        tmp = 0.5 * tmp_far * CONST_G * ac.mbh / (vn.newton_norm * snd_far * snd_far * sx1);
        prs = PresIdealEOS(rho, tmp, MU_NORM);
    }

    /* Solution for g_gamma = 5/3 */
    else {

        /* Flow speed */
        vs1 = snd = -sqrt(CONST_G * ac.mbh / (2 * vn.newton_norm * sx1));

        /* Density, from conservation of mass */
        rho /= -(4 * CONST_PI * sx1 * sx1 * vs1);

        /* Pressure */
        prs = snd * snd * rho / g_gamma;
    }


    /* Copy quantities to primitives array */
    for (nv = 0; nv < NVAR; ++nv) result[nv] = 0;
    result[RHO] = rho;
    result[PRS] = prs;

    EXPAND(result[VX1] = VSPH_1(sx1, sx2, sx3, vs1, 0, 0);,
           result[VX2] = VSPH_2(sx1, sx2, sx3, vs1, 0, 0);,
           result[VX3] = VSPH_3(sx1, sx2, sx3, vs1, 0, 0););


}




/* ************************************************ */
void SphericalFreeflowInternalBoundary(const double ****Vc, int i, int j, int k, const double *x1, const double *x2,
                                       const double *x3, double *result) {
/*!
 * Apply Spherical inward freeflow internal boundary to cells around cell i, j, k.
 * This routine is called from UserDefBoundary
 *
 ************************************************** */

    int nv;
    double vx1, vx2, vx3, vs1;
    double sx1, sx2, sx3;
    double prims[NVAR];

    SphericalFreeflow(prims, Vc, x1, x2, x3, k, j, i);

    /* Remove non spherical velocity components of gas */
    vx1 = vx2 = vx3 = 0;
    EXPAND(vx1 = prims[VX1];,
           vx2 = prims[VX2];,
           vx3 = prims[VX3];);
    sx1 = SPH1(x1[i], x2[j], x3[k]);
    sx2 = SPH2(x1[i], x2[j], x3[k]);
    sx3 = SPH3(x1[i], x2[j], x3[k]);
    vs1 = VSPH1(x1[i], x2[j], x3[k], vx1, vx2, vx3);
    vs1 = MIN(vs1, 0);

    for (nv = 0; nv < NVAR; ++nv) result[nv] = prims[nv];

    EXPAND(result[VX1] = VSPH_1(sx1, sx2, sx3, vs1, 0, 0);,
           result[VX2] = VSPH_2(sx1, sx2, sx3, vs1, 0, 0);,
           result[VX3] = VSPH_3(sx1, sx2, sx3, vs1, 0, 0);
    );

}




/* ************************************************ */
void VacuumInternalBoundary(double *result) {
/*!
 * Apply Vacuum internal boundary conditions (as in Ruffert et al 1994)
 * This routine is called from UserDefBoundary
 *
 * Create a vaccuum accretor.
 * Using this method is not recommended,
 * as it drastically reduces the timestep.
 *
 ************************************************** */

    int nv;

    for (nv = 0; nv < NVAR; ++nv) result[nv] = 0;
    result[RHO] = g_smallDensity;
    result[PRS] = g_smallPressure;

}

