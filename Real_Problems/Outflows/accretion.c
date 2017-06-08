//
// Created by Alexander Y. Wagner on 4/6/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "accretion.h"
#include "grid_geometry.h"
#include "idealEOS.h"
#include "interpolation.h"
#include "outflow.h"
#include "init_tools.h"
#include "hot_halo.h"


/* Global struct for accretion */
Accretion ac;


/* ************************************************ */
void SetAccretionPhysics() {
/*
 * Set values for the accretion structure
 *
 * SetAccretionPhysics is only called once
 * at the beginning of initialization.
 *
 ************************************************** */

    // TODO: Some of these quantities need to be uptdated after restart.

    /* Quantities from input parameters */
    ac.rad = g_inputParam[PAR_ARAD] * ini_code[PAR_ARAD];
    ac.mbh = g_inputParam[PAR_AMBH] * ini_code[PAR_AMBH];
    ac.eff = g_inputParam[PAR_AEFF] * ini_code[PAR_AEFF];
    ac.mld = g_inputParam[PAR_AMLD] * ini_code[PAR_AMLD];
    ac.snk = g_inputParam[PAR_ASNK] * ini_code[PAR_ASNK];

    /* Bondi Accretion parameters */
    double snd_acc = sqrt(g_gamma * g_smallPressure / g_smallDensity);
    ac.accr_rate_bondi = BondiAccretionRate(ac.mbh, g_smallDensity, snd_acc);

    /* Eddington rate */
    ac.edd = EddingtonLuminosity(ac.mbh);

    /* Global accretion rate */
    ac.accr_rate = 0;

    /* Area of accretion surface. */
    ac.area = 4 * CONST_PI * ac.rad * ac.rad;

    /* Number of cells in the spherical surface of radius ac.rad is
     * calculated in Analysis analysis */

    /* Set deboost to 1 to trigger reference value calculations in
     * SetNozzleGeometry and SetOutflowState (used only in FEEDBACK_CYCLE mode) */
    ac.deboost = 1;
    os.is_on = 1;

    /* Set outflow geometry struct with parameters of cone */
    SetNozzleGeometry(&(ac.nzi));

    /* Set outflow physics */
    SetOutflowState(&(ac.osi));

    /* Initial value deboosting factor for outflow (used only in FEEDBACK_CYCLE mode) */
    ac.deboost = 0;
    os.is_on = 0;

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

    /* Sink is be deactivated if ASNK is zero */
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
    // TODO: make consistent with Nozzle (is_two_sided)
#if INTERNAL_BOUNDARY == YES
    D_SELECT(cx1p, cx2p, cx3p) = fabs(D_SELECT(cx1p, cx2p, cx3p));
#endif

    if (nz.is_fan) {

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

    /* Measure accretion rate at a given radius */
    ac.accr_rate = SphericalSampledAccretion(d, grid, ac.rad);

#if MEASURE_BONDI_ACCRETION
    ac.accr_rate_bondi = BondiAccretion(d, grid);
#endif

    /* Increase BH mass by measured accretion rate * dt */
    ac.mbh += ac.accr_rate * g_dt * (1. - ac.eff) / (1. + ac.mld);

    /* Update Eddington luminosity - calculate after increasing BH mass */
    ac.edd = EddingtonLuminosity(ac.mbh);

}



/* ************************************************ */
double SphericalSampledAccretion(const Data *d, Grid *grid, const double radius) {

/*!
 * Calculate the spherical accretion rate through the surface of a sphere.
 *
 ************************************************** */

    double rho, vs1;
    double vx1, vx2, vx3;
    double *x1, *x2, *x3;
    double accr, accr_rate = 0, accr_rate_rss;

    /* Determine number of sampling points to use */
    int npoints;
    double oversample = 30;
    double area, area_per_point;
    double dl_min = grid[IDIR].dl_min;

    /* Area through which accretion rate is measured.
     * If nozzle is not two-sided, the accretion rate represents
     * the accretion rate in one half of the galaxy. */
    area = 4 * CONST_PI * radius * radius;
    if (!nz.is_two_sided) area /= 2.;

    npoints = (int) (area / (dl_min * dl_min) * oversample);
#if DIMENSIONS == 2
    npoints = sqrt(npoints);
#endif

    /* Get npoint points on spherical surface */
    x1 = ARRAY_1D(npoints, double);
    x2 = ARRAY_1D(npoints, double);
    x3 = ARRAY_1D(npoints, double);
    npoints = UniformSamplingSphericalSurface(npoints, radius, x1, x2, x3);
    area_per_point = area / npoints;

    /* Determine which variables to interpolate */
    double v[NVAR];
    int vars[] = {RHO, ARG_EXPAND(VX1, VX2, VX3), -1};

    /* Calculate accretion rate at every point, if it is in the local domain */
    for (int ipoint = 0; ipoint < npoints; ipoint++){

        if( PointInDomain(grid, x1[ipoint], x2[ipoint], x3[ipoint]) ) {

            InterpolateGrid(d, grid, vars, x1[ipoint], x2[ipoint], x3[ipoint], v);

            /* Calculate and sum accretion rate */
            rho = v[RHO];
            vx1 = vx2 = vx3 = 0;
            EXPAND(vx1 = v[VX1];,
                   vx2 = v[VX2];,
                   vx3 = v[VX3];);
            vs1 = VSPH1(x1[ipoint], x2[ipoint], x3[ipoint], vx1, vx2, vx3);
            vs1 = fabs(MIN(vs1, 0));

            accr = rho * vs1 * area_per_point;
            accr_rate += accr;

        }

    }


#ifdef PARALLEL
    MPI_Allreduce(&accr_rate, &accr_rate_rss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    accr_rate_rss = accr_rate;

#endif

    FreeArray1D(x1);
    FreeArray1D(x2);
    FreeArray1D(x3);

    return accr_rate_rss;

}




/* ************************************************ */
double SphericalSelectedAccretion(const Data *d, Grid *grid, const double radius) {
/*!
 * Calculate the spherical accretion rate through the surface of a sphere.
 *
 * Deprecated - use SphericalSampledAccretion instead.
 *
 ************************************************** */


    /* Accretion - General variables */
    int i, j, k;
    double rho, vs1;
    double vx1, vx2, vx3;
    double *x1, *x2, *x3;
    double accr, accr_rate = 0, accr_rate_sel;
    double area, area_per_cell;
    int gcount = 0, lcount = 0;


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
                                                           radius);
#endif

#if SIC_METHOD == SIC_CORNERS || SIC_METHOD == SIC_HYBRID
                sicc = SphereSurfaceIntersectsCellByCorners(x1[i], x2[j], x3[k],
                                                            grid[IDIR].dx[i], grid[JDIR].dx[j], grid[KDIR].dx[k],
                                                            radius);
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

                    accr = rho * vs1;

                    /* Count twice in hybrid mode */
#if SIC_METHOD == SIC_HYBRID
                    if (sicr && sicc) accr *= 2;
#endif

                    accr_rate += accr;

                    /* Count the number of cells actually found */
                    lcount++;

                } // if sicr || sicc

            }  // DOM_LOOP

    /* Averaging in HYBRID mode */
#if SIC_METHOD == SIC_HYBRID
    accr_rate /= 2;

#endif // SIC_HYBRID


    /* Reductions, including MPI reductions,
     * and calculation of other quantities. */

#ifdef PARALLEL
    MPI_Allreduce(&accr_rate, &accr_rate_sel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&lcount, &gcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


#else // If not parallel

    accr_rate_sel = accr_rate;
    gcount = lcount;

#endif // ifdef PARALLEL

    area = 4 * CONST_PI * radius * radius;
    area_per_cell = area / gcount;

    accr_rate_sel *= area_per_cell;

    return accr_rate_sel;

}





/*********************************************************************** */
void FederrathAccretion(const Data *d, Grid *grid) {
/*!
 * Driver routine for removing mass in gas around central BH
 * according to Federrath accretion.
 *
 *********************************************************************** */

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



/*********************************************************************** */
void SphericalAccretionOutput() {
/*!
 * Perform output of spherical accretion calculations
 *
 *********************************************************************** */

    if (prank == 0) {

        char fname[512];
        FILE *fp_acc;
        static double next_output = -1;

        sprintf(fname, "accretion_rate.dat");

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

            /* Append if next output step has been reached */
            if (g_time > next_output) fp_acc = fopen(fname, "a");

        }


        /* Write data */
        if (g_time > next_output) {

            /* Accretion rate in cgs */
            double msun_yr = CONST_Msun / (CONST_ly / CONST_c);
            double accr_rate_msun_yr = ac.accr_rate * vn.mdot_norm / msun_yr;
            double accr_rate_bondi_msun_yr = ac.accr_rate_bondi * vn.mdot_norm / msun_yr;

            fprintf(fp_acc, "%12.6e  %12.6e  %12.6e %12.6e %12.6e %12.6e \n",
                    g_time * vn.t_norm / (CONST_ly / CONST_c),            // time
                    accr_rate_msun_yr,                                    // measured acc rate (flux)
                    accr_rate_bondi_msun_yr,                              // measured bondi accretion rate
                    ac.accr_rate * vn.mdot_norm * CONST_c * CONST_c,      // measured acc power
                    ac.mbh * vn.m_norm / CONST_Msun,                      // black hole mass
                    ac.edd * vn.power_norm);                              // Eddington power

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

    /* Central cell radius */
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
double BondiTurbulentAccretionRate(const double mbh, const double rho_far, const double snd_far, const double vorticity, const double rbondi){
/*!
 * Return Bondi accretion rate
 *
 ************************************************** */

    return  4 * CONST_PI * CONST_G * CONST_G * mbh * mbh * rho_far /
            (snd_far * snd_far * snd_far * vn.newton_norm * vn.newton_norm) *
            0.34 / (1. + pow(vorticity * rbondi / snd_far, 0.9));

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

    static int once = 0;
    static double lambda;

    if (!once) {

        lambda = pow(2., (9. - 7. * g_gamma) / (2. * (g_gamma - 1.))) *
                 pow((5. - 3. * g_gamma), (5. - 3. * g_gamma) / (2. - 2. * g_gamma));

        once = 1;

    }

    return lambda;

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
    double vx1, vx2, vx3, vs1, vs2, vs3;
    double sx1, sx2, sx3;
    double prims[NVAR];

    SphericalFreeflow(prims, Vc, x1, x2, x3, k, j, i);

    /* Limit radial velocity component to negative values. */
    vx1 = vx2 = vx3 = 0;
    EXPAND(vx1 = prims[VX1];,
           vx2 = prims[VX2];,
           vx3 = prims[VX3];);
    sx1 = SPH1(x1[i], x2[j], x3[k]);
    sx2 = SPH2(x1[i], x2[j], x3[k]);
    sx3 = SPH3(x1[i], x2[j], x3[k]);
    vs1 = VSPH1(x1[i], x2[j], x3[k], vx1, vx2, vx3);
    vs2 = VSPH2(x1[i], x2[j], x3[k], vx1, vx2, vx3);
    vs3 = VSPH3(x1[i], x2[j], x3[k], vx1, vx2, vx3);
    vs1 = MIN(vs1, 0);

    for (nv = 0; nv < NVAR; ++nv) result[nv] = prims[nv];

    EXPAND(result[VX1] = VSPH_1(sx1, sx2, sx3, vs1, vs2, vs3);,
           result[VX2] = VSPH_2(sx1, sx2, sx3, vs1, vs2, vs3);,
           result[VX3] = VSPH_3(sx1, sx2, sx3, vs1, vs2, vs3);
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


/* ************************************************ */
double JeansResolvedDensity(const double *prim){
/*!
 * Return the difference between the cell density and maximum density
 * obtained from the criterion that the Jeans length should be resolved
 * within the sink region is satisfied.
 *
 ************************************************** */

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

    EXPAND(double vx1 = prim[VX1];,
           double vx2 = prim[VX2];,
           double vx3 = prim[VX3];);

    double vmag = VMAG(x1, x2, x3, vx1, vx2, vx3);
    double sph1 = SPH1(x1, x2, x3);
    double rho = prim[RHO];
    double prs = prim[PRS];

    double ekin = 0.5 * mass * vmag * vmag;
    double epot = - CONST_G * mass * ac.mbh / (sph1 * vn.newton_norm);
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


/* ************************************************ */
double BondiRadius(double m, double *v){

/*!
 * Calculate Bondi Radius for a given mass and array of primitives.
 *
 ************************************************** */

    double rad, cs2;
    cs2 = g_gamma * v[PRS] / v[RHO];
    rad = CONST_G * m * vn.m_norm / cs2 / vn.l_norm;

    return rad;

}

/* ************************************************ */
double BondiAccretion(const Data *d, Grid *grid) {

/*!
 * Calculate standard average accretion rate within bondi radius.
 *
 ************************************************** */

    double rho, prs, tmp, snd;
    double *x1, *x2, *x3, rad, bondi_radius;
    double accr = 0, accr_g;
    int lcount = 0, gcount;
    double halo[NVAR];

    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    /* Get Bondi radius */
    tmp = g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP];
    halo[RHO] = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
    EXPAND(halo[VX1] = 0;,
           halo[VX2] = 0;,
           halo[VX3] = 0;);
    halo[PRS] = halo[RHO] * tmp / MU_NORM;
    halo[TRC] = 0;
    bondi_radius = BondiRadius(ac.mbh, halo);

#if TURBULENT_BONDI_ACCRETION
    double vort1, vort2, vort3, vorticity, accr_vort;
    EXPAND(double ***vx1 = d->Vc[VX1];,
           double ***vx2 = d->Vc[VX2];,
           double ***vx3 = d->Vc[VX3];);
#endif

    int i, j, k;
    DOM_LOOP(k, j, i) {

                /* Bondi accretion rate (smooth) within Bondi radius */
                rad = SPH1(x1[i], x2[j], x3[k]);
                if (rad < bondi_radius) {

                    rho = d->Vc[RHO][k][j][i];
                    prs = d->Vc[PRS][k][j][i];
                    snd = sqrt(g_gamma * prs / rho);

                    accr = BondiAccretionRate(ac.mbh, rho, snd);

#if TURBULENT_BONDI_ACCRETION

                    /* Curl operation. */
                    // TODO: Turn this into a function (or macro)
                    EXPAND(,
                           vort3 = CDIFF_X1(vx2, k, j, i) - CDIFF_X2(vx1, k, j, i);,
                           vort2 = CDIFF_X3(vx1, k, j, i) - CDIFF_X1(vx3, k, j, i);
                           vort1 = CDIFF_X2(vx3, k, j, i) - CDIFF_X3(vx2, k, j, i););

                    vorticity = sqrt(EXPAND(0, + vort3 * vort3,  + vort2 * vort2 + vort1 * vort1));

                    accr_vort = BondiTurbulentAccretionRate(ac.mbh, rho, snd, vorticity, bondi_radius);

                    accr = 1. / sqrt(1. / (accr * accr) + 1. / (accr_vort * accr_vort));

#endif

                    lcount++;

                }
            }


#if PARALLEL
    /* MPI reduce Bondi accretion parameters */
    MPI_Allreduce(&lcount, &gcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&accr, &accr_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    accr_g /= gcount;

#else

    accr_g = accr / lcount;

#endif

    return accr_g;

}


