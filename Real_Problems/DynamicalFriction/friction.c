//
// Created by Alexander Y. Wagner on 2018/03/23.
//

#include "friction.h"
#include "pluto.h"
#include "pluto_usr.h"
#include "io_tools.h"
#include "init_tools.h"

GasdynamicalFriction gf;



void FrictionForce(const Data *d, Grid *grid) {

    /* Calculate force at (0, 0, 0) componentwise */

    int i, j, k;
    double fg, fg_x1, fg_x2, fg_x3;
    double ***rho;
    double *x1, *x2, *x3;

    /* State variables */
    rho = d->Vc[RHO];

    /* These are the geometrical central points */
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    double cell_vol, rad, rc;

    /* Core radius */
    rc = g_inputParam[PAR_RCORE];

    /* Total gravitational force in each direction */

    fg_x1 = fg_x2 = fg_x3 = 0;

    DOM_LOOP(k, j, i) {

                /* Gravitational force in cartesian coordinates*/
                rad = SPH1(x1[i], x2[j], x3[k]);
                cell_vol = ElevateVolume(grid->dV[k][j][i]);

                fg = -rho[k][j][i] * cell_vol / (rad * rad + rc * rc);
//                fg *= vn.force_norm;
                fg_x1 += VSPH2CART1(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x2 += VSPH2CART2(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x3 += VSPH2CART3(x1[i], x2[j], x3[k], fg, 0, 0);

            }

#ifdef PARALLEL

    MPI_Allreduce(&fg_x1, &gf.fx1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&fg_x2, &gf.fx2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&fg_x3, &gf.fx3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else // If not parallel

    gf.fx1 = fg_x1;
    gf.fx2 = fg_x2;
    gf.fx3 = fg_x3;

#endif // ifdef PARALLEL


}


void FrictionForceRadius(const Data *d, Grid *grid) {

    /* Calculate force at (0, 0, 0) in radial shells component-wise */

    int i, j, k;
    double fg, fg_x1, fg_x2, fg_x3;
    double ***rho;
    double *x1, *x2, *x3;

    /* State variables */
    rho = d->Vc[RHO];

    /* These are the geometrical central points */
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    double cell_vol, rad, rc;

    /* Core radius */
    rc = g_inputParam[PAR_RCORE];


    /* Local arrays */
    double * frx1_loc = ARRAY_1D(gf.nr, double);
    double * frx2_loc = ARRAY_1D(gf.nr, double);
    double * frx3_loc = ARRAY_1D(gf.nr, double);
    int ibin;

    for (ibin = 0; ibin < gf.nr; ibin++) frx1_loc[ibin] = frx2_loc[ibin] = frx3_loc[ibin] = 0.;

#if GEOMETRY == SPHERICAL
    IDOM_LOOP(i) {
        fg_x1 = fg_x2 = fg_x3 = 0;

        KDOM_LOOP(k){
            JDOM_LOOP(j) {

                /* Gravitational force in cartesian coordinates*/
                rad = SPH1(x1[i], x2[j], x3[k]);
                cell_vol = ElevateVolume(grid->dV[k][j][i]);

                fg = -rho[k][j][i] * cell_vol / (rad * rad + rc * rc);
//                fg *= vn.force_norm;
                fg_x1 += VSPH2CART1(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x2 += VSPH2CART2(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x3 += VSPH2CART3(x1[i], x2[j], x3[k], fg, 0, 0);

            }

        }

        /* Find bin and fill array */
        ibin = (int) (log10(rad / gf.rmin) / log10(gf.rmax / gf.rmin) * gf.nr);
        ibin = MIN(ibin, gf.nr - 1);

        frx1_loc[ibin] = fg_x1;
        frx2_loc[ibin] = fg_x2;
        frx3_loc[ibin] = fg_x3;
    }

    /* Reduce for each bin in r */
#ifdef PARALLEL

    MPI_Allreduce(frx1_loc, gf.frx1, gf.nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(frx2_loc, gf.frx2, gf.nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(frx3_loc, gf.frx3, gf.nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else // If not parallel

    for (ibin = 0; ibin < gf.nr; ibin++) {
        gf.frx1[ibin] = frx1_loc[ibin];
        gf.frx2[ibin] = frx2_loc[ibin];
        gf.frx3[ibin] = frx3_loc[ibin];
    }

#endif // ifdef PARALLEL

#else

    /* TODO: Not programmed yet - use random or uniform (?) sampling on sphere */

#endif

    FreeArray1D(frx1_loc);
    FreeArray1D(frx2_loc);
    FreeArray1D(frx3_loc);

}

void FrictionForceSlab(const Data *d, Grid *grid) {

    // TODO: Not programmed yet! Still the same as FrictionForceRadius

    /* Calculate force at (0, 0, 0) in radial shells component-wise */

    int i, j, k;
    double fg, fg_x1, fg_x2, fg_x3;
    double ***rho;
    double *x1, *x2, *x3;

    /* State variables */
    rho = d->Vc[RHO];

    /* These are the geometrical central points */
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    double cell_vol, rad, rc;

    /* Core radius */
    rc = g_inputParam[PAR_RCORE];


    /* Local arrays */
    double * fsx1_loc = ARRAY_1D(gf.nr, double);
    double * fsx2_loc = ARRAY_1D(gf.nr, double);
    double * fsx3_loc = ARRAY_1D(gf.nr, double);
    int ibin;

    for (ibin = 0; ibin < gf.nr; ibin++) fsx1_loc[ibin] = fsx2_loc[ibin] = fsx3_loc[ibin] = 0.;

#if GEOMETRY == SPHERICAL
    IDOM_LOOP(i) {
        fg_x1 = fg_x2 = fg_x3 = 0;

        KDOM_LOOP(k){
            JDOM_LOOP(j) {

                /* Gravitational force in cartesian coordinates*/
                rad = SPH1(x1[i], x2[j], x3[k]);
                cell_vol = ElevateVolume(grid->dV[k][j][i]);

                fg = -rho[k][j][i] * cell_vol / (rad * rad + rc * rc);
                fg_x1 += VSPH2CART1(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x2 += VSPH2CART2(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x3 += VSPH2CART3(x1[i], x2[j], x3[k], fg, 0, 0);

            }

        }

        /* Find bin and fill array */
        ibin = (int) (log10(rad / gf.rmin) / log10(gf.rmax / gf.rmin) * gf.nr);
        ibin = MIN(ibin, gf.nr - 1);

        fsx1_loc[ibin] = fg_x1;
        fsx2_loc[ibin] = fg_x2;
        fsx3_loc[ibin] = fg_x3;
    }

    /* Reduce for each bin in r */
#ifdef PARALLEL

    MPI_Allreduce(fsx1_loc, gf.frx1, gf.nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(fsx2_loc, gf.frx2, gf.nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(fsx3_loc, gf.frx3, gf.nr, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else // If not parallel

    for (ibin = 0; ibin < gf.nr; ibin++) {
        gf.frx1[ibin] = fsx1_loc[ibin];
        gf.frx2[ibin] = fsx2_loc[ibin];
        gf.frx3[ibin] = fsx3_loc[ibin];
    }

#endif // ifdef PARALLEL

#else

    /* TODO: Not programmed yet - use random or uniform (?) sampling on sphere */

#endif

    FreeArray1D(fsx1_loc);
    FreeArray1D(fsx2_loc);
    FreeArray1D(fsx3_loc);

// TODO: FrictionForceSlabOutput

}




void FrictionForcePolar(const Data *d, Grid *grid) {

    /* Calculate force at (0, 0, 0) in polar cones (theta) component-wise */

    int i, j, k;
    double fg, fg_x1, fg_x2, fg_x3;
    double ***rho;
    double *x1, *x2, *x3;

    /* State variables */
    rho = d->Vc[RHO];

    /* These are the geometrical central points */
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    double cell_vol, rad, theta, rc;

    /* Core radius */
    rc = g_inputParam[PAR_RCORE];


    /* Local arrays */
    double * fthx1_loc = ARRAY_1D(gf.nr, double);
    double * fthx2_loc = ARRAY_1D(gf.nr, double);
    double * fthx3_loc = ARRAY_1D(gf.nr, double);
    int ibin;

    for (ibin = 0; ibin < gf.nr; ibin++) fthx1_loc[ibin] = fthx2_loc[ibin] = fthx3_loc[ibin] = 0.;

#if GEOMETRY == SPHERICAL
    IDOM_LOOP(j) {
        fg_x1 = fg_x2 = fg_x3 = 0;

        KDOM_LOOP(i){
            JDOM_LOOP(k) {

                /* Gravitational force in cartesian coordinates*/
                rad = SPH1(x1[i], x2[j], x3[k]);
                cell_vol = ElevateVolume(grid->dV[k][j][i]);

                fg = -rho[k][j][i] * cell_vol / (rad * rad + rc * rc);
                fg_x1 += VSPH2CART1(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x2 += VSPH2CART2(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x3 += VSPH2CART3(x1[i], x2[j], x3[k], fg, 0, 0);

            }

        }

        /* Find bin and fill array */
        theta = SPH2(x1[i], x2[j], x3[k]);
        ibin = (int) ((theta - gf.thmin) / (gf.thmax - gf.thmin) * gf.nth);
        ibin = MIN(ibin, gf.nth - 1);

        fthx1_loc[ibin] = fg_x1;
        fthx2_loc[ibin] = fg_x2;
        fthx3_loc[ibin] = fg_x3;

    }

    // TODO: The reduction here is not correct
    /* Reduce for each bin in r */
#ifdef PARALLEL

    MPI_Allreduce(fthx1_loc, gf.fthx1, gf.nth, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(fthx2_loc, gf.fthx2, gf.nth, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(fthx3_loc, gf.fthx3, gf.nth, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else // If not parallel

    for (ibin = 0; ibin < gf.nth; ibin++) {
        gf.fthx1[ibin] = fthx1_loc[ibin];
        gf.fthx2[ibin] = fthx2_loc[ibin];
        gf.fthx3[ibin] = fthx3_loc[ibin];
    }

#endif // ifdef PARALLEL

#else

    /* TODO: Not programmed yet - use random or uniform (?) sampling on sphere */

#endif

    FreeArray1D(fthx1_loc);
    FreeArray1D(fthx2_loc);
    FreeArray1D(fthx3_loc);

}


void SetFriction(const Grid *grid) {/* Make arrays for force vs radius and force vs polar angle */

    if (gf.r == NULL){

        OptimalRadialBinning(grid, &gf.rmin, &gf.rmax, &gf.nr);

        OptimalPolarBinning(grid, &gf.thmin, &gf.thmax, &gf.nth);

        /* Initialize arrays */
        gf.r = ARRAY_1D(gf.nr, double);
        gf.dr = ARRAY_1D(gf.nr, double);
        gf.frx1 = ARRAY_1D(gf.nr, double);
        gf.frx2 = ARRAY_1D(gf.nr, double);
        gf.frx3 = ARRAY_1D(gf.nr, double);

        gf.th = ARRAY_1D(gf.nth, double);
        gf.fthx1 = ARRAY_1D(gf.nth, double);
        gf.fthx2 = ARRAY_1D(gf.nth, double);
        gf.fthx3 = ARRAY_1D(gf.nth, double);

        MakeLogRadiusArray(gf.r, gf.dr, &gf.dy, gf.rmin, gf.rmax, gf.nr);

        MakeUniformPolarArray(gf.th, &gf.dth, gf.thmin, gf.thmax, gf.nth);

    }
}


void MakeUniformPolarArray(double *th, double *dth, double thmin, double thmax, int nth) {

    /* Create uniform theta array */
    for (int i = 0; i < nth; i++) {
            th[i] = thmin + (i + 0.5) / nth * (thmax - thmin) ;
        }

    *dth = (thmax - thmin) / nth;

}


void MakeLogRadiusArray(double *r, double *dr, double *dy, double rmin, double rmax, int nr) {

    /* Create logarithmic radial array */

    /* Uniform spacing in log space */
    *dy  = log10(rmax / rmin);
    *dy /= (double) (nr);

    /* Real space points */
    double xl, xr, dx;
    xl = rmin;
    for (int i = 0; i < nr; i++) {
            dr[i] = xl * (pow(10.0, *dy) - 1.0);
            xr = xl + dr[i];
            r[i] = 0.5 * (xl + xr);
            if (i < nr - 1) xl = xr;
        }
}


void OptimalRadialBinning(const Grid *grid, double *rmin, double *rmax, int *nr) {

    int nx1_glob = grid->np_int_glob[IDIR];
    int nx2_glob = grid->np_int_glob[JDIR];
    int nx3_glob = grid->np_int_glob[KDIR];

#if GEOMETRY == SPHERICAL
    *nr = nx1_glob;
    *rmin = MAX(g_domBeg[IDIR], grid->dl_min[IDIR]);
    *rmax = g_domEnd[IDIR];

/* Untested */
#elif GEOMETRY == POLAR
    *nr = MIN(nx1_glob / 2, nx3_glob / 2);
    *rmin = MAX(g_domBeg[IDIR], grid->dl_min[IDIR]);
    *rmax = MIN(fabs(g_domBeg[KDIR]), g_domEnd[KDIR]);
    *rmax = MIN(*rmax, g_domEnd[IDIR]);

/* Untested */
#else
    *nr = MIN(nx1_glob, D_SELECT(nx1_glob, nx2_glob, nx3_glob));
    *rmin = grid->dl_min[IDIR];
    D_EXPAND(*rmax =            MIN(fabs(g_domBeg[IDIR]), g_domEnd[IDIR]);,
             *rmax = MIN(*rmax, MIN(fabs(g_domBeg[JDIR]), g_domEnd[JDIR]));,
             *rmax = MIN(*rmax, MIN(fabs(g_domBeg[KDIR]), g_domEnd[KDIR])););

#endif


}


void OptimalPolarBinning(const Grid *grid, double *thmin, double *thmax, int *nth) {

    int nx1_glob = grid->np_int_glob[IDIR];
    int nx2_glob = grid->np_int_glob[JDIR];
    int nx3_glob = grid->np_int_glob[KDIR];

#if GEOMETRY == SPHERICAL
    *nth = nx2_glob;

/* Untested */
#elif GEOMETRY == POLAR
    *nth = MIN(nx1_glob / 2, nx3_glob / 2);
    *nth = MIN(*nth, nx2_glob / 2);

/* Untested */
#else
    *nth = MIN(nx1_glob, D_SELECT(nx1_glob, nx2_glob, nx3_glob));

#endif

    *thmin = 0;
    *thmax = CONST_PI;

}


void FrictionForceOutput() {

    if (prank == 0) {

        FILE *fp_out;
        char fname[512];
        sprintf(fname, "total_force.dat");
        static double next_output = -1;

        next_output = OutputContextEnter(fname, &fp_out, next_output, 0.3);

        /* Write data */
        if (g_time > next_output) {

            fprintf(fp_out, "%12.6e %12.6e %12.6e %12.6e \n",
                    g_time ,                       // time (code units)
                    gf.fx1 ,                       // Total force x1 component
                    gf.fx2 ,                       // Total force x2 component
                    gf.fx3                         // Total force x3 component
            );

        }

        next_output = OutputContextExit(&fp_out, next_output, 0.3);

    }
}

void FrictionForceRadiusOutput() {

    if (prank == 0) {

        FILE *fp_out;
        char fname[512];
        sprintf(fname, "radius_force.dat");
        static double next_output = -1;

        next_output = OutputContextEnter(fname, &fp_out, next_output, 0.3);

        /* Write data */
        if (g_time > next_output) {

            fprintf(fp_out, "%12.6e %12d %12d %12d %12d \n",
                    g_time,  gf.nr, 5, 0, 0        // time (code units), nlines, ncols
            );

            for (int i = 0; i < gf.nr; i++)
            fprintf(fp_out, "%12.6e %12.6e %12.6e %12.6e %12.6e \n",
                    gf.r[i]       ,                // ith bin radius
                    gf.dr[i]      ,                // ith bin width
                    gf.frx1[i]    ,                // Total force x1 component at r[i]
                    gf.frx2[i]    ,                // Total force x2 component at r[i]
                    gf.frx3[i]                     // Total force x3 component at r[i]
            );
        }

        next_output = OutputContextExit(&fp_out, next_output, 0.3);

    }
}


void FrictionForcePolarOutput() {

    if (prank == 0) {

        FILE *fp_out;
        char fname[512];
        sprintf(fname, "polar_force.dat");
        static double next_output = -1;

        next_output = OutputContextEnter(fname, &fp_out, next_output, 0.3);

        /* Write data */
        if (g_time > next_output) {

            fprintf(fp_out, "%12.6e %12d %12d %12d %12d \n",
                    g_time,  gf.nth, 5, 0, 0        // time (code units), nlines, ncols
            );

            for (int i = 0; i < gf.nth; i++)
                fprintf(fp_out, "%12.6e %12.6e %12.6e %12.6e %12.6e \n",
                        gf.th[i]       ,           // ith bin polar angle
                        gf.dth         ,           // ith bin width
                        gf.fthx1[i]    ,           // Total force x1 component at r[i]
                        gf.fthx2[i]    ,           // Total force x2 component at r[i]
                        gf.fthx3[i]                // Total force x3 component at r[i]
                );
        }

        next_output = OutputContextExit(&fp_out, next_output, 0.3);

    }
}
