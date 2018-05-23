//
// Created by Alexander Y. Wagner on 2018/03/23.
//

#include "friction.h"
#include "pluto.h"
#include "pluto_usr.h"
#include "io_tools.h"

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
    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    double cell_vol, rad;
    double *dV1, *dV2, *dV3;
    dV1 = grid[IDIR].dV;
    dV2 = grid[JDIR].dV;
    dV3 = grid[KDIR].dV;

    fg_x1 = fg_x2 = fg_x3 = 0;


    /* Make arrays for force vs radius and force vs polar angle */

    if (gf.r == NULL){

        /* Number of points for r and theta arrays. For POLAR, assume 3D setup.
         * For CARTESIAN, assume NX1 = NX2 (flow direction is KDIR). */

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

        MakeLogRadiusArray(gf.r, gf.dr, gf.rmin, gf.rmax, gf.nr);

        MakeUniformPolarArray(gf.th, &gf.dth, gf.thmin, gf.thmax, gf.nth);

    }


    /* Total gravitational force in each direction */

    DOM_LOOP(k, j, i) {

                /* Gravitational force in cartesian coordinates*/
                rad = SPH1(x1[i], x2[j], x3[k]);
                cell_vol = dV1[i] * dV2[j] * dV3[k];

                fg = -rho[k][j][i] * cell_vol / (rad * rad);
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


    /* TODO: Extract the force as functions of radius and polar angle as two separate functions */

    /* Gravitational force as a function of radius */

#if GEOMETRY == SPHERICAL
    IDOM_LOOP(i) {
        KDOM_LOOP(k){
            JDOM_LOOP(j) {

                /* Gravitational force in cartesian coordinates*/
                rad = SPH1(x1[i], x2[j], x3[k]);
                cell_vol = dV1[i] * dV2[j] * dV3[k];

                fg = -rho[k][j][i] * cell_vol / (rad * rad);
                fg_x1 += VSPH2CART1(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x2 += VSPH2CART2(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x3 += VSPH2CART3(x1[i], x2[j], x3[k], fg, 0, 0);

            }

        }

        /* Reduce for each r */
#ifdef PARALLEL

        MPI_Allreduce(&fg_x1, &gf.frx1[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&fg_x2, &gf.frx2[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&fg_x3, &gf.frx3[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else // If not parallel

        gf.frx1[i] = fg_x1;
        gf.frx2[i] = fg_x2;
        gf.frx3[i] = fg_x3;

#endif // ifdef PARALLEL

    }
#else

    /* TODO: Not programmed yet - use random or uniform (?) sampling on sphere */

#endif



    /* Gravitational force as a function of theta */

#if GEOMETRY == SPHERICAL
    JDOM_LOOP(j) {
        KDOM_LOOP(k){
            IDOM_LOOP(i) {

                /* Gravitational force in cartesian coordinates*/
                rad = SPH1(x1[i], x2[j], x3[k]);
                cell_vol = dV1[i] * dV2[j] * dV3[k];

                fg = -rho[k][j][i] * cell_vol / (rad * rad);
                fg_x1 += VSPH2CART1(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x2 += VSPH2CART2(x1[i], x2[j], x3[k], fg, 0, 0);
                fg_x3 += VSPH2CART3(x1[i], x2[j], x3[k], fg, 0, 0);

            }

        }

        /* Reduce for each r */
#ifdef PARALLEL

        MPI_Allreduce(&fg_x1, &gf.fthx1[j], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&fg_x2, &gf.fthx2[j], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&fg_x3, &gf.fthx3[j], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else // If not parallel

        gf.fthx1[j] = fg_x1;
        gf.fthx2[j] = fg_x2;
        gf.fthx3[j] = fg_x3;

#endif // ifdef PARALLEL

    }
#else

    /* TODO: Not programmed yet - use random or uniform (?) sampling on circle */

#endif



}


void MakeUniformPolarArray(double *th, double *dth, double thmin, double thmax, int nth) {

    /* Create uniform theta array */
    for (int i = 0; i < nth; i++) {
            th[i] = thmin + (i + 0.5) / nth * (thmax - thmin) ;
        }

    *dth = (thmax - thmin) / nth;

}


void MakeLogRadiusArray(double *r, double *dr, double rmin, double rmax, int nr) {

    /* Create logarithmic radial array, and uniform theta array */

    /* Uniform spacing in log space */
    double dy;
    dy  = log10(rmax / rmin);
    dy /= (double) (nr);

    /* Real space points */
    double xl, xr, dx;
    xl = rmin;
    for (int i = 0; i < nr; i++) {
            dr[i] = xl * (pow(10.0, dy) - 1.0);
            xr = xl + dr[i];
            r[i] = 0.5 * (xl + xr);
            if (i < nr - 1) xl = xr;
        }
}


void OptimalRadialBinning(const Grid *grid, double *rmin, double *rmax, int *nr) {

#if GEOMETRY == SPHERICAL
    *nr = NX1;
    *rmin = MAX(g_domBeg[IDIR], grid[IDIR].dl_min);
    *rmax = g_domEnd[IDIR];

/* Untested */
#elif GEOMETRY == POLAR
    *nr = MIN(NX1 / 2, NX3 / 2);
    *rmin = MAX(g_domBeg[IDIR], grid[IDIR].dl_min);
    *rmax = MIN(fabs(g_domBeg[KDIR]), g_domEnd[KDIR]);
    *rmax = MIN(*rmax, g_domEnd[IDIR]);

/* Untested */
#else
    *nr = MIN(NX1, D_SELECT(NX1, NX2, NX3));
    *rmin = grid[IDIR].dl_min;
    D_EXPAND(*rmax =            MIN(fabs(g_domBeg[IDIR]), g_domEnd[IDIR]);,
             *rmax = MIN(*rmax, MIN(fabs(g_domBeg[JDIR]), g_domEnd[JDIR]));,
             *rmax = MIN(*rmax, MIN(fabs(g_domBeg[KDIR]), g_domEnd[KDIR])););

#endif


}


void OptimalPolarBinning(const Grid *grid, double *thmin, double *thmax, int *nth) {

#if GEOMETRY == SPHERICAL
    *nth = NX2;

/* Untested */
#elif GEOMETRY == POLAR
    *nth = MIN(NX1 / 2, NX3 / 2);
    *nth = MIN(*nth, NX2 / 2);

/* Untested */
#else
    *nth = MIN(NX1, D_SELECT(NX1, NX2, NX3));

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

            fprintf(fp_out, "%12.6e ",
                    g_time                         // time (code units)
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
        sprintf(fname, "palar_force.dat");
        static double next_output = -1;

        next_output = OutputContextEnter(fname, &fp_out, next_output, 0.3);

        /* Write data */
        if (g_time > next_output) {

            fprintf(fp_out, "%12.6e ",
                    g_time                         // time (code units)
            );

            for (int i = 0; i < gf.nr; i++)
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
