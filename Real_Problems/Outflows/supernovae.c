//
// Created by Alexander Y. Wagner on 6/21/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "supernovae.h"
#include "init_tools.h"


/* Global struct for accretion */
Supernovae sn;


void SetSupernovaePhysics(){

    sn.pow = g_inputParam[PAR_SPOW] * ini_code[PAR_SPOW];
    sn.dur = g_inputParam[PAR_SDUR] * ini_code[PAR_SDUR];
    sn.rad = g_inputParam[PAR_SRAD] * ini_code[PAR_SRAD];
    sn.hgt = g_inputParam[PAR_SHGT] * ini_code[PAR_SHGT];
    sn.en1 = SUPERNOVA_ENERGY / vn.eint_norm;
    sn.ene = sn.pow * sn.dur;
    sn.num_tot = sn.ene / sn.en1 + 1;
    sn.num_left = sn.num_tot;
    sn.num_dt = 0;

}


/* ****************************************************** */
void StrewSupernovae() {
/*!
 *  Create random SN locations and for a random number of SN
 *
 * ****************************************************** */

    if (prank == 0) {
        /* Inverse probability */
        double nsteps_left = (sn.dur - g_time) / g_dt;

        /* Case of SN */
        sn.num_dt = 0;
        while ((RandomNumber(0, sn.num_left / MAX(nsteps_left, 1) + 1) > 1) &&
               (sn.num_dt < NSTARS_MAX)) {

            /* SN location in cylindrical polar coords.
             * This will not give a uniform distribution in cylindrical coords */
            double sn_rad = RandomNumber(0, sn.rad);
            double sn_hgt = RandomNumber(-sn.hgt, sn.hgt);
            double sn_phi = RandomNumber(0, 2 * CONST_PI);

            /* Convert to current coordinate system */
            double xcur1 = POL_1(sn_rad, sn_phi, sn_hgt);
            double xcur2 = POL_2(sn_rad, sn_phi, sn_hgt);
            double xcur3 = POL_3(sn_rad, sn_phi, sn_hgt);
            sn.x1[sn.num_dt] = xcur1;
            sn.x2[sn.num_dt] = D_SELECT(0, xcur3, xcur2);
            sn.x3[sn.num_dt] = D_SELECT(0, 0, xcur3);


            /* Count SN */
            sn.num_dt += 1;
            sn.num_left -= 1;

        }

    }

#ifdef PARALLEL
    MPI_Bcast(sn.x1, NSTARS_MAX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sn.x2, NSTARS_MAX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sn.x3, NSTARS_MAX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sn.num_dt, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sn.num_left, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

}


/* ****************************************************** */
void InjectSupernovae(double *result, const double x1, const double x2, const double x3, double sn_inj) {
/*!
 *   Inject supernova energy cellwise by increasing pressure.
 *
 * ****************************************************** */

    for (int isn = 0; isn < sn.num_dt; isn++) {

        double r = SPH1(x1 - sn.x1[isn], x2 - sn.x2[isn], x3 - sn.x3[isn]);

        /* Inject supernova(e) in region around (0,0,0) within a given radius */
        if (r < sn_inj) {

            double vol = D_EXPAND(2. * sn_inj, / 2. * sn_inj * CONST_PI, *4. / 3. * sn_inj);

            result[PRS] += sn.en1 / vol * (g_gamma - 1.);
            result[TRC] = 1;
        }

    }
}