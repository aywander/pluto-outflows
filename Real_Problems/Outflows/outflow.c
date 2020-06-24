//
// Created by Alexander Y. Wagner on 4/7/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "outflow.h"
#include "grid_geometry.h"
#include "hot_halo.h"
#include "accretion.h"
#include "init_tools.h"
#include "idealEOS.h"
#include "io_tools.h"
#include "nozzle.h"

/* Global struct for outflow parameters */
OutflowState os;



/* ************************************************ */
void SetOutflowState(OutflowState *ofs) {
/*
 * Set outflow state.
 * Set outflow primitives/parameters structs.
 *
 ************************************************** */

    // TODO: Make sure OutflowPrimitives is not called in a DOM/TOT_LOOP
    // TODO: Change g_inputParam[PAR_OPOW] -> nz.pow, etc, everywhere in code.

    // TODO: Some of these quantities need to be uptdated after restart.

    // TODO: allow the possibilty to set fluxes directly if outflow is aligned with grid.

    // TODO: Create freely expanding, self-similar solution inside nozzle (for improved interpolation).

    NOZZLE_SELECT(SetJetState(ofs), SetUfoState(ofs));

}


/* ************************************************ */
void SetJetState(OutflowState *ofs) {
/*
 * Set outflow state.
 * Set outflow primitives/parameters structs.
 *
 ************************************************** */


    double power, chi, lorentz;
    double speed, dens, pres, gmm1;
    int is_on;

    /* The parameters from ini, and normalize them */
    chi = g_inputParam[PAR_OMDT] * ini_code[PAR_OMDT];
    lorentz = g_inputParam[PAR_OSPD] * ini_code[PAR_OSPD];

    // TODO: Test accretion cycle
#if FEEDBACK_CYCLE == YES

    /* Jet coupling efficiency */
    double oeff = g_inputParam[PAR_OEFF] * ini_code[PAR_OEFF];

    /* AGN power from accretion rate */
    power = (1. - ac.eff) / (1. + ac.mld) * oeff * ac.eff * ac.accr_rate * CONST_c * CONST_c / vn.pot_norm;

    /* Hot phase pressure */
    double hrho = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
    double htmp = g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP];
    double hprs = PresIdealEOS(hrho, htmp, MU_NORM);

    /* Maximum thermal pressure of AGN outflow */
    /* For the jet, we always use the input lorentz factor at the nozzle */
    speed = Lorentz2Speed(lorentz);
    gmm1 = g_gamma / (g_gamma - 1.);
    pres = power / (gmm1 * lorentz * lorentz * speed * ac.nzi.area);

    /* This is for initialization of reference outflow state struct contained in the ac struct,
     * through SetAccretionPhysics at the beginning of Init. It does not affect dynamics of first timestep. */
    if (g_time == 0 && ac.deboost == 1) {
        power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
        speed = Lorentz2Speed(lorentz);
        gmm1 = g_gamma / (g_gamma - 1.);
        pres = power / (gmm1 * lorentz * lorentz * speed * nz.area * (1. + (lorentz - 1.) / lorentz * chi));
        dens = chi * gmm1 * pres;
        is_on = 1;
    }

    /* If time = 0, or insufficient power for pressure equilibrium, shut off AGN */
    else if ((pres < hprs) || (g_time == 0 )) {
        ac.deboost = 0;
        dens = g_smallDensity;
        pres = g_smallPressure;
        lorentz = 0;
        chi = 0;
        power = 0;
        is_on = 0;
    }

        /* If sufficient power, AGN is on */
    else {

        /* Set pressure, density, and speed, adjusting parameters such that
         * pressure is be at least ambient pressure, hprs. */

        /* If sufficient power ... */
        pres = power / (gmm1 * lorentz * lorentz * speed * ac.nzi.area * (1. + (lorentz - 1.) / lorentz * chi));
        dens = chi * gmm1 * pres;
        is_on = 1;

#if FBC_DEBOOST == YES
        /* ...but if insufficient heat with default Ekin, reduce (deboost) ...? */
        if (pres < hprs) {

            /* Ensure pressure equilibrium */
            pres = hprs;

#if FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_0

            /* Deboost area only. Density stays the same. */

            ac.deboost = power / (gmm1 * pres * lorentz * lorentz * speed * ac.nzi.area *
                         (1. + (lorentz - 1.) / lorentz * chi));

#elif FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_1

            /* Deboost chi. Density decreases accordingly. */

            ac.deboost = (power / (gmm1 * pres * ac.nzi.area * speed * lorentz * lorentz) - 1) /
                         ((lorentz - 1) / lorentz * chi)

            dens *= ac.deboost

#elif FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_2

            /* Deboost area and chi.  */

            ac.deboost = (sqrt(1 + 4 * power * chi * (lorentz - 1) /
                         (gmm1 * pres * ac.nzi.area * speed * lorentz * lorentz * lorentz)) - 1) /
                         (2 * chi * (lorentz - 1) / lorentz)

            dens *= ac.deboost;

#elif FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_3

            /* Deboost speed
             *
             * The expression for this is too complicated.
             * */

            print("Error in SetJetState for Feedback cycle FBC_DEBOOST_MODE3 not supported.");
            QUIT_PLUTO(1);

#endif  // FCB_DEBOOST_MODE

        } // Sufficient power for heat to ensure pressure equilibrium, i.e., need to deboost?

#endif  // if FCB_DEBOOST

    } // Is AGN on or off (sufficient power to ensure pressure equilibrium ?)


#else  // if not FEEDBACK_CYCLE

    power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
    speed = Lorentz2Speed(lorentz);
    gmm1 = g_gamma / (g_gamma - 1.);
    pres = power / (gmm1 * lorentz * lorentz * speed * nz.area * (1. + (lorentz - 1.) / lorentz * chi));
    dens = chi * gmm1 * pres;
    is_on = 1;

#endif  // whether FEEDBACK_CYCLE

    ofs->pow = power;
    ofs->mdt = chi;
    // NOTE: DANGER (when testing feedback cycle), changed from lorentz to speed.
    //       Maybe this is OK, and only change chi to actual mdot. Check feedback cycle code.
    ofs->spd = speed;
    ofs->rho = dens;
    ofs->prs = pres;
    ofs->eth = gmm1 * pres * speed * nz.area * lorentz * lorentz;
    ofs->kin = (lorentz - 1.) * lorentz * speed * nz.area * dens * CONST_c * CONST_c / (vn.v_norm * vn.v_norm);
    ofs->is_on = is_on;

}


/* ************************************************ */
void SetUfoState(OutflowState *ofs) {
/*
 * Set outflow state.
 * Set the outflow primitives/parameters struct.
 * Calculate primitives from input parameters
 * power, angle, speed, mdot, and radius and width.
 *
 ************************************************** */

    double power, speed, mdot;
    double heat, dens, pres, kine;
    int is_on;

    /* Input parameters, normalized in code units*/

    speed = g_inputParam[PAR_OSPD] * ini_code[PAR_OSPD];
    mdot = g_inputParam[PAR_OMDT] * ini_code[PAR_OMDT];

    // TODO: Test accretion cycle
#if FEEDBACK_CYCLE == YES

    /* Wind coupling efficiency */
    double oeff = g_inputParam[PAR_OEFF] * ini_code[PAR_OEFF];

    /* AGN power from accretion rate */
    power = (1. - ac.eff) / (1. + ac.mld) * oeff * ac.eff * ac.accr_rate * CONST_c * CONST_c / vn.pot_norm;

    /* Mass outflow rate from accretion rate */
    mdot = ac.accr_rate * g_dt * (1. - ac.eff) / (1. + ac.mld) * ac.mld;

    /* Maximum speed (if all of the power was in kinetic form) = sqrt(2 * power / mdot). */
    double vmax = (2. * oeff * ac.eff / ac.mld) * CONST_c / vn.v_norm;

    /* Make the above speed a ceiling for the speed - this might make deboosting redundant */
    speed = MIN(speed, vmax);


    /* Special case for time = 0 */
    if (g_time == 0) {

        /* This is for initialization of reference outflow state struct contained in the ac struct,
         * through SetAccretionPhysics at the beginning of Init. It does not affect dynamics of first timestep. */
        if (ac.deboost == 1) {

            power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
            heat = power - 0.5 * mdot * speed * speed;
            pres = heat * (g_gamma - 1) / (g_gamma * nz.area * speed);
            dens = mdot / (nz.area * speed);
            is_on = 1;
        }

        /* At time = 0, AGN is off */
        else {

            ac.deboost = 0;
            dens = g_smallDensity;
            pres = g_smallPressure;
            speed = 0;
            mdot = 0;
            power = 0;
            is_on = 0;

        }

    }

    /* If sufficient power, AGN is on */
    else if (FeedbackTrigger(power, speed)){

        /* Set pressure, density, and speed, adjusting parameters such that
         * pressure is be at least ambient pressure, hprs. */

        /* If sufficient power ... */
        kine = 0.5 * mdot * speed * speed;
        heat = power - kine;
        pres = heat * (g_gamma - 1) / (g_gamma * ac.nzi.area * speed);
        dens = mdot / (ac.nzi.area * speed);
        ac.deboost = 1.0;
        is_on = 1;

#if FBC_DEBOOST == YES
        /* ...but if insufficient heat with default Ekin, reduce (deboost) mdot and speed and area */
        ram_pres_max = 2 * kine / (ac.nzi.area * speed);
        pres_max = MAX(ram_pres_max, pres);

        if (pres_max < hprs) {

            /* Ensure pressure equilibrium */
            pres = hprs;

#if FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_0

            /* Adjust only area. The deboost factor simply reduces area, and thus increases the density. */

            ac.deboost = (g_gamma - 1) / gamma * (power - 0.5 * mdot * speed * speed) / (pres * nz.area * speed);

            /* boost density */
            dens /= ac.deboost;

#elif FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_1

            /* Adjust only mdot. The deboost factor simply reduces mdot, and thus the density. */

            ac.deboost = ( power - pres * nz.area * speed * g_gamma / (g_gamma - 1) ) / ( 0.5 * mdot * speed * speed );

            /* deboost mdot and density */
            mdot *= ac.deboost;
            dens *= ac.deboost;


#elif FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_2

            /* Adjust only area and mdot. The deboost factor is chosen such that the density remains constant:
             *
             *   mdot_db = mdot * deboost
             *   area_db = area * deboost
             *
             * */

            ac.deboost = power / ( pres * nz.area * speed * g_gamma / (g_gamma - 1) + 0.5 * mdot * speed * speed );

            /* deboost mdot */
            mdot *= ac.deboost;

#elif FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_3

            /* Adjust speed, area, and mdot. The deboost factor is chosen such that the density remains constant:
             *
             *   mdot_db = mdot * deboost
             *   speed_db = speed * √deboost
             *   area_db = area * √deboost
             *
             * */

            /* deboost should be 0 < deboost < 1 */
            double a = g_gamma * nz.area * pres * speed;
            double b = 2 * mdot * power * speed * speed;
            ac.deboost = (-a + sqrt(a * a + b * (1 - 2 * g_gamma + g_gamma * g_gamma))) /
                         (mdot * speed * speed * (g_gamma - 1));

            /* Deboost speed and mdot */
            speed *= sqrt(ac.deboost);
            mdot *= ac.deboost;

            /* Note that if sufficient heat with default Ekin, the increased pres is just used.
             * Speed and mdot two are effectively upper limits of what can be achieved */

#endif  // FBC_DEBOOST_MODE

        } // Sufficient power for heat to ensure pressure equilibrium, i.e., need to deboost?


#endif // Whether FCB_DEBOOST is ON

    }

    /* Not enough power - AGN is off */
    else {

        ac.deboost = 0;
        dens = g_smallDensity;
        pres = g_smallPressure;
        speed = 0;
        mdot = 0;
        power = 0;
        is_on = 0;

    } // Is AGN on or off (sufficient power to ensure pressure equilibrium ?)




#else   // if not FEEDBACK_CYCLE

    power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
    heat = power - 0.5 * mdot * speed * speed;

    pres = heat * (g_gamma - 1) / (g_gamma * nz.area * speed);
    dens = mdot / (nz.area * speed);

    is_on = 1;

#endif  // if FEEDBACK_CYCLE


    /* Fill outflow struct */
    ofs->pow = power;
    ofs->mdt = mdot;
    ofs->rho = dens;
    ofs->prs = pres;
    ofs->spd = speed;
    ofs->eth = pres * g_gamma / (g_gamma - 1) * speed * nz.area;
    ofs->kin = 0.5 * dens * speed * speed * speed * nz.area;
    ofs->is_on = is_on;

}

/* ************************************************ */
int FeedbackTrigger(double power, double speed) {

/*
 * Trigger function for AGN Feedback in feedback cycle
 *
 ************************************************** */

    /* Hot phase pressure */
//    double hrho = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
//    double htmp = g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP];
//    double hprs = PresIdealEOS(hrho, htmp, MU_NORM);
//
//    /* Maximum ram pressure of AGN outflow (Ram pressure ~ 5 th. pressure, for gamma = 5/3) */
//    /* Note, accretion struct must be initialized beforehand */
//    double ram_pres_max, pres_max;
//    ram_pres_max = 2 * power / ( ac.nzi.area * speed);

//    return (ram_pres_max > hprs);

    return (power * vn.power_norm > 1.e42);

}


/* ************************************************ */
void OutflowPrimitives(double *out_primitives, const double x1, const double x2, const double x3) {
/*
 * Runs the relevant primitives function for nozzle flow.
 *
 * NOTE:
 *  - JetPrimitives and UfoPrimities only differ in their
 *    normalizations and how the paramters are set.
 *
 ************************************************** */


    /* Set primitives array */
    out_primitives[RHO] = os.rho;
    out_primitives[PRS] = os.prs;
    out_primitives[TRC] = 1.0;
#if CLOUDS
    out_primitives[TRC + 1] = 0.0;
#endif
    OutflowVelocity(out_primitives, os.spd, x1, x2, x3);

}

/* ************************************************ */
void OutflowVelocity(double *out_primitives, double speed,
                     const double x1, const double x2, const double x3) {
/*
 * Calculate outflow velocity vector inside out_primitives given
 * a cell location x1, x2, x3, and velocity magnitude speed.
 *
 * This function is used by UfoPrimitives and JetPrimitives
 * and should be generic for any outflow nozzle.
 *
 ************************************************** */

    /* Work in spherical. The velocity vectors are radial along the cone.
     * The cone apex is not necessarily the precession axis or the coordinate origin. */

    /* Grid point in cartesian coordinates */
    double cx1, cx2, cx3;
    cx1 = CART1(x1, x2, x3);
    cx2 = CART2(x1, x2, x3);
    cx3 = CART3(x1, x2, x3);

    /* Rotate grid to align to flow axis */
    double cx1p, cx2p, cx3p;
    RotateGrid2Nozzle(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);

    int mirror_side = 0;

    /* Cartesian velocities in unrotated and rotated frames */
    double cv1, cv2, cv3;
    double cv1p, cv2p, cv3p;

    if (nz.is_fan) {

#if INTERNAL_BOUNDARY == YES
        /* If we're in the counter-nozzle region use mirror symmetric value of cx1p
         NOTE, since we've rotated to flow-axis, the nozzle is cylindrically symmetric
         and we don't have to use rotational symmetry. */

        /* Can't use FLOWAXIS macro though because we are in transformed coords. */

        if (D_SELECT(cx1p, cx2p, cx3p) < 0) {
            mirror_side = 1;
            D_SELECT(cx1p, cx2p, cx3p) *= -1;
        }
#endif

        /* Move cone apex to 0,0,0 */
        /* Can't use FLOWAXIS macro though because we are in transformed coords. */
        D_SELECT(cx1p, cx2p, cx3p) -= nz.cone_apex;

        /* Spherical coordinates of the rotated cartesian coordinates */
        double sx1p, sx2p, sx3p;
        sx1p = CART2SPH1(cx1p, cx2p, cx3p);
        sx2p = CART2SPH2(cx1p, cx2p, cx3p);
        sx3p = CART2SPH3(cx1p, cx2p, cx3p);

        /* Cartesian velocities from spherical velocity */
        EXPAND(cv1p = VSPH2CART1(sx1p, sx2p, sx3p, speed, 0, 0);,
               cv2p = VSPH2CART2(sx1p, sx2p, sx3p, speed, 0, 0);,
               cv3p = VSPH2CART3(sx1p, sx2p, sx3p, speed, 0, 0););

        /* Move cone apex back */
        D_SELECT(cx1p, cx2p, cx3p) += nz.cone_apex;

#if INTERNAL_BOUNDARY == YES
        /* Create mirror-symmetrically opposite velocity vector
         and mirror back cell position */

        if (mirror_side) {
            D_SELECT(cx1p, cx2p, cx3p) *= -1;
            D_SELECT(cv1p, cv2p, cv3p) *= -1;
        }

#endif

    }
    else { // if nozzle is not fan

        cv1p = 0;
        cv2p = 0;
        cv3p = 0;
        D_SELECT(cv1p, cv2p, cv3p) = speed;

#if INTERNAL_BOUNDARY == YES
        /* Create mirror-symmetrically opposite velocity vector */

        if (D_SELECT(cx1p, cx2p, cx3p) < 0) {
            D_SELECT(cv1p, cv2p, cv3p) *= -1;
        }
#endif

    } // whether nozzle is fan

    /* Rotate vector back */
    RotateNozzle2Grid(cv1p, cv2p, cv3p, &cv1, &cv2, &cv3);

    double vx1, vx2, vx3;
    EXPAND(vx1 = VCART_1(cx1, cx2, cx3, cv1, cv2, cv3);,
           vx2 = VCART_2(cx1, cx2, cx3, cv1, cv2, cv3);,
           vx3 = VCART_3(cx1, cx2, cx3, cv1, cv2, cv3););

    EXPAND(out_primitives[VX1] = vx1;,
           out_primitives[VX2] = vx2;,
           out_primitives[VX3] = vx3;);

    return;
}


/*********************************************************************** */
void OutflowOutput() {
/*!
 * Perform output of spherical accretion calculations
 *
 *********************************************************************** */

    if (prank == 0) {

        FILE *fp_out;
        char fname[512];
        sprintf(fname, "outflow.dat");
        static double next_output = -1;

        next_output = OutputContextEnter(fname, &fp_out, next_output, OUTFLOW_OUTPUT_RATE);

        /* Write data */
        if (g_time > next_output) {


            fprintf(fp_out, "%12.6e %2d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e \n",
                    g_time * vn.t_norm / (CONST_ly / CONST_c),            // time (yr)
                    os.is_on,                                             // Whether outflow is on
                    os.pow * vn.power_norm,                               // Outflow power
                    os.kin * vn.power_norm,                               // Outflow kinetic power
                    os.eth * vn.power_norm,                               // Outflow enthalpy (inj. rate)
                    os.rho,                                               // Outflow density (1/cm3)
                    os.prs * vn.pres_norm / CONST_kB,                     // Outflow pressure (p/k)
                    os.spd * vn.v_norm / 1.e5                             // Outflow speed (km/s)
            );

        }

        next_output = OutputContextExit(&fp_out, next_output, OUTFLOW_OUTPUT_RATE);

    }

}


