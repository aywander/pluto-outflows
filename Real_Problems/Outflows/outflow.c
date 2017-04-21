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

/* Global struct for nozzle */
Nozzle nz;

/* Global struct for outflow parameters */
OutflowState os;


/* ************************************************ */
void SetNozzleGeometry(Nozzle * noz) {
/*
 * Set outflow geometry struct with parameters of cone
 *
 *   Nozzle is always a fan shaped region defined by
 * the intersection of a cone with the boundary volume.
 * A cap is included at the top of the cone, making it
 * an ice-cream.
 *   The outer rim of the maximum cone radius touches at the
 * surface of the computational domain. The cone apex is
 * not necessarily at (0,0,0), neither is the tilt rotation
 * axis. The apex and tilt rotation axis are also not
 * necessarily at the same point.
 *    If the half opening angle
 * ang == 0, the ice-cream becomes a bullet.
 *
 * NOTE:
 *  - JetPrimitives and UfoPrimities only differ in their
 *    normalizations and how the paramters are set.
 *
 *
 ************************************************** */

    double ang, rad, sph, dir, dbh, omg, phi;
    double cbh, orig, area, vol, cone_height, cone_apex;
    int is_fan, is_two_sided;

    double small_angle = 1.e-12;
    double large = 1.e30;

    /* First determine if we are dealing with a cone or a parallel nozzle */
    ang = g_inputParam[PAR_OANG] * ini_code[PAR_OANG];
    if (ang > small_angle) is_fan = 1;
    else is_fan = 0;

    /* Determine if nozzle is one-sided */
#if GEOMETRY == SPHERICAL
    if (g_domEnd[JDIR] > 3. * CONST_PI / 4.) is_two_sided = 1;
#else
    if (FLOWAXIS(g_domBeg[IDIR], g_domBeg[JDIR], g_domBeg[KDIR]) < 0.) is_two_sided = 1;
#endif
    else is_two_sided = 0;


#if ACCRETION == YES && FEEDBACK_CYCLE == YES

    /* Accretion struct must be initialized beforehand */

    /* The angle and radius (ORAD) as input parameter defines the maximum angle */
    if (g_time == 0 && ac.deboost == 1) {

        /* Quantities from input parameters */
        rad = g_inputParam[PAR_ORAD] * ini_code[PAR_ORAD];
    }
    else {

        /* Reduction of radius */
        rad = sin(ang) * sqrt( ac.nzi.area * ac.deboost / (2. * CONST_PI * (1. - cos(ang))) ) ;
    }

#else

    /* Quantities from input parameters */
    rad = g_inputParam[PAR_ORAD] * ini_code[PAR_ORAD];

#endif   // whether ACCRETION and FEEDBACK_CYCLE

    /* Quantities from input parameters */
    ang = g_inputParam[PAR_OANG] * ini_code[PAR_OANG];
    sph = g_inputParam[PAR_OSPH] * ini_code[PAR_OSPH];
    dir = g_inputParam[PAR_ODIR] * ini_code[PAR_ODIR];
    dbh = g_inputParam[PAR_ODBH] * ini_code[PAR_ODBH];
    omg = g_inputParam[PAR_OOMG] * ini_code[PAR_OOMG];
    phi = g_inputParam[PAR_OPHI] * ini_code[PAR_OPHI];

    /* Derived quantitites */

#if INTERNAL_BOUNDARY == YES
    cbh = sqrt(sph * sph - rad * rad);
    orig = cbh;
#else
    cbh = (orig - dbh) / cos(dir) + rad * tan(dir);
    orig = g_domBeg[FLOWAXIS(IDIR, JDIR, KDIR)];
#endif

    /* Equations are in terms of rad, rather than sph,
     * so that they are valid also if there are not internal boundaries. */

    if (is_fan) {

        /* cone apex is only valid when cone is aligned with flow axis */
        cone_height = rad / tan(ang);
        cone_apex = orig - cone_height;

        /* Area of nozzle */
        double cone_side = rad / sin(ang);
        area = 2. * CONST_PI * (1. - cos(ang)) * cone_side * cone_side;

        /* Volume of nozzle. It is the revolved volume. */
        // TODO: The volume is only valid in the case of two-sided galaxies, where dbh = 0.
#if INTERNAL_BOUNDARY == YES
        double spherical_cone_volume = area * cone_side / 3.;
        double overlap_rad = rad * (1. - cbh / cone_height);
        double overlap_volume = CONST_PI * overlap_rad * overlap_rad * (cone_height - cbh) / 3.;
        vol = spherical_cone_volume - overlap_volume ;
#else
        vol = 0;
#endif

    }
    else {
        cone_height = large;
        cone_apex = -large;

        /* Area and volume of nozzle */
        area = CONST_PI * rad * rad;

        /* Volume of nozzle */
        // TODO: The volume is only valid in the case of two-sided galaxies, where dbh = 0.
#if INTERNAL_BOUNDARY == YES
        vol = area * cbh;
#else
        vol = 0;
#endif

    }

    /* Power should always be the total power injected into box */
    // TODO: Make sure two-sided and one-sided cases are consistent.
    if (is_two_sided){
        area *= 2.;
        vol *= 2.;
    }

    /* Some consistency checks */

#if INTERNAL_BOUNDARY == YES
    /* Currently we require the BH location to be at 0,0,0  and noz.dbh = 0*/

    if (0 > dbh || dbh > 0) {
        print1("Warning: If INTERNAL_BOUNDARY == YES, ODBH = 0. Setting noz.dbh to 0.");
        dbh = 0;
    }
    if (sph < rad) {
        print1("Error: OSPH must be larger than ORAD");
        QUIT_PLUTO(1);
    }
    /* This check is for avoiding the cone of the nozzle to be buried in the
     outflow boundary for half-galaxy simulations. */
    if (FLOWAXIS(g_domBeg[IDIR], g_domBeg[JDIR], g_domBeg[KDIR]) >= 0.) {

        if (acos(1 - ((sph - cbh) * (sph - cbh) + rad * rad) /
                     (2 * sph * sph)) + dir > CONST_PI / 2) {
            print1("Error: OSPH is too small. It must be at least...");
            QUIT_PLUTO(1);
            /* TODO: Calculate minimum OSPH */
        }
    }

#else // INTERNAL_BOUNDARY

#if NOZZLE_FILL == NF_CONSERVATIVE

    print1("Error: INTERNAL_BOUNDARY must be TRUE if NOZZLE_FILL == NF_CONSERVATIVE.");
    QUIT_PLUTO(1);

#endif

#endif

    noz->ang = ang;
    noz->rad = rad;
    noz->sph = sph;
    noz->dir = dir;
    noz->dbh = dbh;
    noz->omg = omg;
    noz->phi = phi;
    noz->cbh = cbh;
    noz->orig = orig;
    noz->area = area;
    noz->vol = vol;
    noz->cone_height = cone_height;
    noz->cone_apex = cone_apex;
    noz->is_fan = is_fan;
    noz->is_two_sided = is_two_sided;

}


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

    /* The parameters from ini, and normalize them */
    chi = g_inputParam[PAR_OMDT] * ini_code[PAR_OMDT];
    lorentz = g_inputParam[PAR_OSPD] * ini_code[PAR_OSPD];


    // TODO: Test accretion cycle
#if ACCRETION == YES && FEEDBACK_CYCLE == YES

    /* AGN power from accretion rate */
    power = ac.eff * ac.accr_rate * CONST_c * CONST_c / vn.pot_norm;

    /* Hot phase pressure */
    double hrho = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
    double htmp = g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP];
    double hprs = PresIdealEOS(hrho, htmp, MU_NORM);

    /* Maximum thermal pressure of AGN outflow */
    speed = Lorentz2Speed(lorentz);
    gmm1 = g_gamma / (g_gamma - 1.);
    pres = power / (gmm1 * lorentz * lorentz * speed * ac.nzi.area);

    /* This is for initialization of reference outflow state struct
     * in SetAccretionPhysics, contained within ac struct */
    if (g_time == 0 && ac.deboost == 1) {
        power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
        speed = Lorentz2Speed(lorentz);
        gmm1 = g_gamma / (g_gamma - 1.);
        pres = power / (gmm1 * lorentz * lorentz * speed * nz.area * (1. + (lorentz - 1.) / lorentz * chi));
        dens = chi * gmm1 * pres;
    }

    /* If insufficient power for pressure equilibrium, shut off AGN */
    else if ((pres < hprs) || (g_time == 0 )) {
        ac.deboost = 0;
        dens = g_smallDensity;
        pres = g_smallPressure;
        lorentz = 0;
        chi = 0;
        power = 0;
    }

        /* If sufficient power, AGN is on */
    else {

        /* Set pressure, density, and speed, adjusting parameters such that
         * pressure is be at least ambient pressure, hprs. */

        /* If sufficient power ... */
        pres = power / (gmm1 * lorentz * lorentz * speed * ac.nzi.area * (1. + (lorentz - 1.) / lorentz * chi));
        dens = chi * gmm1 * pres;

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

            print1("Error in SetJetState for Feedback cycle FBC_DEBOOST_MODE3 not supported.");
            QUIT_PLUTO(1);

#endif  // FCB_DEBOOST_MODE

        } // Sufficient power for heat to ensure pressure equilibrium, i.e., need to deboost?

    } // Is AGN on or off (sufficient power to ensure pressure equilibrium ?)

#else  // if not ACCRETION && FEEDBACK_CYCLE

    power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
    speed = Lorentz2Speed(lorentz);
    gmm1 = g_gamma / (g_gamma - 1.);
    pres = power / (gmm1 * lorentz * lorentz * speed * nz.area * (1. + (lorentz - 1.) / lorentz * chi));
    dens = chi * gmm1 * pres;

#endif  // whether ACCRETION && FEEDBACK_CYCLE

    ofs->pow = power;
    ofs->mdt = chi;
    // NOTE: DANGER (when testing feedback cycle), changed from lorentz to speed.
    //       Maybe this is OK, and only change chi to actual mdot. Check feedback cycle code.
    ofs->spd = speed;
    ofs->rho = dens;
    ofs->prs = pres;
    ofs->eth = gmm1 * pres * speed * nz.area * lorentz * lorentz;
    ofs->kin = (lorentz - 1.) * lorentz * speed * nz.area * dens * CONST_c * CONST_c / (vn.v_norm * vn.v_norm);

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
    double heat, dens, pres;

    /* Input parameters, normalized in code units*/

    speed = g_inputParam[PAR_OSPD] * ini_code[PAR_OSPD];
    mdot = g_inputParam[PAR_OMDT] * ini_code[PAR_OMDT];

    // TODO: Test accretion cycle
#if ACCRETION == YES && FEEDBACK_CYCLE == YES

    /* AGN power from accretion rate */
    power = ac.eff * ac.accr_rate * CONST_c * CONST_c / vn.pot_norm;

    /* Hot phase pressure */
    double hrho = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
    double htmp = g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP];
    double hprs = PresIdealEOS(hrho, htmp, MU_NORM);

    /* Thermal pressure of AGN outflow */
    /* Note, accretion struct must be initialized beforehand */
    pres = power * (g_gamma - 1) / (g_gamma * ac.nzi.area * speed);

    /* This is for initialization of reference outflow state struct
     * in SetAccretionPhysics, contained within ac struct */
    if (g_time == 0 && ac.deboost == 1) {
        power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
        heat = power - 0.5 * mdot * speed * speed;
        pres = heat * (g_gamma - 1) / (g_gamma * nz.area * speed);
        dens = mdot / (nz.area * speed);

    }

    /* If insufficient power for pressure equilibrium, shut off AGN */
    else if ((pres < hprs) || (g_time == 0)) {
        ac.deboost = 0;
        dens = g_smallDensity;
        pres = g_smallPressure;
        speed = 0;
        mdot = 0;
        power = 0;
    }

    /* If sufficient power, AGN is on */
    else {

        /* Set pressure, density, and speed, adjusting parameters such that
         * pressure is be at least ambient pressure, hprs. */

        /* If sufficient power ... */
        heat = power - 0.5 * mdot * speed * speed;
        pres = heat * (g_gamma - 1) / (g_gamma * ac.nzi.area * speed);
        dens = mdot / (ac.nzi.area * speed);

        /* ...but if insufficient heat with default Ekin, reduce (deboost) mdot and speed and area */
        if (pres < hprs) {

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

    } // Is AGN on or off (sufficient power to ensure pressure equilibrium ?)



#else   // if not ACCRETION and FEEDBACK_CYCLE

    power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
    heat = power - 0.5 * mdot * speed * speed;

    pres = heat * (g_gamma - 1) / (g_gamma * nz.area * speed);
    dens = mdot / (nz.area * speed);

#endif  // if ACCRETION and FEEDBACK_CYCLE


    /* Fill outflow struct */
    ofs->pow = power;
    ofs->mdt = mdot;
    ofs->rho = dens;
    ofs->prs = pres;
    ofs->spd = speed;
    ofs->eth = pres * g_gamma / (g_gamma - 1);
    ofs->kin = 0.5 * dens * speed * speed;

}


/* ************************************************ */
void NozzleFill(const Data *d, const Grid *grid) {
/*
 * Dumps mass, momentum, and energy into nozzle region,
 * according to nozzle parameters.
 *
 ************************************************** */

    RBox *box = GetRBox(DOM, CENTER);

    // Just for initialization, need to convert primitives to conservatives
    // TODO: This may not be needed if this is modularized and attached somewhere else in main.c
    if (g_time <= 0.) PrimToCons3D(d->Vc, d->Uc, box);

    double *x1, *x2, *x3;
    double out_conservatives[NVAR];
    int nv;

    /* zero buffer */
    VAR_LOOP(nv) out_conservatives[nv] = 0.;

    /* These are the geometrical central points */
    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;


    // TODO: Relativistic version
#if (PHYSICS == RHD) || (PHYSICS == RMHD)
    print1("Error: NOZZLE_FILL CONSERVATIVE not supported if PHYSICS is RHD or RMHD.");
    QUIT_PLUTO(1);
#endif

    // TODO: Adjust volume for the case below
#if GEOMETRY == SPHERICAL
    if (g_domBeg[IDIR] > 0.) {
        print1("Error: NOZZLE_FILL CONSERVATIVE not supported if (GEOMETRY == SPHERICAL) && (g_domBeg[IDIR] > 0.).");
        QUIT_PLUTO(1);
    }
#endif

    // TODO: Need both old and new timesteps here. Need to do this from main.c
    double dt_nf = g_dt;

    /* Total energy to dump per unit volume */
    double energy_dump = os.pow * dt_nf / nz.vol;

    /* Total mass to dump per unit volume */
    double mass_dump = os.mdt * dt_nf / nz.vol;

    /* Momentum input per unit volume */
    double momentum_dump = mass_dump * os.spd;

    /* Weights */
    // Do uniform dump first.

    // TODO: Conservative variables are specific energy ?
    // TODO: Volume is incorrectly calculated for spherical setups; nothing is dumped in the excluded x1_beg region.

    /* Apply on all cells in Nozzle region */
    int k, j, i;
    DOM_LOOP(k, j, i) {

                if (InNozzleRegion(x1[i], x2[j], x3[k])) {

                    out_conservatives[RHO] = mass_dump;
                    out_conservatives[ENG] = energy_dump;
                    OutflowVelocity(out_conservatives, momentum_dump, x1[i], x2[j], x3[k]);

                    VAR_LOOP(nv) d->Uc[k][j][i][nv] += out_conservatives[nv];

                }

            }

    /* Update primitives */
    ConsToPrim3D(d->Uc, d->Vc, d->flag, box);
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

    return;

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
        cv1p = VSPH2CART1(sx1p, sx2p, sx3p, speed, 0, 0);
        cv2p = VSPH2CART2(sx1p, sx2p, sx3p, speed, 0, 0);
        cv3p = VSPH2CART3(sx1p, sx2p, sx3p, speed, 0, 0);

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

#if RECONSTRUCT_4VEL == YES
    /* This is the same for all geometries */
    double lorentz;
    lorentz = Speed2Lorentz(speed);
    EXPAND(out_primitives[VX1] *= lorentz;,
           out_primitives[VX2] *= lorentz;,
           out_primitives[VX3] *= lorentz;);
#endif

    return;
}


/* ************************************************ */
int InNozzleCap(const double x1, const double x2, const double x3) {
/*
 * Returns 1 if r is in nozzle cap, 0 if not.
 * Nozzle cap is a hemispherical region above the base of
 * the outflow cone.
 *
 ************************************************** */

#if NOZZLE_CAP == YES

    /* Nozzle is be deactivated if ORAD is zero */
    if (nz.rad == 0) return 0;

    /* Grid point in cartesian coordinates */
    double cx1, cx2, cx3;
    cx1 = CART1(x1, x2, x3);
    cx2 = CART2(x1, x2, x3);
    cx3 = CART3(x1, x2, x3);

    /* Rotate and shift cartesian coords so that base of cone is at (0,0,0)
       This is the same regardless if we are using internal boundary or not.
       But for internal boundary, we're including a mirrored component.
       Inclusion of the mirror component with fabs must occur after rotation
       but before translation. */
    double cx1p, cx2p, cx3p;
    RotateGrid2Nozzle(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);
#if INTERNAL_BOUNDARY
    D_SELECT(cx1p, cx2p, cx3p) = fabs(D_SELECT(cx1p, cx2p, cx3p));
#endif
    D_SELECT(cx1p, cx2p, cx3p) -= nz.dbh + nz.cbh;

    /* Turn into spherical coords */
    double sr, st;
    sr = CART2SPH1(cx1p, cx2p, cx3p);
    st = CART2SPH2(cx1p, cx2p, cx3p);

    if (sr < nz.rad){
       int a = 1;
    }
    /* Do radius and angle check */
    return (sr < nz.rad) && (st < CONST_PI / 2.);

#else

    /* if NOZZLE_CAP == NO */
    return 0;

#endif

}

/* ************************************************ */
int InNozzleRegion(const double x1, const double x2, const double x3) {
/*
 * Returns 1 if r is in outflow region, 0 if not.
 * The outflow region is defined for a fan-shaped (conical)
 * outlet as as the conical region capped with a
 * spherical section with radius of the cone edge.
 *
 ************************************************** */


    /* Nozzle is be deactivated if ORAD is zero */
    if (nz.rad == 0.) return 0;

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

    if (nz.is_fan) {

        /* Shift so that cone apex is at (0,0,0) */
        D_SELECT(cx1p, cx2p, cx3p) -= nz.cone_apex;

        /* Turn into spherical coords */
        double sr, st;
        sr = CART2SPH1(cx1p, cx2p, cx3p);
        st = CART2SPH2(cx1p, cx2p, cx3p);

        return (sr < nz.cone_height / cos(nz.ang)) && (st < nz.ang);
    }

    else {

        /* Shift cartesian coords so that base of cone is at (0,0,0) */
        D_SELECT(cx1p, cx2p, cx3p) -= nz.dbh + nz.cbh;

        /* Turn into cylindrical coords */
        double cr, cz;
        cr = CART2CYL1(cx1p, cx2p, cx3p);
        cz = CART2CYL2(cx1p, cx2p, cx3p);

        return (cr < nz.rad) && (cz < 0);
    }
}

/* ************************************************ */
double Profile(const double x1, const double x2, const double x3)
/*!
  * Some smoothing function, e.g., 1/cosh
  *
 ************************************************** */
{

    if (1) { return 1.0; }
        /* NOTE: Currently we don't use a smoothing profile.
         This won't work well with internal boundaries anyway. */

    else {

        /* Steepness of cosh profile */
        int n = 14;

        /* Grid point in cartesian coordinates */
        double cx1, cx2, cx3, cx1p, cx2p, cx3p;
        cx1 = CART1(x1, x2, x3);
        cx2 = CART2(x1, x2, x3);
        cx3 = CART3(x1, x2, x3);

        /* Rotate cartesian coords and turn into cylindrical coords */
        double cr, cz;
        RotateGrid2Nozzle(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);
        cr = CART2CYL1(cx1p, cx2p, cx3p);
        cz = CART2CYL2(cx1p, cx2p, cx3p);

        /* Return smoothing factor */
        return 1.0 / cosh(pow(cr / nz.rad, n));
    }
}

