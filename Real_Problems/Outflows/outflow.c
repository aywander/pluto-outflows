//
// Created by Alexander Y. Wagner on 4/7/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "outflow.h"
#include "grid_geometry.h"
#include "hot_halo.h"

/* Global struct for nozzle */
Nozzle nz;

/* ************************************************ */
void SetNozzleGeometry() {
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

    double small_angle = 1.e-12;
    double large = 1.e30;

    /* Quantities from input parameters */
    nz.ang = g_inputParam[PAR_OANG] * ini_code[PAR_OANG];
    nz.rad = g_inputParam[PAR_ORAD] * ini_code[PAR_ORAD];
    nz.dir = g_inputParam[PAR_ODIR] * ini_code[PAR_ODIR];
    nz.dbh = g_inputParam[PAR_ODBH] * ini_code[PAR_ODBH];
    nz.omg = g_inputParam[PAR_OOMG] * ini_code[PAR_OOMG];
    nz.phi = g_inputParam[PAR_OPHI] * ini_code[PAR_OPHI];
    nz.sph = g_inputParam[PAR_OSPH] * ini_code[PAR_OSPH];

    /* Derived quantitites */

#if INTERNAL_BOUNDARY == YES
    nz.cbh = sqrt(nz.sph * nz.sph - nz.rad * nz.rad);
    nz.orig = nz.cbh;
#else
    nz.cbh = (nz.orig - nz.dbh) / cos(nz.dir) + nz.rad * tan(nz.dir);
    nz.orig = g_domBeg[FLOWAXIS(IDIR, JDIR, KDIR)];
#endif

    if (nz.ang > small_angle) {
        nz.isfan = 1;
        nz.area = 2. * CONST_PI * (1. - cos(nz.ang)) * pow(nz.rad / sin(nz.ang), 2);

        /* cone apex is only valid when cone is aligned with flow axis */
        nz.cone_height = nz.rad / tan(nz.ang);
        nz.cone_apex = nz.orig - nz.cone_height;
    }
    else {
        nz.isfan = 0;
        nz.area = CONST_PI * nz.rad * nz.rad;

        nz.cone_height = large;
        nz.cone_apex = -large;
    }


    /* Some consistency checks */

#if INTERNAL_BOUNDARY == YES
    /* Currently we require the BH location to be at 0,0,0  and nz.dbh = 0*/

    if (0 > nz.dbh || nz.dbh > 0) {
        print1("Warning: If INTERNAL_BOUNDARY == YES, ODBH = 0. Setting nz.dbh to 0.");
        nz.dbh = 0;
    }
    if (nz.sph < nz.rad) {
        print1("Error: OSPH must be larger than ORAD");
        QUIT_PLUTO(1);
    }
    /* This check is for avoiding the cone of the nozzle to be buried in the
     outflow boundary for half-galaxy simulations. */
    if (FLOWAXIS(g_domBeg[IDIR], g_domBeg[JDIR], g_domBeg[KDIR]) >= 0.) {

        if (acos(1 - ((nz.sph - nz.cbh) * (nz.sph - nz.cbh) + nz.rad * nz.rad) /
                     (2 * nz.sph * nz.sph)) + nz.dir > CONST_PI / 2) {
            print1("Error: OSPH is too small. It must be at least...");
            QUIT_PLUTO(1);
            /* TODO: Calculate minimum OSPH */
        }
    }
#endif

}

/* ************************************************ */
void OutflowPrimitives(double *out_primitives, const double x1, const double x2, const double x3,
                       const double accr_rate) {
/*
 * Runs the relevant primitives function for nozzle flow.
 *
 * NOTE:
 *  - JetPrimitives and UfoPrimities only differ in their
 *    normalizations and how the paramters are set.
 *
 ************************************************** */

    /* Get primitives array depending on outflow type */
    NOZZLE_SELECT(JetPrimitives(out_primitives, x1, x2, x3, accr_rate),
                  UfoPrimitives(out_primitives, x1, x2, x3, accr_rate));

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

    if (nz.isfan) {

#if INTERNAL_BOUNDARY == YES
        /* If we're in the counter-nozzle region use mirror symmetric value of cx1p
         NOTE, since we've rotated to flow-axis, the nozzle is cylindrically symmetric
         and we don't have to use rotational symmetry. */

        if (FLOWAXIS(cx1p, cx2p, cx3p) < 0) {
            mirror_side = 1;
            FLOWAXIS(cx1p, cx2p, cx3p) *= -1;
        }
#endif

        /* Move cone apex to 0,0,0 */
        FLOWAXIS(cx1p, cx2p, cx3p) -= nz.cone_apex;

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
        FLOWAXIS(cx1p, cx2p, cx3p) += nz.cone_apex;

#if INTERNAL_BOUNDARY == YES
        /* Create mirror-symmetrically opposite velocity vector
         and mirror back cell position */

        if (mirror_side) {
            FLOWAXIS(cx1p, cx2p, cx3p) *= -1;
            FLOWAXIS(cv1p, cv2p, cv3p) *= -1;
        }

#endif

    }
    else { // if nozzle is not fan

        cv1p = 0;
        cv2p = 0;
        cv3p = 0;
        FLOWAXIS(cv1p, cv2p, cv3p) = speed;

#if INTERNAL_BOUNDARY == YES
        /* Create mirror-symmetrically opposite velocity vector */

        if (FLOWAXIS(cx1p, cx2p, cx3p) < 0) {
            FLOWAXIS(cv1p, cv2p, cv3p) *= -1;
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

#if USE_FOUR_VELOCITY == YES
    /* This is the same for all geometries */
    double lorentz;
    lorentz = Vel2Lorentz(speed);
    EXPAND(out_primitives[VX1] *= lorentz;,
           out_primitives[VX2] *= lorentz;,
           out_primitives[VX3] *= lorentz;);
#endif

    return;
}

/* ************************************************ */
void JetPrimitives(double *jet_primitives, const double x1, const double x2, const double x3, const double accr_rate) {
/*
 * Returns the array of primitives for jet parameters
 * power, chi, and lorentz, where
 * chi = (gmm-1)rho_j c^2/gmm p_j
 * The power is in erg / s
 * and the radius is in kpc
 *
 * NOTE:
 *  - combine this with UfoPrimities into
 *    OutflowPrimitives. Might not need NOZZLE_SELECT
 *    anymore.
 ************************************************** */

    double power, chi, lorentz;

    /* The parameters from ini, and normalize them */
    power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
    chi = g_inputParam[PAR_OMDT] * ini_code[PAR_OMDT];
    lorentz = g_inputParam[PAR_OSPD] * ini_code[PAR_OSPD];

    /* Some derived quantities */
    double vel, dens, pres, gmm1;
    vel = Lorentz2Vel(lorentz);
    gmm1 = g_gamma / (g_gamma - 1.);
    pres = power / (gmm1 * lorentz * lorentz * vel * nz.area *
                    (1. + (lorentz - 1.) / lorentz * chi));
    dens = chi * gmm1 * pres;

    /* Set primitives array */
    jet_primitives[RHO] = dens;
    jet_primitives[PRS] = pres;
    jet_primitives[TRC] = 1.0;
#if CLOUDS
    jet_primitives[TRC + 1] = 0.0;
#endif
    OutflowVelocity(jet_primitives, vel, x1, x2, x3);

    return;
}

/* ************************************************ */
void UfoPrimitives(double *ufo_primitives, const double x1, const double x2, const double x3, const double accr_rate) {
/*
 * Returns the array of primitives for ufo parameters
 * power, angle, speed, mdot, and radius and width
 *
 ************************************************** */

    double power, speed, mdot;
    double heat, dens, pres;

    /* Input parameters, normalized in code units*/

#if ACCRETION
//    power = ac.eff * accr_rate * CONST_c * CONST_c * vn.pot_norm;
    power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
#else
    power = g_inputParam[PAR_OPOW] * ini_code[PAR_OPOW];
#endif
    speed = g_inputParam[PAR_OSPD] * ini_code[PAR_OSPD];
    mdot = g_inputParam[PAR_OMDT] * ini_code[PAR_OMDT];


    /* Derived quantities */

    heat = power - 0.5 * mdot * speed * speed > 0;
    if (heat > 0) {
        pres = (power - 0.5 * mdot * speed * speed) * (g_gamma - 1) / (g_gamma * nz.area * speed);
        dens = mdot / (nz.area * speed);
        OutflowVelocity(ufo_primitives, speed, x1, x2, x3);
    }
    else{
        double halo[NVAR];
        HotHaloPrimitives(halo, x1, x2, x3);
//        dens = halo[RHO];
//        pres = halo[PRS];
        pres = g_smallPressure;
        dens = g_smallDensity;
        EXPAND(ufo_primitives[VX1] = 0;,
               ufo_primitives[VX2] = 0;,
               ufo_primitives[VX3] = 0;);
    }

    /* Set primitives array */
    ufo_primitives[RHO] = dens;
    ufo_primitives[PRS] = pres;
    ufo_primitives[TRC] = 1.0;
#if CLOUDS
    ufo_primitives[TRC + 1] = 0.0;
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

    if (nz.isfan) {

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