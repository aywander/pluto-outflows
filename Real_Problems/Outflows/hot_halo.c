//
// Created by Alexander Y. Wagner on 4/7/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"
#include "hot_halo.h"
#include "read_grav_table.h"
#include "read_hot_table.h"
#include "abundances.h"
#include "idealEOS.h"
#include "outflow.h"
#include "accretion.h"
#include "grid_geometry.h"


/* ************************************************ */
void HaloOuterBoundary(const int side, const Data *d, int i, int j, int k, const double x1, const double x2,
                       const double x3, int *touch) {/* Get primitives array for hot halo*/
/*
 * Set the values for the outer boundary conditions on sides that do not contain the jet.
 *
 * This routine is called from UserDefBoundary.
 *
 ************************************************** */


    int nv;
    double vmag, vmag_small = 1.e-6;
    double halo_primitives[NVAR];
    double vx1, vx2, vx3;

    HotHaloPrimitives(halo_primitives, x1, x2, x3);

    /* fill array */

    /* Calculate velocity */
    switch(side) {
        case X1_BEG:
            EXPAND(vx1 = d->Vc[VX1][k][j][IBEG];,
                   vx2 = d->Vc[VX2][k][j][IBEG];,
                   vx3 = d->Vc[VX3][k][j][IBEG];);
            break;
        case X1_END:
            EXPAND(vx1 = d->Vc[VX1][k][j][IEND];,
                   vx2 = d->Vc[VX2][k][j][IEND];,
                   vx3 = d->Vc[VX3][k][j][IEND];);
            break;
        case X2_BEG:
            EXPAND(vx1 = d->Vc[VX1][k][JBEG][i];,
                   vx2 = d->Vc[VX2][k][JBEG][i];,
                   vx3 = d->Vc[VX3][k][JBEG][i];);
            break;
        case X2_END:
            EXPAND(vx1 = d->Vc[VX1][k][JEND][i];,
                   vx2 = d->Vc[VX2][k][JEND][i];,
                   vx3 = d->Vc[VX3][k][JEND][i];);
            break;
        case X3_BEG:
            EXPAND(vx1 = d->Vc[VX1][KBEG][j][i];,
                   vx2 = d->Vc[VX2][KBEG][j][i];,
                   vx3 = d->Vc[VX3][KBEG][j][i];);
            break;
        case X3_END:
            EXPAND(vx1 = d->Vc[VX1][KEND][j][i];,
                   vx2 = d->Vc[VX2][KEND][j][i];,
                   vx3 = d->Vc[VX3][KEND][j][i];);
            break;
        default:
            QUIT_PLUTO(1);
    }

    vmag = VMAG(x1, x2, x3, vx1, vx2, vx3);

    /* Set outflow (zero grad) if cell outside boundary has some velocity,
     * else set to Hot halo. */

    if (vmag > vmag_small || *touch == 1) {

        switch(side) {
            case X1_BEG:
                for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
                break;
            case X1_END:
                for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
                break;
            case X2_BEG:
                for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i] = d->Vc[nv][k][JBEG][i];
                break;
            case X2_END:
                for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i] = d->Vc[nv][k][JEND][i];
                break;
            case X3_BEG:
                for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i] = d->Vc[nv][KBEG][j][i];
                break;
            case X3_END:
                for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i] = d->Vc[nv][KEND][j][i];
                break;
            default:
                QUIT_PLUTO(1);
        }

        *touch = 1;
    }

    else
        for (nv = 0; nv < NVAR; ++nv) d->Vc[nv][k][j][i] = halo_primitives[nv];

}

/* ************************************************************** */
void HotHaloPrimitives(double *halo,
                       const double x1, const double x2, const double x3) {
/*
 * Return array of primitives containing Halo quantities
 *
 * double    halo         array of halo primitives
 * double    x1, x2, x3   first, second, third coordinate
 *
 **************************************************************** */

    double vel, scrh;
    double inv_unit_G;
    double r, a, rs, rho0;
    double frac;
    double y0, y1, y2, y3, r1, r2;
    int iv, il, ic;
    static int once01 = 0;


    /* Initialize gravity arrays - as good as any other place to do it */
#ifdef GRAV_TABLE
    if (gr_rad == NULL) {
        ReadGravTable();
    }
    if (hot_rad == NULL) {
        ReadHotTable();
    }
#endif

    /* Consider different distributions */

    /* Hernquist potential (hydrostatic)*/
#if GRAV_POTENTIAL == GRAV_HERNQUIST
    rho0 = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
    r = SPH1(x1, x2, x3);
    a = g_inputParam[PAR_HRAD] * ini_code[PAR_HRAD];
    rs = r / a;
    halo[RHO] = rho0 / (rs * pow((1 + rs), 3));
    halo[PRS] = -2 * CONST_PI * CONST_G / vn.newton_norm * rho0 * rho0 * a * a *
                (pow(1. + rs, -4) / 12. * (25. + 52. * rs + 42. * rs * rs + 12. * rs * rs * rs) -
                log(1. + rs) + log(rs));


    /* This is also the default */
#elif (GRAV_POTENTIAL == GRAV_HOMOGENEOUS) || (GRAV_POTENTIAL == NONE)

    halo[RHO] = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
    halo[PRS] = PresIdealEOS(halo[RHO], g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP], MU_NORM);


#elif defined(GRAV_TABLE)

    /* The density from the table is in units of cm^-3 */
    r = SPH1(x1, x2, x3);
    halo[RHO] = InterpolationWrapper(hot_rad, hot_rho, hot_ndata, r);
    halo[RHO] = halo[RHO] * g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];

#if (GRAV_POTENTIAL == GRAV_SINGLE_ISOTHERMAL) || (GRAV_POTENTIAL == GRAV_DOUBLE_ISOTHERMAL)

    halo[PRS] = PresIdealEOS(halo[RHO], g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP], MU_NORM);

#else

    halo[PRS] = InterpolationWrapper(hot_rad, hot_prs, hot_ndata, r);

#endif // if ISOTHERMAL

#endif // GRAV_POTENTIAL types


    /* Velocities. */

    double vx1, vx2, vx3;

    if (g_inputParam[PAR_HVRD] != 0) {

        double vrad;
        double sx1, sx2, sx3;

        vrad = g_inputParam[PAR_HVRD] * ini_code[PAR_HVRD];

        /* Convert current coordinates to spherical */
        sx1 = SPH1(x1, x2, x3);
        sx2 = SPH2(x1, x2, x3);
        sx3 = SPH3(x1, x2, x3);

        /* Convert spherical velocity to back to current coordinates */
        EXPAND(vx1 = VSPH_1(sx1, sx2, sx3, vrad, 0, 0);,
               vx2 = VSPH_2(sx1, sx2, sx3, vrad, 0, 0);,
               vx3 = VSPH_3(sx1, sx2, sx3, vrad, 0, 0););

    }
    else {

        double vcart1, vcart2, vcart3;
        double xcart1, xcart2, xcart3;

        /* Convert coordinates and velocity vectors to cartesian */

        xcart1 = CART1(x1, x2, x3);
        xcart2 = CART2(x1, x2, x3);
        xcart3 = CART3(x1, x2, x3);

        /* Apply constant velocity (assumed to be in km/s) */
        vcart1 = g_inputParam[PAR_HVX1] * ini_code[PAR_HVX1];
        vcart2 = g_inputParam[PAR_HVX2] * ini_code[PAR_HVX2];
        vcart3 = g_inputParam[PAR_HVX3] * ini_code[PAR_HVX3];


        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(vx1 = VCART_1(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3);,
               vx2 = VCART_2(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3);,
               vx3 = VCART_3(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3););

    }

    EXPAND(halo[VX1] = vx1;, halo[VX2] = vx2;, halo[VX3] = vx3;);

#if USE_FOUR_VELOCITY == YES
    vel = VMAG(x1, x2, x3, halo[VX1], halo[VX2], halo[VX3]);
    scrh = Vel2Lorentz(vel);
    EXPAND(halo[VX1] *= scrh;, halo[VX2] *= scrh;, halo[VX3] *= scrh;);
#endif


    /* Tracers */
    halo[TRC] = 0.0;
#if CLOUDS
    halo[TRC + 1] = 0.0;
#endif


    /* A message for having initialized halo with potential*/
    if (!once01) {
        print1("> Initializing hot halo distribution of type: %d\n\n", GRAV_POTENTIAL);
        once01 = 1;
    }

    return;

}

/* ************************************************ */
int InFlankRegion(const double x1, const double x2, const double x3) {
/*
 * Returns 1 if r is in nearly spherical region excluding nozzle
 * and excluding the region above the nozzle.
 * It also excludes accretion region if accretion is on.
 * Only for INTERNAL_BOUNDARY YES
 * Sphere is assumed to be at (0,0,0) and has a radius
 * of OSPH.
 ************************************************** */

    /* Nozzle is be deactivated if OSPH is zero */
    if (nz.sph == 0.) return 0;

    /* Do radius check. Get spherical r coordinate */
    double sr = SPH1(x1, x2, x3);
    if (sr > nz.sph) return 0;
#if ACCRETION == YES
    if (sr < ac.snk) return 0;
#endif

    /* Grid point in cartesian coordinates */
    double cx1, cx2, cx3;
    cx1 = CART1(x1, x2, x3);
    cx2 = CART2(x1, x2, x3);
    cx3 = CART3(x1, x2, x3);

    /* Rotate cartesian coordinates */
    double cx1p, cx2p, cx3p;
    RotateGrid2Nozzle(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);

    /* we include a mirrored component with fabs.
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