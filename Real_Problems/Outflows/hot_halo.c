//
// Created by Alexander Y. Wagner on 4/7/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"
#include "hot_halo.h"
#include "read_grav_table.h"
#include "read_hot_table.h"
#include "idealEOS.h"
#include "outflow.h"
#include "accretion.h"
#include "grid_geometry.h"
#include "interpolation.h"
#include "nozzle.h"


/* ************************************************ */
void HaloOuterBoundary(const int side, const Data *d, int i, int j, int k, Grid *grid, int *touch) {/* Get primitives array for hot halo*/
/*
 * Set the values for the outer boundary conditions on sides that do not contain the jet.
 *
 * This routine is called from UserDefBoundary.
 *
 ************************************************** */

    // TODO: Also have a density and pressure condition for initiating touch.

    double halo_primitives[NVAR];
    double rho, prs, vx1, vx2, vx3, vmag;
    double rho_ini, prs_ini, vx1_ini, vx2_ini, vx3_ini, vmag_ini;
    double rho_cond, prs_cond, vel_cond;
    double fvel_small = 1.e-3;
    double frho_small = 1.e-3;
    double fprs_small = 1.e-3;
    int nv;

    /* These are the geometrical central points */
    double x1, x2, x3;
    double x1_inside, x2_inside, x3_inside;
    x1 = grid->x[IDIR][i];
    x2 = grid->x[JDIR][j];
    x3 = grid->x[KDIR][k];


    /* Get velocity, density, and prs from data, just inside the domain. */
    switch(side) {
        case X1_BEG:
            EXPAND(vx1 = d->Vc[VX1][k][j][IBEG];,
                   vx2 = d->Vc[VX2][k][j][IBEG];,
                   vx3 = d->Vc[VX3][k][j][IBEG];);
            rho = d->Vc[RHO][k][j][IBEG];
            prs = d->Vc[PRS][k][j][IBEG];
            x1_inside = grid->x[IDIR][IBEG];
            x2_inside = grid->x[JDIR][j];
            x3_inside = grid->x[KDIR][k];
            break;
        case X1_END:
            EXPAND(vx1 = d->Vc[VX1][k][j][IEND];,
                   vx2 = d->Vc[VX2][k][j][IEND];,
                   vx3 = d->Vc[VX3][k][j][IEND];);
            rho = d->Vc[RHO][k][j][IEND];
            prs = d->Vc[PRS][k][j][IEND];
            x1_inside = grid->x[IDIR][IEND];
            x2_inside = grid->x[JDIR][j];
            x3_inside = grid->x[KDIR][k];
            break;
        case X2_BEG:
            EXPAND(vx1 = d->Vc[VX1][k][JBEG][i];,
                   vx2 = d->Vc[VX2][k][JBEG][i];,
                   vx3 = d->Vc[VX3][k][JBEG][i];);
            rho = d->Vc[RHO][k][JBEG][i];
            prs = d->Vc[PRS][k][JBEG][i];
            x1_inside = grid->x[IDIR][i];
            x2_inside = grid->x[JDIR][JBEG];
            x3_inside = grid->x[KDIR][k];
            break;
        case X2_END:
            EXPAND(vx1 = d->Vc[VX1][k][JEND][i];,
                   vx2 = d->Vc[VX2][k][JEND][i];,
                   vx3 = d->Vc[VX3][k][JEND][i];);
            rho = d->Vc[RHO][k][JEND][i];
            prs = d->Vc[PRS][k][JEND][i];
            x1_inside = grid->x[IDIR][i];
            x2_inside = grid->x[JDIR][JEND];
            x3_inside = grid->x[KDIR][k];
            break;
        case X3_BEG:
            EXPAND(vx1 = d->Vc[VX1][KBEG][j][i];,
                   vx2 = d->Vc[VX2][KBEG][j][i];,
                   vx3 = d->Vc[VX3][KBEG][j][i];);
            rho = d->Vc[RHO][KBEG][j][i];
            prs = d->Vc[PRS][KBEG][j][i];
            x1_inside = grid->x[IDIR][i];
            x2_inside = grid->x[JDIR][j];
            x3_inside = grid->x[KDIR][KBEG];
            break;
        case X3_END:
            EXPAND(vx1 = d->Vc[VX1][KEND][j][i];,
                   vx2 = d->Vc[VX2][KEND][j][i];,
                   vx3 = d->Vc[VX3][KEND][j][i];);
            rho = d->Vc[RHO][KEND][j][i];
            prs = d->Vc[PRS][KEND][j][i];
            x1_inside = grid->x[IDIR][i];
            x2_inside = grid->x[JDIR][j];
            x3_inside = grid->x[KDIR][KEND];
            break;
        default:
            QUIT_PLUTO(1);
    }

    vmag = VMAG(x1, x2, x3, vx1, vx2, vx3);


    /* Hot halo primitives, just inside the domain. */
    HotHaloPrimitives(halo_primitives, x1_inside, x2_inside, x3_inside);
    rho_ini = halo_primitives[RHO];
    prs_ini = halo_primitives[PRS];
    EXPAND(vx1_ini = halo_primitives[VX1];,
           vx2_ini = halo_primitives[VX2];,
           vx3_ini = halo_primitives[VX3];);
    vmag_ini = VMAG(x1_inside, x2_inside, x3_inside, vx1_ini, vx2_ini, vx3_ini);

    /* Conditions */
    vel_cond = vmag_ini > 0. ? (vmag - vmag_ini) / vmag_ini : vmag;
    rho_cond = (rho - rho_ini) / rho_ini;
    prs_cond = (prs - prs_ini) / prs_ini;

    /* Set outflow (zero grad) if cell outside boundary has some velocity,
     * else set to Hot halo. */

    if ( ( (vel_cond > fvel_small) && (rho_cond > frho_small) && (prs_cond > fprs_small) ) || *touch == 1 ) {

        switch(side) {
            case X1_BEG:
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
                break;
            case X1_END:
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
                break;
            case X2_BEG:
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][JBEG][i];
                break;
            case X2_END:
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][JEND][i];
                break;
            case X3_BEG:
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][KBEG][j][i];
                break;
            case X3_END:
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][KEND][j][i];
                break;
            default:
                QUIT_PLUTO(1);
        }

        *touch = 1;
    }

    else

        /* Hot halo primitives, in the ghost zones. */
        HotHaloPrimitives(halo_primitives, x1, x2, x3);
        VAR_LOOP(nv) d->Vc[nv][k][j][i] = halo_primitives[nv];

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

#if (GRAV_POTENTIAL == NONE) || (GRAV_POTENTIAL == GRAV_HOMOGENEOUS)

    halo[RHO] = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
    halo[PRS] = PresIdealEOS(halo[RHO], g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP], MU_NORM);


    /* Hernquist potential (isothermal hydrostatic)*/
#elif GRAV_POTENTIAL == GRAV_HERNQUIST

    double r, a, rs, rho0, phi, phi0, te, c2;

    /* Thermodynamic parameters */
    rho0 = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];
    te = g_inputParam[PAR_HTMP] * ini_cgs[PAR_HTMP];
    c2 = g_gamma * CONST_kB * te  / (MU_CONST * CONST_amu) / vn.pot_norm;

    /* Distribution scale */
    a = g_inputParam[PAR_HRAD] * ini_code[PAR_HRAD];
    r = SPH1(x1, x2, x3);
    rs = r / a;

    /* Values for the potential */
    phi = BodyForcePotential(x1, x2, x3);
    phi0 = BodyForcePotential(0, 0, 0);

    halo[RHO] = rho0 * exp((phi - phi0) / c2);
    halo[PRS] = halo[RHO] * c2;


    /* This is also the default */

#elif defined(GRAV_TABLE)

    /* The density from the table is in units of cm^-3
     * The density can be scaled, by PAR_HRHO, but make sure PAR_HTMP
     * is the same as that used to generate the hot halo profile. */
    double r = SPH1(x1, x2, x3);
    halo[RHO] = InterpolationWrapper(hot_rad, hot_rho, hot_ndata, r);
    halo[RHO] *= g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];

    /* Currently, all gravity table profiles are assumed to be isothermal */
    halo[PRS] = PresIdealEOS(halo[RHO], g_inputParam[PAR_HTMP] * ini_code[PAR_HTMP], MU_NORM);

    /* Else we need the following line, and g_inputParam[PAR_HRHO] needs to be 1. */
//    halo[PRS] = InterpolationWrapper(hot_rad, hot_prs, hot_ndata, r);

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

#if RECONSTRUCT_4VEL == YES
    vel = VMAG(x1, x2, x3, halo[VX1], halo[VX2], halo[VX3]);
    scrh = Speed2Lorentz(vel);
    EXPAND(halo[VX1] *= scrh;, halo[VX2] *= scrh;, halo[VX3] *= scrh;);
#endif


    /* Tracers */
    halo[TRC] = 0.0;
#if CLOUDS
    halo[TRC + 1] = 0.0;
#endif


    /* A message for having initialized halo with potential*/
    if (!once01) {
        print("> Initializing hot halo distribution of type: \n\n");
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