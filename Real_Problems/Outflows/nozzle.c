//
// Created by Alexander Y. Wagner on 2019-01-02.
//

#include <pluto.h>
#include "nozzle.h"
#include "outflow.h"
#include "hot_halo.h"
#include "grid_geometry.h"
#include "pluto_usr.h"
#include "init_tools.h"

/* Global struct for nozzle */
Nozzle nz;

/* ************************************************ */
int InNozzleCap(const double x1, const double x2, const double x3) {
/*
 * Returns 1 if r is in nozzle cap, 0 if not.
 * Nozzle cap is a hemispherical region above the base of
 * the outflow cone.
 *
 ************************************************** */

// TODO: Debug this... doesn't seem to work anymore
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
    if (nz.is_two_sided) D_SELECT(cx1p, cx2p, cx3p) = fabs(D_SELECT(cx1p, cx2p, cx3p));

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

        /* Special case when dumping conservative quantitites and when nozzle is one-sided,
         * and aligned with FLOWDIR. We extend nozzle region along flowdir by rad. See SetNozzleGeometry. */
#if NOZZLE_FILL == NF_CONSERVATIVE
        if (!nz.is_two_sided && !(nz.dir > 0) D_SELECT(cx1p, cx2p, cx3p) -= nz.rad;
#endif


        /* Turn into cylindrical coords */
        double cr, cz;
        cr = CART2CYL1(cx1p, cx2p, cx3p);
        cz = CART2CYL2(cx1p, cx2p, cx3p);

        return (cr < nz.rad) && (cz < 0);
    }
}

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
    int is_fan, is_two_sided, is_halved, is_quarter;

    double small_angle = 1.e-12;
    double large = 1.e30;

    /* Get boundary conditions */
    int *left_bound = RuntimeGet()->left_bound;
    int *right_bound = RuntimeGet()->right_bound;
    int lb1 = left_bound[IDIR], rb1 = right_bound[IDIR];
    int lb2 = left_bound[JDIR], rb2 = right_bound[JDIR];
    int lb3 = left_bound[KDIR], rb3 = right_bound[KDIR];


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

    /* One of the most important quantities: the radius of the nozzle */
    rad = g_inputParam[PAR_ORAD] * ini_code[PAR_ORAD];


#if FEEDBACK_CYCLE == YES
    if (os.is_on == 0) {
        rad = 0;
    }
    else {

    /* Reduction of nozzle area via radius for deboosting.
       Accretion struct must be initialized beforehand */
        // TODO: always make ACCRETION TRUE if FEEDBACK_CYCLE is on
#if FBC_DEBOOST == YES && \
    (FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_0 || \
     FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_2 || \
     FBC_DEBOOST_MODE == FBC_DEBOOST_MODE_3)
        rad = sin(ang) * sqrt( ac.nzi.area * ac.deboost / (2. * CONST_PI * (1. - cos(ang))) ) ;
#endif

    }
#endif

    /* Quantities from input parameters */
    ang = g_inputParam[PAR_OANG] * ini_code[PAR_OANG];
    sph = g_inputParam[PAR_OSPH] * ini_code[PAR_OSPH];
    dir = g_inputParam[PAR_ODIR] * ini_code[PAR_ODIR];
    dbh = g_inputParam[PAR_ODBH] * ini_code[PAR_ODBH];
    omg = g_inputParam[PAR_OOMG] * ini_code[PAR_OOMG];
    phi = g_inputParam[PAR_OPHI] * ini_code[PAR_OPHI];

    /* Some useful derived quantitites */
    if (is_two_sided) {
        cbh = sqrt(sph * sph - rad * rad);
        orig = cbh;
    }
    else {
        cbh = (orig - dbh) / cos(dir) + rad * tan(dir);
        orig = g_domBeg[FLOWAXIS(IDIR, JDIR, KDIR)];
    }

    /* Cone geometry for the case of fan-like outflow.
     * Cone apex is only valid when cone is aligned with flow axis */
    if (is_fan) {
        cone_height = rad / tan(ang);
        cone_apex = orig - cone_height;
    }
    else {
        cone_height = large;
        cone_apex = -large;
    }


    /* Area and volume calculations --
     * The volume is the dump volume, into which mass, energy and momentum are dumped in the case of
     * NOZZLE_FILL == CONSERVATIVE. It is always the fully revolved volume.
     * Some special cases:
     * - Note, the volume is tricky to calculate in the case of one-sided fan-like tilted nozzles,
     *   so NOZZLE_FILL CONSERVATIVE is not supported for this case.
     * - In the case of one-sided, plane parallel jets, and NOZZLE_FILL CONSERVATIVE,
     *   we extend nozzle region along FLOWDIR by rad (see InNozzleRegion)
     *   */


    if (is_fan) {

        /* Area of nozzle */
        double cone_side = rad / sin(ang);
        area = 2. * CONST_PI * (1. - cos(ang)) * cone_side * cone_side;

        /* Volume of nozzle. It is the 3D (revolved) volume. */
        double spherical_cone_volume = area * cone_side / 3.;
        double overlap_volume, overlap_rad;

        if (is_two_sided) {
            overlap_rad = rad * (1. - cbh / cone_height);
            overlap_volume = CONST_PI * overlap_rad * overlap_rad * (cone_height - cbh) / 3.;
        }
        else {

            if (dir > 0) {
                vol = 0;
            }
            else {
                overlap_volume = CONST_PI * rad * rad * cone_height / 3.;
            }

        } // is two_sided
        vol = spherical_cone_volume - overlap_volume;


    }
    else {

        /* Area of nozzle */
        area = CONST_PI * rad * rad;

        /* Volume of nozzle */
        if (is_two_sided) {
            vol = area * cbh;
        }
        else {
            if (dir > 0) {
                vol = area * rad * tan(dir);
            }
            else {
                vol = area * rad;
            }
        }

    }

/* Remove inner spherical volume in the case of SPHERICAL geometry */
#if GEOMETRY == SPHERICAL
    if (is_two_sided) vol -= 4. * CONST_PI / 3. * (g_domBeg[IDIR] * g_domBeg[IDIR] * g_domBeg[IDIR]);
    else              vol -= 2. * CONST_PI / 3. * (g_domBeg[IDIR] * g_domBeg[IDIR] * g_domBeg[IDIR]);
#endif

    /* Power should always be the total power injected into box
     * (except for axis-symmetric cases, where volume and area still assume the full phi = 2 pi region) */
    // TODO: Make sure two-sided and one-sided cases are consistent.
    if (is_two_sided){
        area *= 2.;
        vol *= 2.;
    }

    /* -- End of area and volume calculation */


    /* Some consistency checks */

#if INTERNAL_BOUNDARY == YES
    /* Currently we require the BH location to be at 0,0,0  and noz.dbh = 0*/

    if (0 > dbh || dbh > 0) {
        print("Warning: If INTERNAL_BOUNDARY == YES, ODBH = 0. Setting noz.dbh to 0.");
        dbh = 0;
    }
    if (sph < rad) {
        print("Error: OSPH must be larger than ORAD");
        QUIT_PLUTO(1);
    }
    /* This check is for avoiding the cone of the nozzle to be buried in the
     outflow boundary for half-galaxy simulations. */
    if (!is_two_sided) {

        if (acos(1 - ((sph - cbh) * (sph - cbh) + rad * rad) /
                     (2 * sph * sph)) + dir > CONST_PI / 2) {
            print("Error: OSPH is too small. It must be at least...");
            QUIT_PLUTO(1);
            // TODO: Calculate minimum OSPH
        }
    }

#else // INTERNAL_BOUNDARY

    #if NOZZLE_FILL == NF_CONSERVATIVE

    print("Error: INTERNAL_BOUNDARY must be TRUE if NOZZLE_FILL == NF_CONSERVATIVE.");
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
    // TODO: Do we even need vol? If not, we don't need to calculate it
    noz->vol = vol;
    noz->cone_height = cone_height;
    noz->cone_apex = cone_apex;
    noz->is_fan = is_fan;
    noz->is_two_sided = is_two_sided;

}

/* ************************************************ */
void NozzleFill(Data *d, const Grid *grid) {
/*
 * Dumps mass, momentum, and energy into nozzle region,
 * according to nozzle parameters.
 *
 ************************************************** */

    RBox box;
    RBoxDefine (IBEG, IEND, JBEG, JEND, KBEG, KEND, CENTER, &box);

    // Just for initialization, need to convert primitives to conservatives
    // TODO: This may not be needed if this is modularized and attached somewhere else in main.c
    if (g_time <= 0.) PrimToCons3D(d->Vc, d->Uc, &box);

    double *x1, *x2, *x3;
    double out_conservatives[NVAR];
    int nv;

    /* zero buffer */
    VAR_LOOP(nv) out_conservatives[nv] = 0.;

    /* These are the geometrical central points */
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];


    // TODO: Relativistic version
#if (PHYSICS == RHD) || (PHYSICS == RMHD)
    print("Error: NOZZLE_FILL CONSERVATIVE not supported if PHYSICS is RHD or RMHD.");
    QUIT_PLUTO(1);
#endif

    // TODO: Adjust volume for the case below
    // NOTE: Volume is incorrectly calculated for spherical setups; nothing is dumped in the excluded x1_beg region.
//#if GEOMETRY == SPHERICAL
//    if (g_domBeg[IDIR] > 0.) {
//        print("Error: NOZZLE_FILL CONSERVATIVE not supported if (GEOMETRY == SPHERICAL) && (g_domBeg[IDIR] > 0.).");
//        QUIT_PLUTO(1);
//    }
//#endif

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
    ConsToPrim3D(d->Uc, d->Vc, d->flag, &box);
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

/*********************************************************************** */
void NozzleVolume(Data *d, Grid *grid) {
/*!
 * Calculate nozzle volume by summing cell volumes. This is the volume
 * into which mass, momentum, and energy are dumped
 *
 *********************************************************************** */

    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];

    int i, j, k;
    double vol, vol_a = 0;

    DOM_LOOP(k, j, i) {

                if (! InNozzleRegion(x1[i], x2[j], x3[k])) vol += ElevateVolume(grid->dV[k][j][i]);

            }

#ifdef PARALLEL
    MPI_Allreduce(&vol_a, &vol, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    vol = vol_a;

#endif
    nz.vol_counted = vol;
}



/*********************************************************************** */
void InitDomainNozzle(Data *d, Grid *grid){
/*!
 * Initilize the nozzle in the domain with primitive variables.
 * This is done, after initializing clouds and the hot phase.
 *
 *
 * Initialize nozzle if we're in nozzle inlet region, otherwise halo
 *
 *********************************************************************** */

    double *x1 = grid->x[IDIR];
    double *x2 = grid->x[JDIR];
    double *x3 = grid->x[KDIR];

    double out_primitives[NVAR];
    double halo_primitives[NVAR];

    int i, j, k, iv;

    TOT_LOOP(k, j, i) {
        
                if (InNozzleRegion(x1[i], x2[j], x3[k]) || InNozzleCap(x1[i], x2[j], x3[k])) {

                    OutflowPrimitives(out_primitives, x1[i], x2[j], x3[k]);
                    HotHaloPrimitives(halo_primitives, x1[i], x2[j], x3[k]);

                    NVAR_LOOP(iv) {
                        d->Vc[iv][k][j][i] = halo_primitives[iv] +
                                             (out_primitives[iv] - halo_primitives[iv]) * Profile(x1[i], x2[j], x3[k]);
                    }
                }
            }

}


/* ************************************************************** */
void ClearNozzleSurrounding(double *cell, const double *halo, const double x1, const double x2, const double x3) {
/*!
 *  This routine sets the value of the region around the center to hot-halo values,
 *  with a smooth transition. This is useful when switching on the jet midway through
 *  a simulation in cases where the accretion is very high.
 *
 *  Note that CLOUD_EXTRACT_CENTRAL_BUFFER must be YES for this routine to be called.
 *
**************************************************************** */

    double tanhfactor;

    /* Buffer factor around osph */
    double rad  = SPH1(x1, x2, x3);
    double incf = 1.3;    // Radius of central buffer (in units of osph)
    double wsmf = 0.3;    // Width of smoothing region (in units of osph)

    /* Inner hemisphere to keep free */
    double osph = g_inputParam[PAR_OSPH] * ini_code[PAR_OSPH];
    double circ = rad / (incf * osph);

    /* The smoothing region - smooth in logarithm of density ratio*/
    double offset = circ - 1.;
    offset = MIN(offset, 0.5 * wsmf * osph);
    offset = MAX(offset, -0.5 * wsmf * osph);
    tanhfactor = tanh(tan(CONST_PI * offset / (wsmf * osph)));

    double ratio;
    int iv;
    VAR_LOOP(iv) {
        /* Use logarithmic smoothing if possible, else linear smoothing */
        if (cell[iv] * halo[iv] > 0) {
            ratio = cell[iv] / halo[iv];
            cell[iv] = halo[iv] * pow(10, 0.5 * log10(ratio) + 0.5 * tanhfactor * log10(ratio));
        }
        else {
            cell[iv] = 0.5 * (cell[iv] + halo[iv]) + 0.5 * tanhfactor * (cell[iv] - halo[iv]);
        }
    }

}
