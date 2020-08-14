//
// Created by Alexander Y. Wagner on 4/6/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "clouds.h"
#include "idealEOS.h"
#include "hot_halo.h"
#include "interpolation.h"
#include "read_grav_table.h"
#include "init_tools.h"
#include "outflow.h"
#include "grid_geometry.h"
#include "sfr.h"
#include "io_tools.h"
#include "accretion.h"
#include "nozzle.h"

/* Global struct for cloud analytics */
CloudAnalytics ca;

#if CLOUDS == YES

/* ************************************************************** */
int CloudCubePixel(int *el, const double x1,
                   const double x2,
                   const double x3)
/*!
 * Locates cloud cube pixel.
 *
 * \param [out] el  array of size 3 holding cube coordinates
 * \param [in] x1   current zone x1 coordinates (any geometry)
 * \param [in] x2   current zone x2 coordinates (any geometry)
 * \param [in] x3   current zone x3 coordinates (any geometry)
 *
 * Calculate pixel coordinates (here el[]) in fractal cube
 * for given zone. The origin of the pixel coordinates is at
 * the corner with the smallest values of x1, x2, and x3
 * (the lower corner).
 *
 * Returns 1 if we are in a zone that is covered by the
 * fractal cube and 0 otherwise.  *
 *
 **************************************************************** */
{

    double x, y, z;
    double xfrac, yfrac, zfrac;
    double csz1, csz2, csz3;

    /* Convert to cartesian geometry, because clouds cube
     * is always assumed to be in cartesian. */
    D_EXPAND(x = CART1(x1, x2, x3);,
             y = CART2(x1, x2, x3);,
             z = CART3(x1, x2, x3););

    /* The cloud region in simulation domain as obtained
     * in input_data.c */
    D_EXPAND(csz1 = g_idBoxEnd[0] - g_idBoxBeg[0];,
             csz2 = g_idBoxEnd[1] - g_idBoxBeg[1];,
             csz3 = g_idBoxEnd[2] - g_idBoxBeg[2];);

    /* Check whether we're in the fractal cube box */
    D_EXPAND(xfrac = (x - g_idBoxBeg[0]) / csz1;,
             yfrac = (y - g_idBoxBeg[1]) / csz2;,
             zfrac = (z - g_idBoxBeg[2]) / csz3;);

    if (!(D_EXPAND(xfrac > 1., || yfrac > 1., || zfrac > 1.) ||
          D_EXPAND(xfrac < 0., || yfrac < 0., || zfrac < 0.))) {

        /* Fill cube pixel element array */
        D_EXPAND(el[0] = (int) (xfrac * g_idnx1);,
                 el[1] = (int) (yfrac * g_idnx2);,
                 el[2] = (int) (zfrac * g_idnx3););
        return 1;
    }

    else {
        return 0;
    }

}

/* ************************************************************** */
void NormalizeFractalData(double *cloud, const double x1, const double x2, const double x3)
/*!
 *  This routine is just used to noramlize the data
 *  The input density cube has a normalized lognormal density
 *  distribution about 1 (which is later multiplied by WRHO)
 *  and the input velocity cube has a mean of 0
 *  and is assumed to be in units of km / s.
 *
 **************************************************************** */
{

    /* Normalize velocities to code units here */
    EXPAND(cloud[VX1] *= ini_code[PAR_WTRB];,
           cloud[VX2] *= ini_code[PAR_WTRB];,
           cloud[VX3] *= ini_code[PAR_WTRB];);

}

/* ************************************************************** */
void CloudDensity(double *cloud, const double x1, const double x2, const double x3)
/*!
 * Multiply cloud fractal factor (currently in *cloud) with a
 * desired mean density profile. Currently CLOUD_DENSITY =
 *     - CD_HERNQUIST
 *     - CD_KEPLERIAN
 *     - CD_HOMOGENEOUS
 * The default is also CD_HOMOGENEOUS. This is the apodization step. All
 * cells are still cloud material.
 *
 **************************************************************** */
{

    int il;
    double r1, r2, y0, y1, y2, y3, frac;
    double r_sph, r_cyl, z_cyl, phi_rz, phi_r0, phi_00;
    double dens, sigma_g2, wrho, wrho_cgs, wrad_cgs, wtrb, ek2;
    double wtrb_cgs;
    static int once01 = 0;



    /* The following are only some
     * profiles for the warm phase that produce
     * reasonable ratios of mean warm phase density and
     * hot phase density in large parts of the domain. */
#if CLOUD_DENSITY == CD_HERNQUIST

    double a, rs;
    r = SPH1(x1, x2, x3);
    a = g_inputParam[PAR_WRAD] * ini_code[PAR_WRAD];
    rs = r / a;
    wrho = g_inputParam[PAR_WRHO] * ini_code[PAR_WRHO];
    dens = wrho / (rs * pow((1 + rs), 3));


#elif CLOUD_DENSITY == CD_KEPLERIAN

    /* The dense gas mean density input parameter (in code units)*/
    wrho = g_inputParam[PAR_WRHO] * ini_code[PAR_WRHO];

    /* The global dense gas velocity dispersion */
#if CLOUD_SCALE == CS_SCALE_HEIGHT
    wrad_cgs = g_inputParam[PAR_WRAD] * ini_cgs[PAR_WRAD];
    wrho_cgs = g_inputParam[PAR_WRHO] * ini_cgs[PAR_WRHO];
    sigma_g2 = 4 * CONST_PI * CONST_G * wrho_cgs * wrad_cgs * wrad_cgs / 9.;

    /* Total dens gas velocity dispersion (in code units) */
    sigma_g2 /= vn.v_norm * vn.v_norm;

#elif CLOUD_SCALE == CS_VELOCITY_DISPERSION
    /* External initialisation of turbulent velocity dispersion */
    wtrb = g_inputParam[PAR_WTRB] * ini_code[PAR_WTRB];
    sigma_g2 = wtrb * wtrb;

    /* Print out scale height of clouds */
    if (!once01) {
        wtrb_cgs = g_inputParam[PAR_WTRB] * ini_cgs[PAR_WTRB];
        wrho_cgs = g_inputParam[PAR_WRHO] * ini_cgs[PAR_WRHO];
        wrad_cgs = 1.5 * wtrb_cgs / sqrt(CONST_PI * CONST_G * wrho_cgs);
        print("  Cloud scale height is %12.6f pc.\n\n", wrad_cgs / (1000. * CONST_pc));

        once01 = 1;
    }

#endif

    /* Gravitational potential */
    phi_rz = BodyForcePotential(x1, x2, x3);
//    phi_00 = BodyForcePotential(0, 0, 0); // Cannot interpolate this point, if central BH is present.
    phi_00 = 0.;

    /* Is the potential non-spherical? */
    ek2 = g_inputParam[PAR_WROT] * g_inputParam[PAR_WROT];

    if (ek2 > 0) {

        /* Now, the same for the cylindrical radius (in cgs) */
        r_cyl = fabs(CYL1(x1, x2, x3));
        phi_r0 = BodyForcePotential(r_cyl, 0, 0);
    }
    else {
        phi_r0 = 0;

    }

    /* The profile */
    // TODO: Check if this still works if a BH potential is included.
    dens = wrho * exp((-phi_rz + phi_r0 * ek2 + phi_00 * (1. - ek2)) / sigma_g2);

#elif CLOUD_DENSITY == CD_MILKY_WAY_PJM

    /* Parameters of McMillan model for HI disc, without central hole */
    double Rd = 7. * CONST_pc * 1.e3 / vn.l_norm;
    double zd = 0.085 * CONST_pc * 1.e3 / vn.l_norm;
    double Sigma0 = 53.1 * CONST_Msun / (CONST_pc * CONST_pc) / (vn.dens_norm * vn.l_norm);

    r_cyl = fabs(CYL1(x1, x2, x3));
    z_cyl = fabs(CYL2(x1, x2, x3));

    /* It is assumed that cloud[RHO] = 1. for this setup. */
    double sech = 1. / cosh(z_cyl / (2. * zd));
    dens = Sigma0 / (4. * zd) * exp(-r_cyl / Rd) * sech * sech;

#elif CLOUD_DENSITY == CD_HOMOGENEOUS
    /* Homogeneous halo density, but with gravity */
    dens = g_inputParam[PAR_WRHO] * ini_code[PAR_WRHO];

#else
    /* Homogeneous halo, gravity off. */
    dens = g_inputParam[PAR_WRHO] * ini_code[PAR_WRHO];

#endif

    /* Apodize! (code units) */
    cloud[RHO] = cloud[RHO] * dens;

}

/* ************************************************************** */
double CloudExtractEllipsoid(double fdratio, const double x1, const double x2, const double x3) {
/*!
 *
**************************************************************** */

    static int once01 = 0;
    double r_cyl, rz;
    double tanhfactor;
    double ellipse;
    double wrot, wrad, wsmf;


    /* Grid point in cartesian coordinates */
    double cx1, cx2, cx3;
    cx1 = CART1(x1, x2, x3);
    cx2 = CART2(x1, x2, x3);
    cx3 = CART3(x1, x2, x3);

    double cx1p, cx2p, cx3p;
    RotateGrid2Disc(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);

    /* The distances in physical space */
    r_cyl = CYL1(cx1p, cx2p, cx3p);
    rz    = CYL2(cx1p, cx2p, cx3p);

    /* The rotational parameter */
    wrot = g_inputParam[PAR_WROT] * ini_code[PAR_WROT];

    /* Always just use WRAD parameter here rather than velocity dispersion,
     * so that we can control the apodization and ellipsoid extraction separately. */
    wrad = g_inputParam[PAR_WRAD] * ini_code[PAR_WRAD];

    /* ellipse = 1 is the region beyond which the extraction (with smoothing) should take place */
    ellipse = (r_cyl * r_cyl * (1 - wrot * wrot) + rz * rz) / (wrad * wrad);

    /* Smoothing region scale, hardcoded here. */
    /* For discs a factor 1/e was good, for spherical distributions, 0.1.
     * We interpolate between these. */
    wsmf = 0.1 + (exp(-1) - 0.1) * wrot;

    /* The smoothing region - smooth in logarithm of density ratio*/
    double offset = sqrt(ellipse) - 1.;
    offset = MIN(offset, 0.5 * wsmf);
    offset = MAX(offset, -0.5 * wsmf);
    tanhfactor = tanh(tan(-CONST_PI * offset / wsmf));
    fdratio = pow(10, 0.5 * log10(fdratio) + 0.5 * tanhfactor * log10(fdratio));


    /* Take into account tilt of disc in extraction */
    double dir = g_inputParam[PAR_WDIR] * ini_code[PAR_WDIR];

    /* Add additional cylindrical extraction, if disc radial scale > half domain extent */
    if (wrad * wrad / (1 - wrot * wrot) / cos(dir) > g_domEnd[IDIR] * g_domEnd[IDIR]) {

        /* Change meaning of wrad to mean half domain extent here.
         * But take into account rotation too - the disc height will go out of the grid, otherwise */
        double cylinder = r_cyl / (g_domEnd[IDIR] * cos(dir));

        double offset = cylinder - 1.;
        offset = MIN(offset, 0.5 * wsmf);
        offset = MAX(offset, -0.5 * wsmf);
        tanhfactor = tanh(tan(-CONST_PI * offset / wsmf));
        fdratio = pow(10, 0.5 * log10(fdratio) + 0.5 * tanhfactor * log10(fdratio));

    }


    if (!once01){
      print("> Cloud extraction: CLOUD_EXTRACTION_ELLIPSOID.\n\n");
      once01 = 1;
    }

    return fdratio;

}

/* ************************************************************** */
double CloudExtractCentralBuffer(double fdratio, const double x1, const double x2, const double x3) {
/*!
 *
**************************************************************** */

    double tanhfactor;

    /* Buffer factor around osph */
    double rad  = SPH1(x1, x2, x3);
    double incf = 2.0;    // Radius of central buffer (in units of osph)
    double wsmf = 1.0;    // Width of smoothing region

    /* Inner hemisphere to keep free */
    double osph = g_inputParam[PAR_OSPH] * ini_code[PAR_OSPH];
    double circ = rad / (incf * osph);

    /* The smoothing region - smooth in logarithm of density ratio*/
    double offset = circ - 1.;
    offset = MIN(offset, 0.5 * wsmf);
    offset = MAX(offset, -0.5 * wsmf);
    tanhfactor = tanh(tan(CONST_PI * offset / wsmf));
    fdratio = pow(10, 0.5 * log10(fdratio) + 0.5 * tanhfactor * log10(fdratio));

    return fdratio;
}


/* ************************************************************** */
int CloudExtract(double *cloud,
                 const double *halo,
                 const int *pixel,
                 const double x1, const double x2, const double x3)
/*!
 *
 * This function extracts the cloud from the box in different ways,
 * as selected by the CLOUD_EXTRACT. Extraction is any method except
 * the generation of porosity through the critical temperature.
 *
 * Here, cloud[RHO] contains the cloud density but for convenience,
 * we work with fdratio, the ratio of apodized cloud to halo
 * density. After this routine cloud contains the
 * updated cloud value.
 *
 * The function returns True if pixel is a cloud cell and False if
 * pixel is not a cloud cell, but has instead been replaced with a
 * halo cell due to the extraction process
 *
 **************************************************************** */
{


    double fdratio;

    fdratio = cloud[RHO] / halo[RHO];

#if CLOUD_EXTRACT_ELLIPSOID == TRUE
    fdratio = CloudExtractEllipsoid(fdratio, x1, x2, x3);
#endif

#if CLOUD_EXTRACT_CENTRAL_BUFFER == TRUE
    fdratio = CloudExtractCentralBuffer(fdratio, x1, x2, x3);
#endif


    /* The ratio of cloud to halo density for a cell
     * shouldn't really ever be < 1. */
    fdratio = MAX(fdratio, 1.);

    /* Set to cloud density */
    cloud[RHO] = fdratio * halo[RHO];

    return fdratio > 1. ? 1 : 0;
}

/* ************************************************************** */
void CloudVelocity(double *cloud, double *halo,
                   const double x1, const double x2, const double x3)
/*!
 * This function fills in the cloud velocity for the clouds primitives
 * array.
 *   - CV_KEPLERIAN: read in virial velocity and add mean keplerian vphi
 *   - CV_ZERO:: Set all velocities to zero
 *
 * In addition, constant velocities are applied according to the values of
 * runtime parameters:
 *   - PAR_WVRD   Radial velocity (km/s)
 *   - PAR_WVPL   Parallel velocity (km/s)
 *   - PAR_WVPP   Perpendicular velocity (km/s)
 *   - PAR_WVAN   Angular velocity (km/s)
 *
 **************************************************************** */
{

    double v1, v2, v3;

#if CLOUD_VELOCITY == CV_KEPLERIAN

    double ek;
    double r_cyl, rad, r1, r2, frac;
    double y0, y1, y2, y3, dphidr;
    double gvec[COMPONENTS];
    int il;
    double vpol1, vpol2, vpol3;
    double xpol1, xpol2, xpol3;


    /* The rotational parameter */
    ek = g_inputParam[PAR_WROT] * ini_code[PAR_WROT];

    EXPAND(v1 = cloud[VX1];,
           v2 = cloud[VX2];,
           v3 = cloud[VX3];);

    /* If ek == 0, then it's spherical) and velocities are just those read in.
     * If ek > 0., there is a Keplerian component that needs to be added. */
    if (ek > 0.) {

        /* Grid point in cartesian coordinates */
        double cx1, cx2, cx3;
        cx1 = CART1(x1, x2, x3);
        cx2 = CART2(x1, x2, x3);
        cx3 = CART3(x1, x2, x3);

        double cx1p, cx2p, cx3p;
        RotateGrid2Disc(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);

        /* Convert coordinates and velocity vectors to cylindrical polars */
        xpol1 = POL1(cx1p, cx2p, cx3p);
        xpol2 = POL2(cx1p, cx2p, cx3p);
        xpol3 = POL3(cx1p, cx2p, cx3p);

        vpol1 = VPOL1(cx1p, cx2p, cx3p, v1, v2, v3);
        vpol2 = VPOL2(cx1p, cx2p, cx3p, v1, v2, v3);
        vpol3 = VPOL3(cx1p, cx2p, cx3p, v1, v2, v3);

        /* TODO: (NOTE) it is not guaranteed that the potential is also rotated */
        /* Use local gradient, rather than mid plane potential.
         * At least for N-body initializations of thick discs in axisymmetric potentials,
         * this gives a better results, according to Miki et al 2017, MAGI paper. */
        r_cyl = CYL1(cx1p, cx2p, cx3p);
        BodyForceVector(cloud, gvec, cx1p, cx2p, cx3p);  // Acceleration vector (pointing inward)
        dphidr = -VPOL1(cx1p, cx2p, cx3p, gvec[IDIR], gvec[JDIR], gvec[KDIR]);

        /* The angular (linear) velocity */
        // TODO: Rather than use ek, deduct turbulent and thermal support from Keplerian velocity.
        //       Turbulent support  = sigma_g * r / rd where rd is the disc radial scale.
        //       This should also be possible for isothermal turbulent discs, i.e.,
        //       one doesn't need to specify WROT (discuss this with Geoff).
        //       But what is the scale length of an isothermal disc, then?
        //       Maybe just  zd / (1. - ek), where zd is the scale height of the disc.
        //       There are also some typos in Geoff's derivation for the MW_PJM case,
        //       for which the assumption that T is constant is also not correct.

//#if CLOUD_DENSITY == CD_MILKY_WAY_PJM
//        double wtrb = g_inputParam[PAR_WTRB] * ini_code[PAR_WTRB];
//        double wth2 = halo[PRS] / halo[RHO];
//
//        vpol2 += sqrt(r_cyl * dphidr - wtrb * wtrb - wth2);
//#else
        vpol2 += ek * sqrt(r_cyl * dphidr);
//#endif

        /* Get vector in cartesian coordinates */
        double vcx1p, vcx2p, vcx3p;
        EXPAND(vcx1p = VPOL2CART1(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3);,
               vcx2p = VPOL2CART2(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3);,
               vcx3p = VPOL2CART3(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3););

        /* Rotate vector back */
        double cv1, cv2, cv3;
        RotateDisc2Grid(vcx1p, vcx2p, vcx3p, &cv1, &cv2, &cv3);

        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(v1 = VCART_1(cx1, cx2, cx3, cv1, cv2, cv3);,
               v2 = VCART_2(cx1, cx2, cx3, cv1, cv2, cv3);,
               v3 = VCART_3(cx1, cx2, cx3, cv1, cv2, cv3););

    }


#elif (CLOUD_VELOCITY == CV_ZERO) || (CLOUD_VELOCITY == NONE)
    EXPAND(v1 = 0;,
           v2 = 0;,
           v3 = 0;);

#endif


    double vcart1, vcart2, vcart3;
    double xcart1, xcart2, xcart3;
    double vsph1, vsph2, vsph3;
    double xsph1, xsph2, xsph3;
    double vcyl1, vcyl2, vcyl3;
    double xcyl1, xcyl2, xcyl3;

    if (fabs(g_inputParam[PAR_WVRD]) > 0) {

        /* Convert coordinates and velocity vectors to spherical */

        xsph1 = SPH1(x1, x2, x3);
        xsph2 = SPH2(x1, x2, x3);
        xsph3 = SPH3(x1, x2, x3);

        vsph1 = VSPH1(x1, x2, x3, v1, v2, v3);
        vsph2 = VSPH2(x1, x2, x3, v1, v2, v3);
        vsph3 = VSPH3(x1, x2, x3, v1, v2, v3);

        /* Apply change to radial component (assumed to be in km/s) */
        vsph1 += g_inputParam[PAR_WVRD] * ini_code[PAR_WVRD];

        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(v1 = VSPH_1(xsph1, xsph2, xsph3, vsph1, vsph2, vsph3);,
               v2 = VSPH_2(xsph1, xsph2, xsph3, vsph1, vsph2, vsph3);,
               v3 = VSPH_3(xsph1, xsph2, xsph3, vsph1, vsph2, vsph3););

    }

    if (fabs(g_inputParam[PAR_WVPL]) > 0) {

        /* Convert coordinates and velocity vectors to cartesian */

        xcart1 = CART1(x1, x2, x3);
        xcart2 = CART2(x1, x2, x3);
        xcart3 = CART3(x1, x2, x3);


        vcart1 = VCART1(x1, x2, x3, v1, v2, v3);
        vcart2 = VCART2(x1, x2, x3, v1, v2, v3);
        vcart3 = VCART3(x1, x2, x3, v1, v2, v3);

        /* Apply change to component parallel to flow axis (assumed to be in km/s).
         * Can't use FLOWAXIS macro though because we are in transformed coords. */
        D_SELECT(vcart1, vcart2, vcart3) += g_inputParam[PAR_WVPL] * ini_code[PAR_WVPL];

        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(v1 = VCART_1(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3);,
               v2 = VCART_2(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3);,
               v3 = VCART_3(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3););

    }

    if (fabs(g_inputParam[PAR_WVPP]) > 0) {

        /* Convert coordinates and velocity vectors to cartesian */

        xcart1 = CART1(x1, x2, x3);
        xcart2 = CART2(x1, x2, x3);
        xcart3 = CART3(x1, x2, x3);

        vcart1 = VCART1(x1, x2, x3, v1, v2, v3);
        vcart2 = VCART2(x1, x2, x3, v1, v2, v3);
        vcart3 = VCART3(x1, x2, x3, v1, v2, v3);

        /* Apply change to component perpendicular to flow axis (assumed to be in km/s).
        * Can't use FLOWAXIS macro though because we are in transformed coords. */
        D_SELECT(vcart1, vcart1, vcart2) += g_inputParam[PAR_WVPP] * ini_code[PAR_WVPP];

        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(v1 = VCART_1(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3);,
               v2 = VCART_2(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3);,
               v3 = VCART_3(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3););

    }

    if (fabs(g_inputParam[PAR_WVAN]) > 0) {

        /* Convert coordinates and velocity vectors to cylindrical */

        xcyl1 = POL1(x1, x2, x3);
        xcyl2 = POL2(x1, x2, x3);
        xcyl3 = POL3(x1, x2, x3);

        vcyl1 = VPOL1(x1, x2, x3, v1, v2, v3);
        vcyl2 = VPOL2(x1, x2, x3, v1, v2, v3);
        vcyl3 = VPOL3(x1, x2, x3, v1, v2, v3);

        /* Apply change to polar component (assumed to be in km/s) */
        vcyl2 += g_inputParam[PAR_WVAN] * ini_code[PAR_WVAN];

        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(v1 = VPOL_1(xcyl1, xcyl2, xcyl3, vcyl1, vcyl2, vcyl3);,
               v2 = VPOL_2(xcyl1, xcyl2, xcyl3, vcyl1, vcyl2, vcyl3);,
               v3 = VPOL_3(xcyl1, xcyl2, xcyl3, vcyl1, vcyl2, vcyl3););

    }


    /* Put back in to cloud array */
    EXPAND(cloud[VX1] = v1;, cloud[VX2] = v2;, cloud[VX3] = v3;);

}

/* ************************************************************** */
int CloudPrimitives(double *cloud,
                    const double x1, const double x2, const double x3)
/*!
 * This function returns 1 if cell is cloud material, 0 if not.
 *
 * The array *cloud is filled in the process that involves three steps:
 * 1) Check that we're inside fractal cube domain. CloudCubePixel.
 * 2) Get the cloud values for the current coordinates. GetFractalData.
 * 3) Apodization step. Currently supported are
 *    CLOUD_DENSITY =
 *        - CD_HOMOGENEOUS
 *        - CD_HERNQUIST
 *        - CD_KEPLERIAN
 *    The default is also CD_HOMOGENEOUS. All cells are still cloud material.
 *    We prefer to work with the density ratio, rho_cloud/rho_halo, mainly
 *    for the next extraction part. It gives us a good handle on how close
 *    we are to the halo density.
 * 4) Use a form of cloud extraction by calling function CloudExtract.
 *    Some cells will not be cloud material anymore.
 * 5) Select, based on temperature criterion, whether cell is cloud material
 *    or not.
 *
 *
 **************************************************************** */
{

    double halo[NVAR], vel[COMPONENTS], scrh;
    int nv, cube_pixel[DIMENSIONS];
    int is_cloud = 0;

    /* Get cloud pixel coordinates */
    if (CloudCubePixel(cube_pixel, x1, x2, x3)) {

        /* Get the fractal factor for this cell */
        NormalizeFractalData(cloud, x1, x2, x3);

        /* Cloud density profile. Apodize with a mean density profile */
        CloudDensity(cloud, x1, x2, x3);

        /* Cloud velocity. */
        CloudVelocity(cloud, halo, x1, x2, x3);

        /* Extract w.r.t. hot halo. Halo primitives are required
         * for this step. */
        HotHaloPrimitives(halo, x1, x2, x3);

        if (CloudExtract(cloud, halo, cube_pixel, x1, x2, x3)) {

            /* Cloud pressure
             * Underpressure the clouds slightly, so that they don't emit sound waves. */
            cloud[PRS] = halo[PRS] * CLOUD_UNDERPRESSURE;

            /* Tracers */
            cloud[TRC] = 0.0;
            cloud[TRC + 1] = 1.0;

            /* Final test - is cloud pixel thermally stable? */
            is_cloud = WarmTcrit(cloud);

        }
    }

    /* Fill cloud array with halo primitves if not a cloud cell. This is not
     * strictly necessary, since it is done outside CloudPrimitives, but we
     * do it anyway for completeness. */
    if (is_cloud == 0) NVAR_LOOP(nv) cloud[nv] = halo[nv];

    return is_cloud;
}

/* ********************************************************************* */
int WarmTcrit(double *const warm)
/*!
 * This routine returns 1 if cloud is still cloud, after thermal
 * stability criterion was set.
 *
 *********************************************************************** */
{
    double mu, wtemp;
    int nv;

#ifndef CLOUD_MUCRIT
    fputs("Error: CLOUD_MUCRIT not defined.\n", stderr); QUIT_PLUTO(1);
#endif

#ifndef CLOUD_TCRIT
    fputs("Error: CLOUD_TCRIT not defined.\n", stderr); QUIT_PLUTO(1);
#endif

    /* Only a cloud pixel if wtemp is below critical temperature
     * of thermal instability */
    if (warm[PRS] / warm[RHO] > CLOUD_TCRIT / CLOUD_MUCRIT / KELVIN) return 0;
    else return 1;

}


#endif


// TODO: Remove the tracer dependence in some of the below

/* ********************************************************************* */
void CloudAnalysis(Data *d, Grid *grid) {
/*!
 * Analysis of warm phase.
 *
 *********************************************************************** */

    /* Radii at which to measure outflow rates */
    ca.nrad = 3;
    ca.radii[0] = 0.1;
    ca.radii[1] = 0.3;
    ca.radii[2] = 1.;

    WarmOutflowRates(d, grid);
    WarmPhaseMass(d, grid);
    WarmPhasePorosity(d, grid);

}

/* ********************************************************************* */
void WarmOutflowRates(Data *d, Grid *grid) {
/*!
 * Calculate warm phase mass outlfow rate, outflow power, and outflow momentum.
 *
 *********************************************************************** */

    double rho, vs1, tr2;
    double vx1, vx2, vx3;
    double *x1, *x2, *x3;
    double edot, pdot, mdot, vrho, wgt;
    double edot_a, pdot_a, mdot_a, vrho_a, wgt_a;


    /* Determine number of sampling points to use */
    int npoints;
    double oversample = 30;
    double area, area_per_point;
    double dl_min = grid->dl_min[IDIR];

    double radius;

    for (int irad = 0; irad < ca.nrad; irad ++) {

        radius = ca.radii[irad];

        /* Area through which accretion rate is measured.
         * If nozzle is not two-sided, the accretion rate represents
         * the accretion rate in one half of the galaxy. */
        area = 4 * CONST_PI * radius * radius;
        if (!nz.is_two_sided) area /= 2.;

        npoints = (int) (area / (dl_min * dl_min) * oversample);
#if DIMENSIONS == 2
        npoints = (int) sqrt(npoints);
#endif

        /* Get npoint points on spherical surface */
        x1 = ARRAY_1D(npoints, double);
        x2 = ARRAY_1D(npoints, double);
        x3 = ARRAY_1D(npoints, double);
        npoints = UniformSamplingSphericalSurface(npoints, radius, x1, x2, x3);
        area_per_point = area / npoints;

        /* Determine which variables to interpolate */
        double v[NVAR];
        int vars[] = {RHO, ARG_EXPAND(VX1, VX2, VX3), TRC+1, -1};

        /* Zero cumulative quantities */
        edot_a = pdot_a = mdot_a = vrho_a = wgt_a = 0;

        /* Calculate outflow rate at every point, if it is in the local domain */
        for (int ipoint = 0; ipoint < npoints; ipoint++) {

            if (PointInDomain(grid, x1[ipoint], x2[ipoint], x3[ipoint])) {

                InterpolateGrid(d, grid, vars, x1[ipoint], x2[ipoint], x3[ipoint], v);

                /* Calculate and sum outflow rate */
                rho = v[RHO];
                vx1 = vx2 = vx3 = 0;
                EXPAND(vx1 = v[VX1];,
                       vx2 = v[VX2];,
                       vx3 = v[VX3];);
                vs1 = VSPH1(x1[ipoint], x2[ipoint], x3[ipoint], vx1, vx2, vx3);
                vs1 = fabs(MAX(vs1, 0));
                tr2 = v[TRC+1];

                wgt = tr2 * rho;
                vrho = wgt * vs1;
                mdot = vrho * area_per_point;
                pdot = vrho * vs1 * area_per_point;
                edot = 0.5 * vrho * vs1 * vs1 * area_per_point;

                wgt_a += wgt;
                vrho_a += vrho;
                mdot_a += mdot;
                pdot_a += pdot;
                edot_a += edot;

            }

        }

#ifdef PARALLEL
        MPI_Allreduce(&wgt_a,  &wgt,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&vrho_a, &vrho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mdot_a, &mdot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&pdot_a, &pdot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&edot_a, &edot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

        wgt = wgt_a;
        vrho = vrho_a;
        mdot = mdot_a;
        pdot = pdot_a;
        edot = edot_a;

#endif

        vrho /= wgt;

        ca.vrho_out[irad] = vrho;
        ca.mdot_out[irad] = mdot;
        ca.pdot_out[irad] = pdot;
        ca.edot_out[irad] = edot;

        FreeArray1D(x1);
        FreeArray1D(x2);
        FreeArray1D(x3);

    }


    /* Total domain based estimates */

    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    double cell_vol, cell_area, rad;

    /* Zero cumulative quantities */
    edot_a = pdot_a = mdot_a = vrho_a = wgt_a = 0;

    int i, j, k;
    DOM_LOOP(k, j, i) {

                if (! InSinkRegion(x1[i], x2[j], x3[k])) {

                    rho = d->Vc[RHO][k][j][i];
                    vx1 = vx2 = vx3 = 0;
                    EXPAND(vx1 = d->Vc[VX1][k][j][i], ;
                            vx2 = d->Vc[VX2][k][j][i], ;
                                   vx3 = d->Vc[VX3][k][j][i];);
                    vs1 = VSPH1(x1[i], x2[j], x3[k], vx1, vx2, vx3);
                    vs1 = fabs(MAX(vs1, 0));
                    tr2 = d->Vc[TRC + 1][k][j][i];

                    rad = SPH1(x1[i], x2[j], x3[k]);
                    cell_vol = ElevateVolume(grid->dV[k][j][i]);
                    cell_area = pow(0.75 * cell_vol / CONST_PI, 2. / 3.) * CONST_PI;

                    wgt = tr2 * rho;
                    vrho = wgt * vs1;
                    mdot = vrho * cell_area;
                    pdot = vrho * vs1 * cell_area;
                    edot = 0.5 * vrho * vs1 * vs1 * cell_area;

                    wgt_a += wgt;
                    vrho_a += vrho;
                    mdot_a += mdot;
                    pdot_a += pdot;
                    edot_a += edot;
                }

            }

#ifdef PARALLEL
        MPI_Allreduce(&wgt_a,  &wgt,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&vrho_a, &vrho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mdot_a, &mdot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&pdot_a, &pdot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&edot_a, &edot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

        wgt = wgt_a;
        vrho = vrho_a;
        mdot = mdot_a;
        pdot = pdot_a;
        edot = edot_a;

#endif

        vrho /= wgt;

        ca.vrho_out_dom = vrho;
        ca.mdot_out_dom = mdot;
        ca.pdot_out_dom = pdot;
        ca.edot_out_dom = edot;

}


/* ********************************************************************* */
void WarmPhaseMass(Data *d, Grid *grid) {
/*!
 * Calculate warm phase mass.
 *
 *********************************************************************** */

    double rho, tr2, te;
    double mass_tr2, mass_tr2_a = 0;
    double mass_rtc, mass_rtc_a = 0;

    /* Region in phase to consider as warm */
    double rho_c = 1.e1;
    double te_c = CLOUD_TCRIT;


    /* Total mass in domain */

    double *x1, *x2, *x3;
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    double cell_vol;

    int i, j, k;
    DOM_LOOP(k, j, i) {

                if (! InSinkRegion(x1[i], x2[j], x3[k])) {

                    /* Primitives */
                    rho = d->Vc[RHO][k][j][i];
                    tr2 = d->Vc[TRC + 1][k][j][i];
                    te = d->Vc[PRS][k][j][i] / rho * MU_NORM * KELVIN;

                    cell_vol = ElevateVolume(grid->dV[k][j][i]);

                    /* Mass by tracer */
                    mass_tr2 = tr2 * rho * cell_vol;
                    mass_tr2_a += mass_tr2;

                    /* Mass by rho T cut */
                    if (te < te_c && rho < rho_c) {
                        mass_rtc = rho * cell_vol;
                        mass_rtc_a += mass_rtc;
                    }

                }

            }

#ifdef PARALLEL
        MPI_Allreduce(&mass_tr2_a, &mass_tr2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mass_rtc_a, &mass_rtc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

        mass_tr2 = mass_tr2_a;
        mass_rtc = mass_rtc_a;

#endif
    ca.mass_tr2_dom = mass_tr2;
    ca.mass_rtc_dom = mass_rtc;

}




/* ********************************************************************* */
void WarmPhaseStarFormationRate(Data *d, Grid *grid) {
/*!
 * Calculate total star-formation rate in clouds.
 * Remove mass from cells.
 *
 *********************************************************************** */

    double sfr, sfr_a = 0;

    double ***delta_m;
    delta_m = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);

    double *x1, *x2, *x3;
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    double cell_vol;

    int i, j, k;

    DOM_LOOP(k, j, i) delta_m[k][j][i] = 0;

    DOM_LOOP(k, j, i) {

                if (! InSinkRegion(x1[i], x2[j], x3[k])) {

                    sfr = StarFormationRateDensity(d, grid, i, j, k);

                    cell_vol = ElevateVolume(grid->dV[k][j][i]);

                    sfr_a += sfr * cell_vol;

                    delta_m[k][j][i] = sfr * g_dt;

                }

            }

    /* Remove mass */
    DOM_LOOP(k, j, i) {
                d->Uc[k][j][i][RHO] -= delta_m[k][j][i];
            }

    FreeArray3D((void *) delta_m);


#ifdef PARALLEL
    MPI_Allreduce(&sfr_a, &sfr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    sfr = sfr_a;

#endif

    ca.sfr = sfr;

}


/*********************************************************************** */
void WarmPhasePorosity(Data *d, Grid *grid) {
/*!
 * Calculate the porosity globally!
 *
 *********************************************************************** */

    /* Region in phase to consider as warm */
    double rho_c, te_c, tr2_c;
    WarmPhaseConditions(&rho_c, &te_c, &tr2_c);

    // TODO: Script this up.

    ca.porosity = 1.;

}



/*********************************************************************** */
void WarmPhaseConditions(double *rho_c, double *te_c, double *tr2_c) {
/*!
 * The critical values for clouds and SFR... not really thought through yet
 *
 *********************************************************************** */

    *rho_c = 1.e1;
    *te_c = CLOUD_TCRIT;
    *tr2_c = 0.98;

}


/*********************************************************************** */
void CloudOutput() {
/*!
 * Perform output of spherical accretion calculations
 *
 *********************************************************************** */

    if (prank == 0) {

        FILE *fp_acc;
        char fname[512];
        sprintf(fname, "warm_phase.dat");
        static double next_output = -1;

        next_output = OutputContextEnter(fname, &fp_acc, next_output, CLOUD_OUTPUT_RATE);

        double year = CONST_ly / CONST_c;

        /* Write data */
        if (g_time > next_output) {

            fprintf(fp_acc, "%12.6e  "
                            "%12.6e  %12.6e  %12.6e  %12.6e  "
                            "%12.6e  %12.6e  %12.6e  %12.6e  "
                            "%12.6e  %12.6e  %12.6e  %12.6e  "
                            "%12.6e  %12.6e  %12.6e  %12.6e  "
                            "%12.6e  %12.6e  "
                            "%12.6e  %12.6e  "
                            "\n",
                    g_time * vn.t_norm / year,
                    ca.vrho_out[0] * vn.v_norm / 1.e5,
                    ca.vrho_out[1] * vn.v_norm / 1.e5,
                    ca.vrho_out[2] * vn.v_norm / 1.e5,
                    ca.vrho_out_dom * vn.v_norm / 1.e5,
                    ca.mdot_out[0] * vn.mdot_norm / (CONST_Msun / year),
                    ca.mdot_out[1] * vn.mdot_norm / (CONST_Msun / year),
                    ca.mdot_out[2] * vn.mdot_norm / (CONST_Msun / year),
                    ca.mdot_out_dom * vn.mdot_norm / (CONST_Msun / year),
                    ca.pdot_out[0] * vn.power_norm / vn.v_norm,
                    ca.pdot_out[1] * vn.power_norm / vn.v_norm,
                    ca.pdot_out[2] * vn.power_norm / vn.v_norm,
                    ca.pdot_out_dom * vn.power_norm / vn.v_norm,
                    ca.edot_out[0] * vn.power_norm,
                    ca.edot_out[1] * vn.power_norm,
                    ca.edot_out[2] * vn.power_norm,
                    ca.edot_out_dom * vn.power_norm,
                    ca.mass_tr2_dom * vn.m_norm / CONST_Msun,
                    ca.mass_rtc_dom * vn.m_norm / CONST_Msun,
                    ca.porosity,
                    ca.sfr * vn.mdot_norm / (CONST_Msun / year)
            );

        }

        next_output = OutputContextExit(&fp_acc, next_output, CLOUD_OUTPUT_RATE);

    }

}


/*********************************************************************** */
void InputDataClouds(const Data *d, const Grid *grid) {
/*!
 * Routine that puts clouds into grid using PLUTO's input data mechanism.
 *
 *********************************************************************** */

// TODO: Hot phase has been initialized already everywhere. Check whether some of the
//       replacement routines are necessary in setting up the clouds.

    int i, j, k, nv, id, vol;
    double *x1, *x2, *x3;
    double halo_primitives[NVAR], cloud_primitives[NVAR];

    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    id = InputDataOpen("./input-rho.flt", "./grid_in.out", CUBE_ENDIANNESS, 0);
    TOT_LOOP(k,j,i) d->Vc[RHO][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
    InputDataClose(id);

#if CLOUD_VELOCITY

    id = InputDataOpen("./input-vx1.flt", "./grid_in.out", CUBE_ENDIANNESS, 0);
    TOT_LOOP(k,j,i) d->Vc[VX1][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
    InputDataClose(id);

    id = InputDataOpen("./input-vx2.flt", "./grid_in.out", CUBE_ENDIANNESS, 0);
    TOT_LOOP(k,j,i) d->Vc[VX2][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
    InputDataClose(id);

    id = InputDataOpen("./input-vx3.flt", "./grid_in.out", CUBE_ENDIANNESS, 0);
    TOT_LOOP(k,j,i) d->Vc[VX3][k][j][i] = InputDataInterpolate(id, x1[i], x2[j], x3[k]);
    InputDataClose(id);

    double ***v1 = d->Vc[VX1];
    double ***v2 = d->Vc[VX2];
    double ***v3 = d->Vc[VX3];
    InputDataCoordTransformVector(id, x1, x2, x3, v1, v2, v3);

#endif

    /* Do Clouds apodization */
    TOT_LOOP(k, j, i) {

                /* First get primitives array for hot halo */
                HotHaloPrimitives(halo_primitives, x1[i], x2[j], x3[k]);

                /* Fill cloud_primitives array */
                NVAR_LOOP(nv) cloud_primitives[nv] = d->Vc[nv][k][j][i];

                // TODO: Redo and check CLOUDS_MULTI

                /* If we're in the domain of the clouds cube */
                if (CloudPrimitives(cloud_primitives, x1[i], x2[j], x3[k])){
                    NVAR_LOOP(nv) d->Vc[nv][k][j][i] = cloud_primitives[nv];
                }
                    /* If not a cloud pixel then use hot halo primitives*/
                else{
                    NVAR_LOOP(nv) d->Vc[nv][k][j][i] = halo_primitives[nv];
                }
            }
}
