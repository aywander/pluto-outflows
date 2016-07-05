//
// Created by Alexander Y. Wagner on 4/6/16.
//

#include "pluto.h"
#include "pluto_usr.h"
#include "clouds.h"
#include "idealEOS.h"
#include "hot_halo.h"

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
void ReadFractalData()
/*!
 * This routine reads the Fractal cube data once.
 *
 **************************************************************** */
{

//---DM 22feb, 2015: update definition of get_var for cloud velocities---

    static int once01 = 0;

    if (!once01) {

#if CLOUD_VELOCITY != NONE
        int get_var[] = {RHO, ARG_EXPAND(VX1, VX2, VX3), -1};

#else

        int get_var[] = {RHO, -1};

#endif

        /* Read cloud data from external file */
        InputDataSet("./grid_in.out", get_var);
        InputDataRead("./input.flt", CUBE_ENDIANNESS);
        once01 = 1;
    }

}

/* ************************************************************** */
void GetFractalData(double *cloud, const double x1, const double x2, const double x3)
/*!
 * This routine returns the value of the fractal multiplication factor "fd"
 * at point (x1, x2, x3). The array "cloud" is filled in the process. For
 * the values specified in ReadFractalData.
 *
 **************************************************************** */
{

    double x, y, z;
    int nv, inv;

    /* Cloud data is in cartesian coordiantes
       InputDataInterpolate */
    x = CART1(x1, x2, x3);
    y = CART2(x1, x2, x3);
    z = CART3(x1, x2, x3);
    InputDataInterpolate(cloud, x, y, z);

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
    double r1, r2, y0, y1, y2, y3, frac, r_sph, r_cyl, phi, phi_cyl;
    double dens, sigma_g2, wrho, wrho_cgs, wrad_cgs, wtrb, ek;

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

    /* Must use gravitational potential */
#if BODY_FORCE != POTENTIAL
    fputs("Error: Must use BODY_FORCE = POTENTIAL\n", stderr);
    fputs("with CLOUD_DENSITY = CD_KEPLERIAN.\n", stderr);
    QUIT_PLUTO(1);
#endif

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
    wtrb = g_inputParam[PAR_WTRB]*ini_code[PAR_WTRB];
    sigma_g2 = wtrb*wtrb;

#endif

    /* Gravitational potential (cgs from table, then turn into code units) */
    r_sph = SPH1(x1, x2, x3);
    phi = InterpolationWrapper(gr_rad, gr_phi, gr_ndata, r_sph);

    /* Is the potential non-spherical? */
    ek = g_inputParam[PAR_WROT] * g_inputParam[PAR_WROT];

    if (ek > 0) {

        /* Now, the same for the cylindrical radius (in cgs) */
        r_cyl = CYL1(x1, x2, x3);
        phi_cyl = InterpolationWrapper(gr_rad, gr_phi, gr_ndata, r_cyl);
    }

    /* The profile */
    // TODO: we used to assume that phi(0,0) = 0 Eqn 15. Sutherland & Bicknell (2007), but this is no longer the case if a BH potential is included.
    dens = wrho * exp((-phi + phi_cyl * ek) / sigma_g2);


#elif CLOUD_DENSITY == CD_HOMOGENEOUS
    /* Homogeneous halo density, but with gravity */
    dens = g_inputParam[PAR_WRHO] * ini_code[PAR_WRHO];

#else
    /* Homogeneous halo, gravity off. */
    dens = g_inputParam[PAR_WRHO];

#endif

    /* Apodize! (code units) */
    cloud[RHO] = cloud[RHO] * dens;

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

    int is_cloud = 0;
    double fdratio;
    static int once01 = 0;

#if CLOUD_EXTRACT == CE_ELLIPSOID

    double rad, r_cyl, rz;
    double tanhfactor;
    double ellipse, ellipse_i, inner_circ;
    double wrot, wrad, osph, wsmf, incf;

    fdratio = cloud[RHO] / halo[RHO];

    /* The distances in physical space */
    rad   = SPH1(x1, x2, x3);
    r_cyl = CYL1(x1, x2, x3);
    rz    = CYL2(x1, x2, x3);

    /* The ellipse equation */
    wrot = g_inputParam[PAR_WROT] * ini_code[PAR_WROT];
    wrad = g_inputParam[PAR_WRAD] * ini_code[PAR_WRAD];
    ellipse = (r_cyl * r_cyl * (1 - wrot * wrot) + rz * rz) /
              (exp(-1) * wrad * wrad);

    /* Smoothing region scale, wsmf, and buffer factor around osph */
    wsmf = 0.2;
    incf = 1.2;

    /* The inner ellipse equation, beyond which smoothing region begins */
    ellipse_i = (r_cyl * r_cyl * (1 - wrot * wrot) + rz * rz) /
                ((1. - wsmf) * exp(-1) * wrad * wrad);

    /* Inner hemisphere to keep free */
    osph = g_inputParam[PAR_OSPH] * ini_code[PAR_OSPH];
    inner_circ = rad / (incf * osph);

    /* Exclude zone outside ellipsoid */
    if ( ellipse > 1. || inner_circ < 1.){
      fdratio = 1;
      is_cloud = 0;
    }

    /* The smoothing region */
    else if (!((ellipse_i < 1.) || (ellipse > 1.))){
      tanhfactor = tanh(tan(-CONST_PI * (sqrt(ellipse) - (1. - 0.5 * wsmf)) / wsmf));
      fdratio = 0.5 * (fdratio + 1) + 0.5 * tanhfactor * (fdratio - 1.);
      is_cloud = 1;
    }

    /* Inside the fractal sphere */
    else{
      is_cloud = 1;
    }

    if (!once01){
      print1("> Cloud extraction: CE_ELLIPSOID.\n\n");
      once01 = 1;
    }


#elif CLOUD_EXTRACT == CE_DENS
    /* A density based extraction method */

    double dens_lim = 0.0
    if (fd < dens_lim) {
      fdratio = 1.;
      is_cloud = 0;
    }
    fdratio = MAX(dens_lim, fdratio);
    is_cloud = 1;

#else
    is_cloud = 1;
    fdratio = cloud[RHO] / halo[RHO];

    if (!once01) {
        print1("> Cloud extraction: No special cloud extraction.\n\n");
        once01 = 1;
    }

#endif

    /* The ratio of cloud to halo density for a cell
     * shouldn't really ever be < 1. */
    fdratio = MAX(fdratio, 1.);

    /* Set to cloud density */
    cloud[RHO] = fdratio * halo[RHO];

    return is_cloud;
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

    /* The coulds at this point are in CGS units */
    double v1, v2, v3;

#if CLOUD_VELOCITY == CV_KEPLERIAN

    double ek;
    double r_cyl, r1, r2, frac;
    double y0, y1, y2, y3, phiderv_cyl;
    int il;
    double vpol1, vpol2, vpol3;
    double xpol1, xpol2, xpol3;

    /* The rotational parameter */
    ek = g_inputParam[PAR_WROT];

    EXPAND(v1 = cloud[VX1];,
           v2 = cloud[VX2];,
           v3 = cloud[VX3];);

    /* If ek == 0, then it's spherical) and velocities are just those read in.
     * If ek > 0., there is a keplerian component that needs to be added. */
    if (ek > 0.) {

        /* Convert coordinates and velocity vectors to cylindrical polars */
        xpol1 = POL1(x1, x2, x3);
        xpol2 = POL2(x1, x2, x3);
        xpol3 = POL3(x1, x2, x3);

        vpol1 = VPOL1(x1, x2, x3, v1, v2, v3);
        vpol2 = VPOL2(x1, x2, x3, v1, v2, v3);
        vpol3 = VPOL3(x1, x2, x3, v1, v2, v3);


        /* cycldrical radius */
        r_cyl = CYL1(x1, x2, x3);
        phiderv_cyl = InterpolationWrapper(gr_rad, gr_dphidr, gr_ndata, r_cyl);

        /* The angular velocity */
        vpol2 += g_inputParam[PAR_WROT] * sqrt(r_cyl * phiderv_cyl);

        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(v1 = VPOL_1(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3);,
               v2 = VPOL_2(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3);,
               v3 = VPOL_3(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3););


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

        /* Apply change to component paralell to flow axis (assumed to be in km/s)
         * Can't use FLOWAXIS macro though because we are in transformed coords. */
        SELECT(vcart1, vcart2, vcart3) += g_inputParam[PAR_WVPL] * ini_code[PAR_WVPL];

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

        /* Apply change to radial component (assumed to be in km/s) */
        vcart2 += g_inputParam[PAR_WVPP] * ini_code[PAR_WVPP];

        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(v1 = VCART_1(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3);,
               v2 = VCART_2(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3);,
               v3 = VCART_3(xcart1, xcart2, xcart3, vcart1, vcart2, vcart3););

    }

    if (fabs(g_inputParam[PAR_WVAN]) > 0) {

        /* Convert coordinates and velocity vectors to cylerical */

        xcyl1 = POL1(x1, x2, x3);
        xcyl2 = POL2(x1, x2, x3);
        xcyl3 = POL3(x1, x2, x3);

        vcyl1 = VPOL1(x1, x2, x3, v1, v2, v3);
        vcyl2 = VPOL2(x1, x2, x3, v1, v2, v3);
        vcyl3 = VPOL3(x1, x2, x3, v1, v2, v3);

        /* Apply change to radial component (assumed to be in km/s) */
        vcyl2 += g_inputParam[PAR_WVAN] * ini_code[PAR_WVAN];

        /* Convert velocity vectors back to the current coordinate system */
        EXPAND(v1 = VPOL_1(xcyl1, xcyl2, xcyl3, vcyl1, vcyl2, vcyl3);,
               v2 = VPOL_2(xcyl1, xcyl2, xcyl3, vcyl1, vcyl2, vcyl3);,
               v3 = VPOL_3(xcyl1, xcyl2, xcyl3, vcyl1, vcyl2, vcyl3););

    }


    /* Convert to code units here */
    EXPAND(v1 *= ini_code[PAR_WTRB];,
           v2 *= ini_code[PAR_WTRB];,
           v3 *= ini_code[PAR_WTRB];);
    EXPAND(cloud[VX1] = v1;, cloud[VX2] = v2;, cloud[VX3] = v3;);


#if USE_FOUR_VELOCITY == YES
    double vel = VMAG(x1, x2, x3, cloud[VX1], cloud[VX2], cloud[VX3]);
    double scrh = Vel2Lorentz(vel);
    EXPAND(cloud[VX1] *= scrh;, cloud[VX2] *= scrh;, cloud[VX3] *= scrh;);
#endif
}

/* ************************************************************** */
int CloudPrimitives(double *cloud,
                    const double x1, const double x2, const double x3)
/*!
 * This function returns 1 if cell is cloud material, 0 if not.
 *
 * The array *cloud is filled in the process that involves three steps:
 * 0) Reading cloud data once. ReadFractalData.
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

    /* Read in fractal data (once) and put into memory.
     * This needs to happen here because CloudCubePixel
     * requires cube data domain extents. */
    ReadFractalData();

    /* Get cloud pixel coordinates */
    if (CloudCubePixel(cube_pixel, x1, x2, x3)) {

        /* Get the fractal factor for this cell */
        GetFractalData(cloud, x1, x2, x3);

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
    if (is_cloud == 0) for (nv = 0; nv < NVAR; ++nv) cloud[nv] = halo[nv];

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