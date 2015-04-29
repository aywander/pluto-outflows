/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains functions to assist problem initialization in init.c.

  The init_tools.c file contains helper functions used mainly in init.c. This is to keep init.c relatively clean, with only functions that were originally defined. 

  \author AYW 
  \date   2012-12-07 14:41 JST
*/
/* ///////////////////////////////////////////////////////////////////// */


#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"
#include "rGravTable.h"
#include "rHotTable.h"
#include "idealEOS.h"
#include "abundances.h"
#include "interpolation.h"

/* Global struct and arrays for normalization */
VarNorm vn;
double ini_cgs[32];
double ini_code[32];

/* Global struct for nozzle */
Nozzle nz;


/* Functions */


void PrintInitData01(const double * out_primitives, 
                     const double * halo_primitives){

      /* Print some additional data during initialization */

      print1("\n");
      print1("> Conditions at (0, 0, 0):\n");
      print1("\n");
      print1("        %14s  %14s  %14s\n", 
          "Nozzle primitives","Hot halo prim","ratio");
      print1("  rho   %14g  %14g  %14g\n", 
          out_primitives[RHO], halo_primitives[RHO], 
          out_primitives[RHO]/halo_primitives[RHO]);
      print1("  pr    %14g  %14g  %14g\n", 
          out_primitives[PRS], halo_primitives[PRS],
          out_primitives[PRS]/halo_primitives[PRS]);
      print1("\n");
}


/* ************************************************************** */
void SetBaseNormalization() {
/*
 * Sets initializes VarNorm struct with derived normalizations
 * eint_norm is the (mass) specific internal energy [erg/g]
 *
 **************************************************************** */
  double primitives[NVAR];
  double mu;
  double year;

  /* Set primitives to get mu */
  primitives[RHO] = 1.;
  EXPAND(primitives[VX1] = 0.;, primitives[VX2] = 0.;, primitives[VX3] = 0.;);
  primitives[PRS] = 1.e-7;

  /* Set base normalization 
   * density: mean mass per particle
   * length: kpc
   * time: kyr
   * */

  vn.l_norm       = UNIT_LENGTH;
  vn.dens_norm    = UNIT_DENSITY;
  vn.v_norm       = UNIT_VELOCITY;
  vn.temp_norm    = KELVIN;

  /* Derived normalizations */
  vn.t_norm       =   vn.l_norm/vn.v_norm;
  vn.area_norm    =   vn.l_norm*vn.l_norm;
  vn.pres_norm    =   vn.dens_norm*vn.v_norm*vn.v_norm;
  vn.power_norm   =   vn.pres_norm*vn.v_norm*vn.area_norm;
  vn.eflux_norm   =   vn.pres_norm*vn.v_norm;
  vn.eint_norm    =   vn.pres_norm/vn.dens_norm;
  vn.mdot_norm    =   vn.dens_norm*pow(vn.l_norm,3)/vn.t_norm;
  vn.newton_norm  =   vn.t_norm*vn.t_norm*vn.dens_norm;

  print1("> Base normalization struct initialized.\n\n");

  return;
}


/* ************************************************************** */
void SetIniNormalization() {
/*
 * Sets noramlizations for ini input variables.
 * ini_cgs converts ini parameters to cgs upon multiplication,
 * whereas ini_code converts ini parameters to code units.
 *
 **************************************************************** */

  double year, degrad;
  year = CONST_ly/CONST_c;
  degrad = CONST_PI/180.;

  ini_cgs[PAR_OPOW] = 1.;                                  ini_code[PAR_OPOW] = ini_cgs[PAR_OPOW]/vn.power_norm;
  ini_cgs[PAR_OSPD] = NOZZLE_SELECT(1.,CONST_c);           ini_code[PAR_OSPD] = ini_cgs[PAR_OSPD]/NOZZLE_SELECT(1.,vn.v_norm);
  ini_cgs[PAR_OMDT] = NOZZLE_SELECT(1.,CONST_Msun/year);   ini_code[PAR_OMDT] = ini_cgs[PAR_OMDT]/NOZZLE_SELECT(1.,vn.mdot_norm);
  ini_cgs[PAR_OANG] = degrad;                              ini_code[PAR_OANG] = ini_cgs[PAR_OANG];
  ini_cgs[PAR_ORAD] = vn.l_norm;                           ini_code[PAR_ORAD] = ini_cgs[PAR_ORAD]/vn.l_norm;
  ini_cgs[PAR_ODIR] = degrad;                              ini_code[PAR_ODIR] = ini_cgs[PAR_ODIR];
  ini_cgs[PAR_OOMG] = degrad/(1.e6*year);                  ini_code[PAR_OOMG] = ini_cgs[PAR_OOMG]*vn.t_norm;
  ini_cgs[PAR_OPHI] = degrad;                              ini_code[PAR_OPHI] = ini_cgs[PAR_OPHI];
  ini_cgs[PAR_ODBH] = vn.l_norm;                           ini_code[PAR_ODBH] = ini_cgs[PAR_ODBH]/vn.l_norm;
  ini_cgs[PAR_OSPH] = vn.l_norm;                           ini_code[PAR_OSPH] = ini_cgs[PAR_OSPH]/vn.l_norm;
  ini_cgs[PAR_HRHO] = vn.dens_norm;                        ini_code[PAR_HRHO] = ini_cgs[PAR_HRHO]/vn.dens_norm;
  ini_cgs[PAR_HTMP] = 1;                                   ini_code[PAR_HTMP] = ini_cgs[PAR_HTMP]/vn.temp_norm;
  ini_cgs[PAR_HVBG] = DIMDIM3(1,1.e5);                     ini_code[PAR_HVBG] = ini_cgs[PAR_HVBG]/DIMDIM3(1,vn.v_norm);
  ini_cgs[PAR_HRAD] = vn.l_norm;                           ini_code[PAR_HRAD] = ini_cgs[PAR_HRAD]/vn.l_norm;
  ini_cgs[PAR_WRHO] = vn.dens_norm;                        ini_code[PAR_WRHO] = ini_cgs[PAR_WRHO]/vn.dens_norm;
  ini_cgs[PAR_WTRB] = 1.e5;                                ini_code[PAR_WTRB] = ini_cgs[PAR_WTRB]/vn.v_norm;
  ini_cgs[PAR_WRAD] = vn.l_norm;                           ini_code[PAR_WRAD] = ini_cgs[PAR_WRAD]/vn.l_norm;
  ini_cgs[PAR_WROT] = 1.;                                  ini_code[PAR_WROT] = ini_cgs[PAR_WROT];
  ini_cgs[PAR_WSMF] = 1.;                                  ini_code[PAR_WSMF] = ini_cgs[PAR_WSMF];
  ini_cgs[PAR_LEV1] = 1;                                   ini_code[PAR_LEV1] = ini_cgs[PAR_LEV1];

  print1("> Ini parameter normalization array initialized.\n\n");

  return;
}


/* ************************************************ */
void SetNozzleConeGeometry(){
/*
 * Set outflow geometry struct with parameters of cone 
 *
 * NOTE:
 *  - JetPrimitives and UfoPrimities only differ in their
 *    normalizations and how the paramters are set.
 *
 ************************************************** */

  double small_angle = 1.e-12;
  double large = 1.e30;

  /* Quantitites from input parameters */
  nz.ang = g_inputParam[PAR_OANG]*ini_code[PAR_OANG];
  nz.rad = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  nz.dir = g_inputParam[PAR_ODIR]*ini_code[PAR_ODIR];
  nz.dbh = g_inputParam[PAR_ODBH]*ini_code[PAR_ODBH];
  nz.omg = g_inputParam[PAR_OOMG]*ini_code[PAR_OOMG];
  nz.phi = g_inputParam[PAR_OPHI]*ini_code[PAR_OPHI];
  nz.sph = g_inputParam[PAR_OSPH]*ini_code[PAR_OSPH];
  
  /* Derived quantitites */

#if INTERNAL_BOUNDARY == YES
  nz.cbh = sqrt(nz.sph*nz.sph - nz.rad*nz.rad);
  nz.orig = nz.cbh;
#else
  nz.cbh = (nz.orig - nz.dbh)/cos(nz.dir) + nz.rad*tan(nz.dir);
  nz.orig = g_domBeg[FLOWAXIS(IDIR,JDIR,KDIR)];
#endif
    
  if (nz.ang > small_angle) { 
    nz.isfan = 1;
    nz.area = 2.*CONST_PI*(1. - cos(nz.ang))*pow(nz.rad/sin(nz.ang),2);

    /* cone apex is only valid when cone is aligned with flow axis */
    nz.cone_height = nz.rad/tan(nz.ang);
    nz.cone_apex = - (nz.cone_height - nz.orig);
  }
  else { 
    nz.isfan = 0; 
    nz.area = CONST_PI*nz.rad*nz.rad;
    nz.cone_height = nz.cbh;
    nz.cone_apex = -nz.cbh;
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
    /* This check is for avoiding the cone of the nozzle to be burried in the 
     outflow boundary for half-galaxy simulations. */
    if (acos(1 - ((nz.sph - nz.cbh)*(nz.sph - nz.cbh) + nz.rad*nz.rad)/
        (2*nz.sph*nz.sph)) + nz.dir > CONST_PI/2) {
        print1("Error: OSPH is too small. It must be at least...");
        QUIT_PLUTO(1);
        /* TODO: Calculate minimum OSPH */
    }
#endif

}


/* ************************************************ */
void OutflowPrimitives(double* out_primitives, 
    const double x1, const double x2, const double x3) {
/*
 * Runs the relevant primitives function for nozzle flow.
 *
 * NOTE:
 *  - JetPrimitives and UfoPrimities only differ in their
 *    normalizations and how the paramters are set.
 *
 ************************************************** */

  /* Get primitives array depending on outflow type */
  NOZZLE_SELECT(JetPrimitives(out_primitives, x1, x2, x3),
                UfoPrimitives(out_primitives, x1, x2, x3));

}


/* ************************************************ */
void OutflowVelocity(double * out_primitives, double speed,
    const double x1, const double x2, const double x3){
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

  double vx1, vx2, vx3;
  double cvx1, cvx2, cvx3;
  int mirror_side = 0;

  if (nz.isfan) {

#if INTERNAL_BOUNDARY == YES
      /* If we're in the counter-nozzle region use mirror symmetric value of cx1p
       NOTE, since we've rotated to flow-axis, the nozzle is cylindrically symmetric
       and we don't have to use rotational symmetry. */
      
      if (FLOWAXIS(cx1p, cx2p, cx3p) < 0){
          mirror_side = 1;
          FLOWAXIS(cx1p, cx2p, cx3p) *= -1;
      }
#endif
      
    FLOWAXIS(cx1p, cx2p, cx3p) -= nz.cone_apex;

    double sx1, sx2, sx3;
    sx1 = SPH1(cx1p, cx2p, cx3p);
    sx2 = SPH2(cx1p, cx2p, cx3p);
    sx3 = SPH3(cx1p, cx2p, cx3p);

    double cur_vx1, cur_vx2, cur_vx3;
    cur_vx1 = VSPH_1(sx1, sx2, sx3, speed, 0, 0);
    cur_vx2 = VSPH_2(sx1, sx2, sx3, speed, 0, 0);
    cur_vx3 = VSPH_3(sx1, sx2, sx3, speed, 0, 0);
      
    cvx1 = VCART1(cx1p, cx2p, cx3p, cur_vx1, cur_vx2, cur_vx3);
    cvx2 = VCART2(cx1p, cx2p, cx3p, cur_vx1, cur_vx2, cur_vx3);
    cvx3 = VCART3(cx1p, cx2p, cx3p, cur_vx1, cur_vx2, cur_vx3);

    FLOWAXIS(cx1p, cx2p, cx3p) += nz.cone_apex;
      
#if INTERNAL_BOUNDARY == YES
      /* Create mirror-symmetrically opposite velocity vector
       and mirror back cell position */
      
      if (mirror_side){
          FLOWAXIS(cx1p, cx2p, cx3p) *= -1;
          FLOWAXIS(cvx1, cvx2, cvx3) *= -1;
      }

#endif
      
  }
  else { // if nozzle is not fan

    cvx1 = 0; cvx2 = 0; cvx3 = 0;
    FLOWAXIS(cvx1, cvx2, cvx3) = speed;
      
#if INTERNAL_BOUNDARY == YES
      /* Create mirror-symmetrically opposite velocity vector */
      
      if (FLOWAXIS(cx1p, cx2p, cx3p) < 0){
          FLOWAXIS(cvx1, cvx2, cvx3) *= -1;

      }
      
#endif
  } // if nozzle is fan

  /* Rotate vector back */
  double cvx1p, cvx2p, cvx3p; 
  RotateNozzle2Grid(cvx1, cvx2, cvx3, &cvx1p, &cvx2p, &cvx3p);

  EXPAND(vx1 = VCART_1(cx1, cx2, cx3, cvx1p, cvx2p, cvx3p);,
         vx2 = VCART_2(cx1, cx2, cx3, cvx1p, cvx2p, cvx3p);,
         vx3 = VCART_3(cx1, cx2, cx3, cvx1p, cvx2p, cvx3p););

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
void JetPrimitives(double * jet_primitives, 
    const double x1, const double x2, const double x3) {
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
  power   = g_inputParam[PAR_OPOW]*ini_code[PAR_OPOW];
  chi     = g_inputParam[PAR_OMDT]*ini_code[PAR_OMDT];
  lorentz = g_inputParam[PAR_OSPD]*ini_code[PAR_OSPD];

  /* Some derived quantities */
  double vel, dens, pres, gmm1;
  vel = Lorentz2Vel(lorentz);
  gmm1 = g_gamma/(g_gamma - 1.);
  pres = power/(gmm1*lorentz*lorentz*vel*nz.area*
      (1. + (lorentz - 1.)/lorentz*chi));
  dens = chi*gmm1*pres;

  /* Set primitives array */
  jet_primitives[RHO] = dens;
  jet_primitives[PRS] = pres;
  jet_primitives[TRC] = 1.0;
#if CLOUDS
  jet_primitives[TRC+1] = 0.0;
#endif
  OutflowVelocity(jet_primitives, vel, x1, x2, x3);

  return;
}




/* ************************************************ */
void UfoPrimitives(double * ufo_primitives, 
    const double x1, const double x2, const double x3) {
/*
 * Returns the array of primitives for ufo parameters
 * power, angle, speed, mdot, and radius and width 
 *
 ************************************************** */

  double power, speed, mdot;

  /* Input parameters, normalized in code units*/
  power = g_inputParam[PAR_OPOW]*ini_code[PAR_OPOW];
  speed = g_inputParam[PAR_OSPD]*ini_code[PAR_OSPD];
  mdot  = g_inputParam[PAR_OMDT]*ini_code[PAR_OMDT];

  /* Derived quantities */
  double dens, pres;
  pres = (power - 0.5*mdot*speed*speed)*(g_gamma-1)/(g_gamma*nz.area*speed);
  dens = mdot/(nz.area*speed);


  /* Set primitives array */
  ufo_primitives[RHO] = dens;
  ufo_primitives[PRS] = pres;
  ufo_primitives[TRC] = 1.0;
#if CLOUDS
  ufo_primitives[TRC+1] = 0.0;
#endif

  OutflowVelocity(ufo_primitives, speed, x1, x2, x3);

  return;
}



/* ************************************************************** */
void HotHaloPrimitives(double * halo, 
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
  if (gr_rad == NULL){
    readGravTable();
  }
#endif

  /* Consider different distributions */

  /* Hernquist potential (hydrostatic)*/
#if HOT_DISTR == HOT_HERNQUIST
  rho0 = g_inputParam[PAR_HRHO]*ini_code[PAR_HRHO];
  r = SPH1(x1, x2, x3);
  a = g_inputParam[PAR_HRAD]*ini_code[PAR_HRAD];
  rs = r/a;
  halo[RHO] = rho0/(rs*pow((1 + rs), 3));
  halo[PRS] = -2*CONST_PI*CONST_G/vn.newton_norm*rho0*rho0*a*a*
    (pow(1.+rs, -4)/12.*(25. + 52.*rs + 42.*rs*rs + 12.*rs*rs*rs) -
    log(1. + rs) + log(rs));


  /* Any spherical potential read from a table */
#elif HOT_DISTR == HOT_EXT_DATA
  
  if (hot_rad == NULL){
    readHotTable();
  }

  /* Radius of current cell */
  r = SPH1(x1, x2, x3);
  //r = sqrt(D_EXPAND(x1*x1, +x2*x2, +x3*x3));

  /* Find cell left index to interpolate at */
  il = hunter(hot_rad, hot_ndata, r);

  /* Linear fractional location of r in cell */
  r1 = hot_rad[il];
  r2 = hot_rad[il+1];
  frac = (r - r1)/(r2 - r1);

  /* Density interpolation */
  y0 = hot_rho[il-1];
  y1 = hot_rho[il];
  y2 = hot_rho[il+1];
  y3 = hot_rho[il+2];
  halo[RHO] = CubicCatmullRomInterpolate(y0, y1, y2, y3, frac);

  /* Pressure interpolation */
  y0 = hot_prs[il-1];
  y1 = hot_prs[il];
  y2 = hot_prs[il+1];
  y3 = hot_prs[il+2];
  halo[PRS] = CubicCatmullRomInterpolate(y0, y1, y2, y3, frac);

  /* Homogeneous denisty */
#elif HOT_DISTR == HOT_FLAT
  /* Flat thermodynamic profile but gravity is on 
   * The potential, thus, is a parabola, and not flat.
   * The initial conditions are the same as for the case
   * without gravity (else clause).*/

  halo[RHO] = g_inputParam[PAR_HRHO]*ini_code[PAR_HRHO];
  halo[PRS] = PresIdealEOS(halo[RHO], g_inputParam[PAR_HTMP]*ini_code[PAR_HTMP], MU_NORM);


  /* Default is flat */
#else

  halo[RHO] = g_inputParam[PAR_HRHO]*ini_code[PAR_HRHO];
  halo[PRS] = PresIdealEOS(halo[RHO], g_inputParam[PAR_HTMP]*ini_code[PAR_HTMP], MU_NORM);

#endif


  /* Velocities. */

  /* TODO: At some point we may want to generalize the 
   * background velocity assignment to three components.
   * This is mainly for simulations of cluster radio galaxies. 
   * */
  EXPAND(halo[VX1] = 0;,halo[VX2] = 0;, halo[VX3] = 0;);
  FLOWAXIS(halo[VX3], halo[VX1], halo[VX2]) = g_inputParam[PAR_HVBG]*ini_code[PAR_HVBG];
#if USE_FOUR_VELOCITY == YES
  vel = VMAG(x1, x2, x3, halo[VX1], halo[VX2], halo[VX3]);
  scrh = Vel2Lorentz(vel);
  EXPAND(halo[VX1] *= scrh;, halo[VX2] *= scrh;, halo[VX3] *= scrh;);
#endif


  /* Tracers */
  halo[TRC] = 0.0;
#if CLOUDS
  halo[TRC+1] = 0.0;
#endif


   /* A message for having initialized halo with potential*/
  if (!once01){
    print1("> Initializing hot halo distribution of type: %d\n\n", HOT_DISTR);
    once01 = 1;
  }

  return;

}


#if CLOUDS
/* ************************************************************** */
int CloudCubePixel(int * el, const double x1, 
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
  D_EXPAND(xfrac = (x - g_idBoxBeg[0])/csz1;,
           yfrac = (y - g_idBoxBeg[1])/csz2;,
           zfrac = (z - g_idBoxBeg[2])/csz3;);

  if (!(D_EXPAND(xfrac > 1., || yfrac > 1., || zfrac > 1.) || 
        D_EXPAND(xfrac < 0., || yfrac < 0., || zfrac < 0.))){

    /* Fill cube pixel element array */ 
    D_EXPAND(el[0] = (int) (xfrac*g_idnx1);,
             el[1] = (int) (yfrac*g_idnx2);,
             el[2] = (int) (zfrac*g_idnx3););
    return 1;
  }

  else{
    return 0;
  }

}
#endif


#if CLOUDS
/* ************************************************************** */
void ReadFractalData()
/*!
 * This routine reads the Fractal cube data once. 
 *
 **************************************************************** */
{

#if CLOUD_VELOCITY && ((CLOUD_VEL_DISTR == CV_KEPLERIAN_FRAC) || (CLOUD_VEL_DISTR == CV_VIRIAL_FRAC))

#define ID_MAX_NVAR (1 + COMPONENTS + 1)

  int get_var[ID_MAX_NVAR];
  get_var[0] = RHO;
  EXPAND(get_var[1] = VX1;,
         get_var[2] = VX2;,
         get_var[3] = VX3;);
  get_var[ID_MAX_NVAR-1] = -1;

#else

  int get_var[] = {RHO, -1};

#endif

  static int once01 = 0;

  /* Read cloud data from external file */
  if (!once01) {
    InputDataSet("./grid_in.out", get_var);
    InputDataRead("./input.flt", CUBE_ENDIANNESS);
    once01 = 1;
  }

}
#endif


#if CLOUDS
/* ************************************************************** */
void GetFractalData(double* cloud, const double x1, const double x2, const double x3)
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

  /* Assume input cloud velocities to be in km/s.
   * Convert to code units here */
  EXPAND(cloud[VX1] *= 1.e5 / vn.v_norm;, 
         cloud[VX2] *= 1.e5 / vn.v_norm;, 
         cloud[VX3] *= 1.e5 / vn.v_norm;)

}
#endif


#if CLOUDS
/* ************************************************************** */
void CloudApodize(double* cloud, const double x1, const double x2, const double x3)
/*!
 * Multiply cloud fractal factor (currently in *cloud) with a 
 * desired mean density profile. Currently CLOUD_DISTR = 
 *     - CD_HERNQUIST
 *     - CD_TURB_ISOTH_HYDROSTATIC
 *     - CD_TURB_KEPLERIAN_DISC
 *     - CD_FLAT
 * The default is also CD_FLAT. This is the apodization step. All 
 * cells are still cloud material.
 *
 * For the cloud velocity, currently CLOUD_VEL_DISTR
 *     - CV_UNIFORM
 *     - CV_RADIAL
 *     - CV_ROT_KEPLER_FRAC
 *     - CV_VIRIAL_FRAC
 * Each mode is associated with a different meaning for CV_VALUE
 * The default is CV_UNIFORM with CV_VALUE = 0
 *
 **************************************************************** */
{

  double z, r, r_cyl, a, rs; 
  int il;
  double r1, r2, y0, y1, y2, y3, frac, phi, phi_cyl;
  double mu, dens;
  double sigma_g2, csound2;
  double wtemp_av, wrad, wrad_cgs, wrho, wrho_cgs;
  double wtrb, wtrb_cgs, wrot;

  /* The following are only some 
   * profiles for the warm phase that produce
   * reasonable ratios of mean warm phase density and 
   * hot phase density in large parts of the domain. */
#if CLOUD_DISTR == CD_HERNQUIST
  r = SPH1(x1, x2, x3);
  a = g_inputParam[PAR_WRAD]*ini_code[PAR_WRAD];
  rs = r/a;
  wrho = g_inputParam[PAR_WRHO]*ini_code[PAR_WRHO];
  dens = wrho/(rs*pow((1 + rs), 3));

#elif CLOUD_DISTR == CD_TURB_ISOTH_HYDROSTATIC 

  /* Must use gravitational potential */
#if BODY_FORCE != POTENTIAL
  fputs("Error: Must use BODY_FORCE = POTENTIAL.\n", stderr); 
  fputs("with CLOUD_DISTR = CD_TURB_ISOTH_HYDROSTATIC.\n", stderr); 
  QUIT_PLUTO(1);
#endif

  /* The dense gas mean density input parameter */
  wrho_cgs = g_inputParam[PAR_WRHO]*ini_cgs[PAR_WRHO];
  wrho = g_inputParam[PAR_WRHO]*ini_code[PAR_WRHO];

  /* Relationship between warm phase temperature, the turbulent 
   * velocity dispersion (wtrb), the total velocity dispersion^2 
   * (sigma_g2), and the scale height of the warm phase (wrad). 
   * Everything is in cgs in this block */
#if CLOUD_SCALE == CS_WTRB

#if MU_CALC != MU_CONST
  /* Note, in case a variable mu is used, we cannot 
   * calculate a variable mu because cloud is * imcomplete. */
  fputs("Error: Must use MU_CALC = MU_CONST\n", stderr); 
  fputs("with CLOUD_SCALE = CS_WTRB as cloud\n", stderr); 
  fputs("primitives are not yet defined.\n", stderr); 
  QUIT_PLUTO(1);
#endif

  wtrb_cgs = g_inputParam[PAR_WTRB]*ini_cgs[PAR_WTRB];
  mu = MeanMolecularWeight(cloud);
  // Problem here, because we have neither cloud[PRS] yet nor 
  // a reference pressure in this routine
  wtemp_av = TempIdealEOS(wrho, cloud[PRS], mu)*vn.temp_norm;
  csound2  = CONST_kB*wtemp_av/(mu*CONST_amu);
  sigma_g2 = wtrb*wtrb + csound2;
  wrad_cgs = sqrt(9.*sigma_g2/(4.*CONST_PI*CONST_G*wrho_cgs));

#else
  wrad_cgs = g_inputParam[PAR_WRAD]*ini_cgs[PAR_WRAD];
  sigma_g2 = 4*CONST_PI*CONST_G*wrho_cgs*pow(wrad_cgs, 2)/9.;
  mu = MeanMolecularWeight(cloud);
  // Problem here, because we have neither cloud[PRS] yet nor
  // a reference pressure in this routine
  wtemp_av = TempIdealEOS(wrho, cloud[PRS], mu)*vn.temp_norm;
  csound2  = CONST_kB*wtemp_av/(mu*CONST_amu);
  wtrb_cgs = sqrt(sigma_g2 - csound2);
  // But ok, because we have sigma_g2 already.

#endif

  /* Code units for total dens gas velocity dispersion */
  sigma_g2 /= vn.v_norm*vn.v_norm;

  /* Gravitational potential */
  r = SPH1(x1, x2, x3);

  /* Find cell left index to interpolate at */
  il = hunter(gr_rad, gr_ndata, r);

  /* Linear fractional location of r in cell */
  r1 = gr_rad[il];
  r2 = gr_rad[il+1];
  frac = (r - r1)/(r2 - r1);

  /* Do potential interpolation */
  y0 = gr_phi[il-1];
  y1 = gr_phi[il];
  y2 = gr_phi[il+1];
  y3 = gr_phi[il+2];
  phi = CubicCatmullRomInterpolate(y0, y1, y2, y3, frac);

  /* The profile */
  dens = wrho*exp(-phi/sigma_g2);


#elif CLOUD_DISTR == CD_TURB_KEPLERIAN_DISC

  /* Must use gravitational potential */
#if BODY_FORCE != POTENTIAL
  fputs("Error: Must use BODY_FORCE = POTENTIAL\n", stderr); 
  fputs("with CLOUD_DISTR = CD_TURB_KEPLERIAN_DISC.\n", stderr); 
  QUIT_PLUTO(1);
#endif


  /* The dense gas mean density input parameter */
  wrho_cgs = g_inputParam[PAR_WRHO]*ini_cgs[PAR_WRHO];
  wrho = g_inputParam[PAR_WRHO]*ini_code[PAR_WRHO];

  /* Relationship between warm phase temperature, the turbulent 
   * velocity dispersion (wtrb), the total velocity dispersion^2 
   * (sigma_g2), and the scale height of the warm phase (wrad). 
   * Everything is in cgs in this block */
#if CLOUD_SCALE == CS_WTRB

#if MU_CALC != MU_CONST
  /* Note, in case a variable mu is used, we cannot 
   * calculate a variable mu because cloud is * imcomplete. */
  fputs("Error: Must use MU_CALC = MU_CONST\n", stderr); 
  fputs("with CLOUD_SCALE = CS_WTRB as cloud\n", stderr); 
  fputs("primitives are not yet defined.\n", stderr); 
  QUIT_PLUTO(1);
#endif

  wtrb_cgs = g_inputParam[PAR_WTRB]*ini_cgs[PAR_WTRB];
  mu = MeanMolecularWeight(cloud);
  // Problem here, because we have neither cloud[PRS] yet nor 
  // a reference pressure in this routine
  wtemp_av = TempIdealEOS(wrho, cloud[PRS], mu)*vn.temp_norm;
  csound2  = CONST_kB*wtemp_av/(mu*CONST_amu);
  sigma_g2 = wtrb*wtrb + csound2;
  wrad_cgs = sqrt(9.*sigma_g2/(4.*CONST_PI*CONST_G*wrho_cgs));

#else
  wrad_cgs = g_inputParam[PAR_WRAD]*ini_cgs[PAR_WRAD];
  sigma_g2 = 4*CONST_PI*CONST_G*wrho_cgs*pow(wrad_cgs, 2)/9.;
  mu = MeanMolecularWeight(cloud);
  // Problem here, because we have neither cloud[PRS] yet nor
  // a reference pressure in this routine
  wtemp_av = TempIdealEOS(wrho, cloud[PRS], mu)*vn.temp_norm;
  csound2  = CONST_kB*wtemp_av/(mu*CONST_amu);
  wtrb_cgs = sqrt(sigma_g2 - csound2);
  // But ok, because we have sigma_g2 already.

#endif

  /* Code units for total dense gas velocity dispersion */
  sigma_g2 /= vn.v_norm*vn.v_norm;

  /* Gravitational potential*/
  r     = SPH1(x1, x2, x3);
  r_cyl = CYL1(x1, x2, x3);
  z     = CYL2(x1, x2, x3);

  /* Find cell left index to interpolate at */
  il = hunter(gr_rad, gr_ndata, r);

  /* Linear fractional location of r in cell */
  r1 = gr_rad[il];
  r2 = gr_rad[il+1];
  frac = (r - r1)/(r2 - r1);

  /* Do potential interpolation */
  y0 = gr_phi[il-1];
  y1 = gr_phi[il];
  y2 = gr_phi[il+1];
  y3 = gr_phi[il+2];
  phi = CubicCatmullRomInterpolate(y0, y1, y2, y3, frac);

  /* Now, the same for the cylindrical radius */
  il = hunter(gr_rad, gr_ndata, r_cyl);
  r1 = gr_rad[il];
  r2 = gr_rad[il+1];
  frac = (r_cyl - r1)/(r2 - r1);
  y0 = gr_phi[il-1]; 
  y1 = gr_phi[il]; 
  y2 = gr_phi[il+1]; 
  y3 = gr_phi[il+2];
  phi_cyl = CubicCatmullRomInterpolate(y0, y1, y2, y3, frac);

  /* The profile */
  wrot = g_inputParam[PAR_WROT]*g_inputParam[PAR_WROT];
  dens = wrho*exp((-phi + phi_cyl*wrot)/sigma_g2);


#elif CLOUD_DISTR == CD_FLAT
  /* Homogeneous halo density, but with gravity */
  dens = g_inputParam[PAR_WRHO];

#else
  /* Homogeneous halo, gravity off. */
  dens = g_inputParam[PAR_WRHO];

#endif

  /* Apodize! */
  cloud[RHO] = cloud[RHO]*dens;

}
#endif



#if CLOUDS
/* ************************************************************** */
int CloudExtract(double* cloud,
                 const double* halo, 
                 const int* pixel, 
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
  double fdratio, rad, rrad, rr_cyl, z, rz;
  double dens_lim;
  double rw, tanhfactor;
  double ellipse, ellipse_i, inner_circ; 
  double wrot, wrad, hemf, urad, uthk, wsmf;
  double rcx, rcy, rcz;
  int crs;
  int is_cloud = 0;
  static int once01 = 0;

  fdratio = cloud[RHO]/halo[RHO];

#if CLOUD_EXTRACT == CE_SPHERE_RIM
  /* Radius of cube pixel from center of cube 
   * rad               the radius in pixel space 
   * rw                rim width in cells which is the smoothing region
   * 0.5*(crs - rw)    mid point between sphere within fractal cube
   *                   and sphere radius - rw. 
   * The origin of the pixel coordinates (here pixel[]) is at the corner 
   * with the smallest values of x1, x2, and x3.
   */

  /* Clouds cube is always assumed to be a data file in cartesian */
  D_EXPAND(rcx = pixel[0] - 0.5*g_idnx1;,
           rcy = pixel[1] - 0.5*g_idnx2;,
           rcz = pixel[2] - 0.5*g_idnx3;);
  rad = SPH1(rcx, rcy, rcz);

  /* The number of pixels of the smoothing region 
   * and the inner exclusion zone fraction is hard-coded */
  rw = 12; 
  hemf = 1.2;

  /* The radius in physical space */
  rrad = SPH1(x1, x2, x3)

  /* Inner hemisphere to keep free */
  urad = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  //uthk = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
  inner_circ = rrad/(hemf*(urad + uthk));

  /* Exclude zone outside interior sphere */
  wrad = g_inputParam[PAR_WRAD]*ini_code[PAR_WRAD];
  if (rad > 0.5*crs || rrad > wrad || inner_circ < 1.){
    fdratio = 1;
    is_cloud = 0;
  }
  /* The smoothing region */
  else if (!((rad < 0.5*crs - rw )|| (rad > 0.5*crs))){
    tanhfactor = tanh(tan(CONST_PI*(-rad + 0.5*(crs - rw))/((double) rw)));
    fdratio = 0.5*(fdratio + 1) + 0.5*tanhfactor*(fdratio - 1.);
    is_cloud = 1;
  }
  /* Inside the fractal sphere */
  else{
    is_cloud = 1;
  }

  if (!once01){
    print1("> Cloud extraction: CE_SPHERE_RIM.\n\n");
    once01 = 1;
  }


#elif CLOUD_EXTRACT == CE_HEMISPHERE_RIM
  /* Radius of cube pixel from center of cube at face x=0
   * rad               the radius in pixel space
   * rw                rim width in cells
   * 0.5*(crs - rw)    mid point between sphere within fractal cube
   *                   and sphere radius - rw. 
   * The origin of the pixel coordinates (here pixel[]) is at the corner 
   * with the smallest values of x1, x2, and x3.
   */
  D_EXPAND(rcx = pixel[0] - 0.5*g_idnx1;,
           rcy = pixel[1] - 0.5*g_idnx2;,
           rcz = pixel[2];);
  rad = SPH1(rcx, rcy, rcz);

  /* The number of pixels of the smoothing region 
   * and the inner exclusion zone fraction is hard-coded */
  rw = 12;
  hemf = 1.2;

  /* The radius in physical space */
  rrad = SPH1(x1, x2, x3);

  /* Inner hemisphere to keep free */
  urad = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  //uthk = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
  inner_circ = rrad/(hemf*(urad + 0.5*uthk));

  /* Exclude zone outside interior sphere */
  wrad = g_inputParam[PAR_WRAD]*ini_code[PAR_WRAD];
  if (rad > 0.5*crs || rrad > wrad || inner_circ < 1.){
    fdratio = 1;
    is_cloud = 0;
  }

  /* The smoothing region */
  else if (!((rad < 0.5*crs - rw) || (rad > 0.5*crs))){
    tanhfactor = tanh(tan(CONST_PI*(-rad + 0.5*(crs - rw))/((double) rw)));
    fdratio = 0.5*(fdratio + 1) + 0.5*tanhfactor*(fdratio - 1.);
    is_cloud = 1;
  }

  /* Inside the fractal sphere */
  else{
    is_cloud = 1;
  }

  if (!once01){
    print1("> Cloud extraction: CE_HEMISPHERE_RIM.\n\n");
    once01 = 1;
  }


#elif CLOUD_EXTRACT == CE_ELLIPSOID
  /* Radius of cube pixel from center of cube at face x=0
   * z                 the x1 distance
   * rad               the cylindrical distance
   * rw                rim width in cells
   * 0.5*(crs - rw)    mid point between sphere within fractal cube
   *                   and sphere radius - rw. 
   * The origin of the pixel coordinates (here pixel[]) is at the corner 
   * with the smallest values of x1, x2, and x3.
   */
  D_EXPAND(rcx = pixel[0] - 0.5*g_idnx1;,
           rcy = pixel[1] - 0.5*g_idnx2;,
           rcz = pixel[2];);
  rad = CYL1(rcx, rcy, rcz);
  z   = CYL2(rcx, rcy, rcz);

  /* The number of pixels of the smoothing region 
   * and the inner exclusion zone fraction is hard-coded */
  rw  = 12;
  hemf = 1.2;

  /* The distances in physical space */
  rrad   = SPH1(x1, x2, x3);
  rr_cyl = CYL1(x1, x2, x3);
  rz     = CYL2(x1, x2, x3);

  /* The ellipse equation, lhs, with a smoothing factor */
  wrot = g_inputParam[PAR_WROT]*ini_code[PAR_WROT];
  wrad = g_inputParam[PAR_WRAD]*ini_code[PAR_WRAD];
  ellipse = (rr_cyl*rr_cyl*(1 - pow(wrot, 2)) + rz*rz)/
            (exp(-1)*pow(wrad,2));

  /* Inner hemisphere to keep free */
  urad = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  //uthk = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
  inner_circ = rrad/(hemf*(urad + 0.5*uthk));

  /* The inner ellipse equation, lhs, beyond which smoothing region begins */
  wsmf = g_inputParam[PAR_WSMF]*ini_code[PAR_WSMF];
  ellipse_i = (rr_cyl*rr_cyl*(1 - pow(wrot, 2)) + rz*rz)/
              ((1. - wsmf)*exp(-1)*pow(wrad,2));

  /* Exclude zone outside ellipsoid */
  if ( ellipse > 1. || inner_circ < 1.){
    fdratio = 1;
    is_cloud = 0;
  }

  /* The smoothing region */ 
  else if (!((ellipse_i < 1.) || (ellipse > 1.))){
    tanhfactor = tanh(tan(-CONST_PI*(sqrt(ellipse) - (1. - 0.5*wsmf))/wsmf));
    fdratio = 0.5*(fdratio + 1) + 0.5*tanhfactor*(fdratio - 1.);
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
  dens_lim = 0.0
  if (fd < dens_lim) {
    fdratio = 1.;
    is_cloud = 0;
  }
  fdratio = MAX(dens_lim, fdratio);
  is_cloud = 1;

#else
  is_cloud = 1;

  if (!once01){
    print1("> Cloud extraction: No special cloud extraction.\n\n");
    once01 = 1;
  }

#endif

  /* The ratio of cloud to halo density for a cell
   * shouldn't really ever be < 1. */
  fdratio = MAX(fdratio, 1.);

  /* Set to cloud density */
  cloud[RHO] = fdratio*halo[RHO];

  return is_cloud;
}

#endif


#if CLOUD_VELOCITY
/* ************************************************************** */
void CloudVelocity(double* cloud, double* halo, 
                   const double x1, const double x2, const double x3)
/*!
 * This function fills in the cloud velocity for the clouds primitives
 * array. If CLOUD_VELOCITY == NO, then there is no cube for velocities.
 * and initial velocity is, at first, 0 everywhere.
 *
 * On top of this, One of several modifications are applied. Currently
 * supported for CLOUD_VEL_DISTR are 
 *   - CV_UNIFORM: a constant, plane parallel velocity for each cell
 *   - CV_RADIAL: a constant, radial velocity for each cell
 *   - CV_KEPLERIAN_FRAC: a fraction CV_VALUE of the Keplerian velocity is added
 *   - CV_VIRIAL_FRAC: a fraction CV_VALUE of the virial velocity is added
 * 
 **************************************************************** */
{

  double vel, scrh;

#if CLOUD_VEL_DISTR == CV_UNIFORM 
  EXPAND(cloud[VX1] += CV_VALUE;,
         cloud[VX2] += CV_VALUE;,
         cloud[VX3] += CV_VALUE;);

#elif CLOUD_VEL_DISTR == CV_RADIAL 
  double v1, v2, v3;
  double vsph1, vsph2, vsph3;
  double xsph1, xsph2, xsph3;

  /* Convert coordinates and velocity vectors to spherical */

  EXPAND(xsph1 = SPH1(x1, x2, x3),
         xsph2 = SPH2(x1, x2, x3),
         xsph3 = SPH3(x1, x2, x3));

  v1 = cloud[VX1]; v2 = cloud[VX2]; v3 = cloud[VX3];

  EXPAND(vsph1 = VSPH1(x1, x2, x3, v1, v2, v3),
         vsph2 = VSPH2(x1, x2, x3, v1, v2, v3),
         vsph3 = VSPH3(x1, x2, x3, v1, v2, v3));
  
  /* Apply change to radial component */
  vsph1 += CV_VALUE;

  /* Convert velocity vectors back to the current coordinate system */
  EXPAND(v1 = VSPH_1(xsph1, xsph2, xsph3, vsph1, vsph2, vsph3),
         v2 = VSPH_2(xsph1, xsph2, xsph3, vsph1, vsph2, vsph3),
         v3 = VSPH_3(xsph1, xsph2, xsph3, vsph1, vsph2, vsph3));

  cloud[VX1] = v1; cloud[VX2] = v2; cloud[VX3] = v3;

#elif CLOUD_VEL_DISTR == CV_KEPLERIAN_FRAC 

  double v1, v2, v3;
  double vpol1, vpol2, vpol3;
  double xpol1, xpol2, xpol3;

  int il;
  double r_cyl, r1, r2, frac, y0, y1, y2, y3, phiderv_cyl;

  /* Convert coordinates and velocity vectors to polar */
  EXPAND(xpol1 = POL1(x1, x2, x3);,
         xpol2 = POL2(x1, x2, x3);,
         xpol3 = POL3(x1, x2, x3););

  v1 = cloud[VX1]; v2 = cloud[VX2]; v3 = cloud[VX3];

  EXPAND(vpol1 = VPOL1(x1, x2, x3, v1, v2, v3);,
         vpol2 = VPOL2(x1, x2, x3, v1, v2, v3);,
         vpol3 = VPOL3(x1, x2, x3, v1, v2, v3););
  
  /* Apply change to radial component */
  //vpol2 += CV_VALUE;

  /* Set mean vphi to keplerian: vphi = wrot * sqrt(rdphi / dr) */

  /* cycldrical radius */
  r_cyl = CYL1(x1, x2, x3);
 
  /* Find cell left index to interpolate at */
  il = hunter(gr_rad, gr_ndata, r_cyl);

  /* Find fraction of radius within cell */
  r1 = gr_rad[il];
  r2 = gr_rad[il+1];
  frac = (r_cyl - r1)/(r2 - r1);

  /* Interpolate */
  y0 = gr_phiderv[il-1]; 
  y1 = gr_phiderv[il];
  y2 = gr_phiderv[il+1];
  y3 = gr_phiderv[il+2];
  phiderv_cyl = CubicCatmullRomInterpolate(y0, y1, y2, y3, frac);

  /* The angular velocity */
  vpol2 += g_inputParam[PAR_WROT]*sqrt(r_cyl * phiderv_cyl);

  /* Convert velocity vectors back to the current coordinate system */
  EXPAND(v1 = VPOL_1(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3);,
         v2 = VPOL_2(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3);,
         v3 = VPOL_3(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3););

  EXPAND(cloud[VX1] = v1;, cloud[VX2] = v2;, cloud[VX3] = v3;);


#elif CLOUD_VEL_DISTR == CV_VIRIAL_FRAC
  /* Not yet programmed */
#endif

  /* Restrict cloud velocity to a fraction of halo sound speed */
  double cs_frac = 0.8;
  double cs = sqrt(SoundSpeed2IdealEOS(halo[RHO], halo[PRS]));
  double cs_lim = cs_frac*cs;
  EXPAND(cloud[VX1] = ABS_MIN(cs_lim, cloud[VX1]);,
         cloud[VX2] = ABS_MIN(cs_lim, cloud[VX2]);,
         cloud[VX3] = ABS_MIN(cs_lim, cloud[VX3]););

#if USE_FOUR_VELOCITY == YES
      vel = VMAG(x1, x2, x3, cloud[VX1], cloud[VX2], cloud[VX3]);
      scrh = Vel2Lorentz(vel);
      EXPAND(cloud[VX1] *= scrh;, cloud[VX2] *= scrh;, cloud[VX3] *= scrh;);
#endif
}
#endif


#if CLOUDS
/* ************************************************************** */
int CloudPrimitives(double* cloud, 
                    const double x1, const double x2, const double x3)
/*!
 * This function returns 1 if cell is cloud material, 0 if not.
 *
 * The array *cloud is filled in the process that involves three steps:
 * 0) Reading cloud data once. This fills cloud and gives us the fractal 
 *    factor "fd". Check that we're inside fractal cube domain.
 * 1) Multiply "fd" with a desired mean density profile. Currently
 *    supported are CLOUD_DISTR = 
 *        - CD_HERNQUIST
 *        - CD_TURB_ISOTH_HYDROSTATIC
 *        - CD_TURB_KEPLERIAN_DISC
 *        - CD_FLAT
 *    The default is also CD_FLAT. This is the apodization step. All 
 *    cells are still cloud material.
 *        We prefer to work with the density ratio, rho_cloud/rho_halo, mainly
 *    for the next extraction part. It gives us a good handle on how close
 *    we are to the halo density. 
 * 2) Use a form of cloud extraction by calling function CloudExtract.
 *    Some cells will not be cloud material anymore.
 * 3) Select, based on temperature criterion, whether cell is cloud material
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
  if (CloudCubePixel(cube_pixel, x1, x2, x3)){

    /* Get the fractal factor for this cell */
    GetFractalData(cloud, x1, x2, x3);

    /* Apodize cloud with a background mean density profile */
    CloudApodize(cloud, x1, x2, x3);

    /* Extract w.r.t. hot halo. Halo primitives are required 
     * for this step. */
    HotHaloPrimitives(halo, x1, x2, x3);
    if (CloudExtract(cloud, halo, cube_pixel, x1, x2, x3)){

      /* Calculate cloud primitives so that we can get temperature */

      /* Cloud velocity */
#if CLOUD_VELOCITY
      CloudVelocity(cloud, halo, x1, x2, x3);
#else
      EXPAND(cloud[VX1] = 0;,
             cloud[VX2] = 0;,
             cloud[VX3] = 0;);
#endif

      /* Cloud pressure
       * Underpressure the clouds slightly, so that they 
       * don't emit sound waves. 
       *   Note that the -xnjet option is available and that
       * the factor here should be compared the pressure gradient
       * calculated in SET_JET_DOMAIN jet_domain.c
       * */
      cloud[PRS] = halo[PRS]*CLOUD_UNDERPRESSURE;

      /* Tracers */
      cloud[TRC] = 0.0;
      cloud[TRC+1] = 1.0;

      /* Final test - is cloud pixel thermally stable? */
      is_cloud = WarmTcrit(cloud);

    }
  }

  /* Fill cloud array with halo primitves if not a cloud cell. This is not 
   * strictly necessary, since it is done outside CloudPrimitives, but we 
   * do it anyway for completeness. */
  if (is_cloud == 0) { for (nv = 0; nv < NVAR; ++nv) cloud[nv] = halo[nv]; }

  return is_cloud;
}
#endif


#if CLOUDS
/* ********************************************************************* */
int WarmTcrit(double * const warm)
/*! 
 * This routine returns 1 if cloud is still cloud, after thermal
 * stability criterion was set.
 *
 *********************************************************************** */
{
  double mu, wtemp;
  int nv;

  /* Warm phase temperature */
  mu = MeanMolecularWeight(warm);
  wtemp = TempIdealEOS(warm[RHO], warm[PRS], mu)*vn.temp_norm;


  /* Only a cloud pixel if wtemp is below critical temperature
   * of thermal instability */
#ifndef CLOUD_TCRIT
  fputs("Error: CLOUD_TCRIT not defined.\n", stderr); QUIT_PLUTO(1);
#endif

  if (wtemp > CLOUD_TCRIT){ 
    return 0;
  }
  else{
    return 1;
  }

}
#endif


/* ************************************************ */
int RotateGrid2Nozzle(
    double const cx1, double const cx2, double const cx3,
    double* cx1p, double* cx2p, double* cx3p){
/*
 * This function assumes 3D cartesian coor-
 * dinates. First rotates according to any precession 
 * through phi + omg*t. Then translates a point by dbh.
 * Then rotates grid around y-axis by angle dir so that
 * it is parallel to the nozzle cone.
 * Then instead of translating point back by dbh again,
 * a translation along the flow-axis by cbh is done to 
 * make the cap base flush with the computational domain.
 * 
 *
 * Thus, the transformation matrix is defined as
 * 
 *  / cx1p \              / cx1 \
 * |  cx2p  | = iT R T P |  cx2  |
 * |  cx3p  |            |  cx3  |
 *  \  1   /              \  1  /
 *
 *  where 
 *
 *  Inverse rotation about Z axis
 *
 *       / cos(pre)   sin(pre)  0  0 \
 *  P = |  -sin(pre)  cos(pre)  0  0  |
 *      |   0            0      1  0  |
 *       \  0            0      0  1 /
 *
 *       where pre = phi + omg*g_time
 *
 *  Inverse rotation about Y axis
 *
 *       /  cos(dir)   0  -sin(dir)  0 \
 *  R = |      0       1      0      0  |
 *      |   sin(dir)   0   cos(dir)  0  |
 *       \     0       0      0      1 /
 *
 *  Translation along z
 *
 *       / 1   0  0   0 \
 *  T = |  0   1  0   0  |
 *      |  0   0  1 -dbh |
 *       \ 0   0  0   1 /
 *
 *       / 1   0  0   0 \
 * iT = |  0   1  0   0  |
 *      |  0   0  1  dbh |
 *       \ 0   0  0   1 /
 *
 ************************************************** */

  /* Precession angle */
  double pre;
  pre = nz.phi + nz.omg*g_time;

  /* The transformations */
  // Sign for nz.dbh translations may be the wrong way around.
  
  SELECT(*cx1p = cx1;,
         *cx1p = cx1*cos(nz.dir) - sin(nz.dir)*(cx2 - nz.dbh);,
         *cx1p = (cx1*cos(pre) - cx2*sin(pre))*cos(nz.dir) - sin(nz.dir)*(cx3 - nz.dbh););

  SELECT(*cx2p = cx2;,
         *cx2p = cx1*sin(nz.dir) + nz.dbh + cos(nz.dir)*(cx2 - nz.dbh);,
         *cx2p = cx2*cos(pre) + cx1*sin(pre););
         
  SELECT(*cx3p = cx3;,
         *cx3p = cx3;,
         *cx3p = (cx1*cos(pre) - cx2*sin(pre))*sin(nz.dir) + nz.dbh + cos(nz.dir)*(cx3 - nz.dbh););

  return 0;
}

/* ************************************************ */
int RotateNozzle2Grid(
    double const cx1, double const cx2, double const cx3,
    double* cx1p, double* cx2p, double* cx3p){
/*
 * This Function does the inverse transformation of
 * RotateGrid2Nozzle
 *
 ************************************************** */

  /* Precession angle */
  double pre;
  pre = nz.phi + nz.omg*g_time;

  SELECT(*cx1p = cx1;,
         *cx1p = cx1*cos(nz.dir) + sin(nz.dir)*(cx2 + nz.dbh);,
         *cx1p = cx2*sin(pre) + cos(pre)*(cx1*cos(nz.dir) + sin(nz.dir)*(cx3 + nz.dbh)););
    
  SELECT(*cx2p = cx2;,
         *cx2p = -(cx1*sin(nz.dir)) - nz.dbh + cos(nz.dir)*(cx2 + nz.dbh);,
         *cx2p = cx2*cos(pre) - sin(pre)*(cx1*cos(nz.dir) + sin(nz.dir)*(cx3 + nz.dbh)););
    
  SELECT(*cx3p = cx3;,
         *cx3p = cx3;,
         *cx3p = -(cx1*sin(nz.dir)) - nz.dbh + cos(nz.dir)*(cx3 + nz.dbh););

  return 0;
}



/* ************************************************ */
int InNozzleRegion(double const x1, double const x2, double const x3){
/*
 * Returns 1 if r is in outflow region, 0 if not.
 * Nozzle is always a fan shaped region defined by
 * the intersection of a cone with the boundary volume.
 * A cap is included at the top of the cone, making it
 * an ice-cream.
 *   The outer rim of the maximum cone radius is at the 
 * surface of the computational domain. The cone apex is 
 * not necessarily at (0,0,0), neither is the tilt rotation
 * axis. The apex and tilt rotation axis are also not 
 * necessarily at the same point. 
 *    If the half opening angle
 * ang == 0, the ice-cream becomes a bullet.
 *
 ************************************************** */

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
    FLOWAXIS(cx1p, cx2p, cx3p) = fabs(FLOWAXIS(cx1p, cx2p, cx3p));
#endif
  FLOWAXIS(cx1p, cx2p, cx3p) -= nz.dbh + nz.cbh;

  /* Turn into cylindrical coords */
  double cr, cz;
  cr = CYL1(cx1p, cx2p, cx3p);
  cz = CYL2(cx1p, cx2p, cx3p);

  /* Spherical cap height */
  double cap_hgt;
#if INTERNAL_BOUNDARY == YES && HEMISPHERICAL_CAP == NO
    cap_hgt = sqrt(nz.sph*nz.sph - cr*cr) - nz.cbh;
#else
  cap_hgt = sqrt(nz.rad*nz.rad - cr*cr);
#endif

  int innozzle;

  /* Condition for Bullet */
  innozzle = cz < cap_hgt;

  /* Condition for ice-cream cone */
  if (nz.isfan) { 
    innozzle = innozzle && (cr < (cz + nz.cone_height)*tan(nz.ang));
  }
    

  return innozzle;

}


#if INTERNAL_BOUNDARY == YES
/* ************************************************ */
int InNozzleSphere(double const x1, double const x2, double const x3){
/*
 * Returns 1 if r is in sphere containinng nozzle.
 * Only for INTERNAL_BOUNDARY YES
 * Sphere is assumed to be at (0,0,0) and has a radius
 * of ODBH. The base of the nozzle has a radius of ORAD.
 ************************************************** */
    
    /* Grid point in cartesian coordinates */
    double cx1, cx2, cx3;
    cx1 = CART1(x1, x2, x3);
    cx2 = CART2(x1, x2, x3);
    cx3 = CART3(x1, x2, x3);
    
    /* Rotate cartesian coords so that base of nozzle cone is at (0,0,0) */
    double cx1p, cx2p, cx3p;
    RotateGrid2Nozzle(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);
    
    /* Turn into spherical coords */
    double sr;
    sr = SPH1(cx1p, cx2p, cx3p);
    
    int insphere;
    insphere = (sr < nz.sph);
    
    return insphere;
    
}
#endif


#if INTERNAL_BOUNDARY == YES
/* ************************************************ */
int SphereIntersectsDomain(Grid *grid){
/*
 * Returns 1 if sphere intersects local domain.
 * Only for INTERNAL_BOUNDARY YES
 *
 * grid is the pointer to the grid struct of the local domain.
 *
 * Sphere is assumed to be at (0,0,0) and has a radius
 * of ODBH. The base of the nozzle has a radius of ORAD.
 ************************************************** */
    

    /* The local domain limits */
    double x1i, x1f, x2i, x2f, x3i, x3f, r;
    
    /* Location of center of sphere - currently this is hardcoded here */
    double s1, s2, s3;
    s1 = s2 = s3 = 0;
    
    x1i = grid[IDIR].xi;
    x1f = grid[IDIR].xf;
    x2i = grid[JDIR].xi;
    x2f = grid[JDIR].xf;
    x3i = grid[KDIR].xi;
    x3f = grid[KDIR].xf;
    r = nz.sph;
    
    /* Check first whether the center of the sphere
     is in the extended (by r) domain box or
     the center of the sphere lies in any
     of the eight rounded corners */
    
    if ((x1i - r < s1 && s1 < x1f + r &&
         x2i - r < s2 && s2 < x2f + r &&
         x3i - r < s3 && s3 < x3f + r) ||
        ((x1i*x1i + x2i*x2i + x3i*x3i < r*r && x1i > s1 && x2i > s2 && x3i > s3) ||
         (x1i*x1i + x2i*x2i + x3f*x3f < r*r && x1i > s1 && x2i > s2 && x3f < s3) ||
         (x1i*x1i + x2f*x2f + x3f*x3f < r*r && x1i > s1 && x2f < s2 && x3f < s3) ||
         (x1i*x1i + x2f*x2f + x3i*x3i < r*r && x1i > s1 && x2f < s2 && x3i > s3) ||
         (x1f*x1f + x2i*x2i + x3i*x3i < r*r && x1f < s1 && x2i > s2 && x3i > s3) ||
         (x1f*x1f + x2i*x2i + x3f*x3f < r*r && x1f < s1 && x2i > s2 && x3f < s3) ||
         (x1f*x1f + x2f*x2f + x3i*x3i < r*r && x1f < s1 && x2f < s2 && x3i > s3) ||
         (x1f*x1f + x2f*x2f + x3f*x3f < r*r && x1f < s1 && x2f < s2 && x3f < s3))) {
            return 1;
    }
    else return 0;
}
#endif



/* ************************************************ */
double Vel2Lorentz(const double vel)
/*!
 * Return Lorentz factor from Lorentz factor
 *
 ************************************************** */
{
  double vsqr = 0, clight2;

  clight2 = pow(CONST_c/UNIT_VELOCITY, 2);

  vsqr = vel*vel;
  return 1.0/sqrt(1.0 - vsqr/clight2);
}


/* ************************************************ */
double Lorentz2Vel(const double lorentz) 
/*!
 * Return velocity from Lorentz factor
 *
 ************************************************** */
{
  double clight;
  clight = CONST_c/UNIT_VELOCITY;
  return sqrt(1. - 1./(lorentz*lorentz))*clight;
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
  cr = CYL1(cx1p, cx2p, cx3p);
  cz = CYL2(cx1p, cx2p, cx3p);

  /* Return smoothing factor */
  return 1.0 / cosh(pow(cr / nz.rad, n));
  }
}


/* ************************************************ */
double Profile_cap(const double x1, const double x2, const double x3)
/*! 
  * 1/cosh smoothing function
  *
 ************************************************** */
{
  /* Steepness of cosh profile */
  int n = 5.;
  int prof;

  /* Grid point in cartesian coordinates */
  double cx1, cx2, cx3, cx1p, cx2p, cx3p;
  cx1 = CART1(x1, x2, x3);
  cx2 = CART2(x1, x2, x3);
  cx3 = CART3(x1, x2, x3);

  /* Rotate cartesian coords and turn into cylindrical coords */
  double cr, cz;
  RotateGrid2Nozzle(cx1, cx2, cx3, &cx1p, &cx2p, &cx3p);

  /* TODO: Assumes origin at 0,0,0. Need to fix */
  cr = SPH1(cx1p, cx2p, cx3p);

  /* Return smoothing factor */
  return 1.0 / cosh(pow(cr / nz.rad, n));
  return prof;
}

