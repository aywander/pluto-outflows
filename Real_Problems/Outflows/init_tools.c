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
#include "nrEOS.h"
#include "abundances.h"
#include "interpolation.h"

/* Base normalization struct */
VarNorm vn;
double ini_cgs[32];
double ini_code[32];

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
 * Sets g_unit*  and initializes VarNorm struct with derived normalizations
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
  vn.pres_norm    =   vn.dens_norm* vn.v_norm*vn.v_norm;
  vn.power_norm   =   vn.pres_norm*vn.v_norm*vn.area_norm;
  vn.eflux_norm   =   vn.pres_norm*vn.v_norm;
  vn.eint_norm    =   vn.pres_norm/vn.dens_norm;
  vn.mdot_norm    =   vn.dens_norm*pow(vn.l_norm,3)/vn.t_norm;

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

  double year;

  year = CONST_ly/CONST_c;

  ini_cgs[PAR_OPOW] = 1.;                                  ini_code[PAR_OPOW] = ini_cgs[PAR_OPOW]/vn.power_norm;
  ini_cgs[PAR_OANG] = CONST_PI/180.;                       ini_code[PAR_OANG] = ini_cgs[PAR_OANG];
  ini_cgs[PAR_OSPD] = NOZZLE_SELECT(1.,CONST_c);           ini_code[PAR_OSPD] = ini_cgs[PAR_OSPD]/vn.v_norm;
  ini_cgs[PAR_OMDT] = NOZZLE_SELECT(1.,CONST_Msun/year);   ini_code[PAR_OMDT] = ini_cgs[PAR_OMDT]/vn.mdot_norm;
  ini_cgs[PAR_ORAD] = vn.l_norm;                           ini_code[PAR_ORAD] = ini_cgs[PAR_ORAD]/vn.l_norm;
  ini_cgs[PAR_OTHK] = vn.l_norm;                           ini_code[PAR_OTHK] = ini_cgs[PAR_OTHK]/vn.l_norm;
  ini_cgs[PAR_ODIR] = 1;                                   ini_code[PAR_ODIR] = ini_cgs[PAR_ODIR];
  ini_cgs[PAR_HRHO] = vn.dens_norm;                        ini_code[PAR_HRHO] = ini_cgs[PAR_HRHO]/vn.dens_norm;
  ini_cgs[PAR_HTE ] = 1;                                   ini_code[PAR_HTE ] = ini_cgs[PAR_HTE ]/vn.temp_norm;
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

  NOZZLE_SELECT(JetPrimitives(out_primitives, x1, x2, x3),
                UfoPrimitives(out_primitives, x1, x2, x3));

}


/* ************************************************ */
void OutflowVelocity(double * out_primitives, double speed,
    const double x1, const double x2, const double x3, const double th){
/*
 * Calculate outflow velocity vector inside out_primitives given 
 *
 * - a cell location x1, x2, x3
 * - angle of cell velocity vector theta (not necessarily = OANG parameter)
 *
 * This function is used by UfoPrimitives and JetPrimitives
 * and should be generic for any outflow nozzle.
 *
 ************************************************** */

  double xs1, xs2, xs3;
  double vs1, vs2, vs3;
#if USE_FOUR_VELOCITY == YES
  double lorentz, vel;
  int iv;
#endif


#if USE_FOUR_VELOCITY == YES
  lorentz = Vel2Lorentz(speed);
#endif


  /* This is for both NOZZLE_SHAPE = ANN and FAN 
   * Work in spherical. Assume the coordinate system describing the
   * velocity is shifted downward by w/tan(th). The azimuthal angle 
   * doesn't change. Special treatment is necessary for geometries
   * that don't have vector components that are invariatne under 
   * translation along z axis, e.g. spherical. */
  xs3 = SPH3(x1, x2, x3);

  /* Convert to current geomtetry */
#if GEOMETRY == SPHERICAL

 

  double vca1, vca2, vca3;
  double xca1, xca2, xca3;

  /*
  EXPAND(vca1 = VSPH2CART1(0, th, xs3, speed, 0, 0);,
         vca2 = VSPH2CART2(0, th, xs3, speed, 0, 0);,
         vca3 = VSPH2CART3(0, th, xs3, speed, 0, 0););

  xca1 = CART1(x1, x2, x3);
  xca2 = CART2(x1, x2, x3);
  xca3 = CART3(x1, x2, x3);

  EXPAND(out_primitives[VX1] = VCART2SPH1(xca1, xca2, xca3, vca1, vca2, vca3);,
         out_primitives[VX2] = VCART2SPH2(xca1, xca2, xca3, vca1, vca2, vca3);,
         out_primitives[VX3] = VCART2SPH3(xca1, xca2, xca3, vca1, vca2, vca3);)
 
   */
    

  EXPAND(out_primitives[VX1] =  speed*cos(x2 - th);,
         out_primitives[VX2] =  speed*sin(x2 - th);,
         out_primitives[VX3] = 0;)


#else
    /* OK in all geometries except spherical and 2D polar */
  EXPAND(out_primitives[VX1] = VSPH_1(1, th, xs3, speed, 0, 0);,
         out_primitives[VX2] = VSPH_2(1, th, xs3, speed, 0, 0);,
         out_primitives[VX3] = VSPH_3(1, th, xs3, speed, 0, 0););
#endif


#if USE_FOUR_VELOCITY == YES
  /* This is the same for all geometries */
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

  double power, chi, lorentz, radius; 
  double dens, pres, vel;
  double area, gmm1;

  /* The parameters from ini, and normalize them */
  power   = g_inputParam[PAR_OPOW]*ini_code[PAR_OPOW];
  chi     = g_inputParam[PAR_OMDT]*ini_code[PAR_OMDT];
  lorentz = g_inputParam[PAR_OSPD]*ini_code[PAR_OSPD];
  radius  = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];

  /* Some needed quantities */
  vel = Lorentz2Vel(lorentz);
  area = CONST_PI*radius*radius;

  /* m_gamma is global in pluto.h */
  gmm1 = g_gamma/(g_gamma - 1.);
  pres = power/(gmm1*lorentz*lorentz*vel*area*
      (1. + (lorentz - 1.)/lorentz*chi));
  dens = chi*gmm1*pres;


  /* Set primitives array */
  jet_primitives[RHO] = dens;
  jet_primitives[PRS] = pres;
  jet_primitives[TRC] = 1.0;
#if CLOUDS
  jet_primitives[TRC+1] = 0.0;
#endif
  OutflowVelocity(jet_primitives, vel, x1, x2, x3, 0.0);

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

  double power, speed, mdot, r0, w, th; 
  double dens, pres;
  double area;

  /* Input parameters, normalized */
  r0    = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  w     = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
  th    = g_inputParam[PAR_OANG]*ini_code[PAR_OANG];
  power = g_inputParam[PAR_OPOW]*ini_code[PAR_OPOW];
  speed = g_inputParam[PAR_OSPD]*ini_code[PAR_OSPD];
  mdot  = g_inputParam[PAR_OMDT]*ini_code[PAR_OMDT];

  /* Derived quantities */
  area = 2*CONST_PI*r0*w;
  pres = (power - 0.5*mdot*speed*speed)*(g_gamma-1)/(g_gamma*area*speed);
  dens = mdot/(area*speed);


  /* Set primitives array */
  ufo_primitives[RHO] = dens;
  ufo_primitives[PRS] = pres;
  ufo_primitives[TRC] = 1.0;
#if CLOUDS
  ufo_primitives[TRC+1] = 0.0;
#endif

#if NOZZLE_SHAPE == FAN
  /* w is assumed not to be angled */
  double theta, c1, c2, k;
  c1 = CYL1(x1, x2, x3);
  c2 = CYL1(x1, x2, x3);
  k = w/(2*tan(th));
  theta = atan(c1/(c2 + k));
  OutflowVelocity(ufo_primitives, speed, x1, x2, x3, theta);
#elif NOZZLE_SHAPE == ANN
  OutflowVelocity(ufo_primitives, speed, x1, x2, x3, th);
#endif

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
  double unit_time, inv_unit_G;
  double r, a, rs, rho0;
  double frac;
  double y0, y1, y2, y3, r1, r2; 
  int iv, il, ic;
  static int once01 = 0;
  //static int once02 = 0;
  //static int once03 = 0;


  /* Initialize gravity arrays - as good as any other place to do it */
#ifdef GRAV_TABLE
  //if (!once03){
  if (gr_rad == NULL){
    readGravTable();
    //once03 = 1;
  }
#endif

  /* Consider different distributions */

  /* Hernquist potential (hydrostatic)*/
#if HOT_DISTR == HOT_HERNQUIST
  rho0 = g_inputParam[PAR_HRHO]*ini_code[PAR_HRHO];
  r = SPH1(x1, x2, x3);
  a = g_inputParam[PAR_HRAD]*ini_code[PAR_HRAD];
  rs = r/a;
  unit_time = UNIT_LENGTH/UNIT_VELOCITY;
  inv_unit_G = unit_time*unit_time*UNIT_DENSITY;
  halo[RHO] = rho0/(rs*pow((1 + rs), 3));
  halo[PRS] = -2*CONST_PI*CONST_G*inv_unit_G*rho0*rho0*a*a*
    (pow(1.+rs, -4)/12.*(25. + 52.*rs + 42.*rs*rs + 12.*rs*rs*rs) -
    log(1. + rs) + log(rs));


  /* Any spherical potential read from a table */
#elif HOT_DISTR == HOT_EXT_DATA
  
  if (hot_rad == NULL){
  //if (!once02){
    readHotTable();
    //once02 = 1;
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
  halo[PRS] = PresNrEOS(halo[RHO], g_inputParam[PAR_HTE]*ini_code[PAR_HTE], MU_NORM);


  /* Default is flat */
#else

  halo[RHO] = g_inputParam[PAR_HRHO]*ini_code[PAR_HRHO];
  halo[PRS] = PresNrEOS(halo[RHO], g_inputParam[PAR_HTE]*ini_code[PAR_HTE], MU_NORM);

#endif

   /* A message for having initialized halo with potential*/
  if (!once01){
    print1("> Initializing hot halo distribution of type: %d\n\n", HOT_DISTR);
    once01 = 1;
  }

  /* Velocities. 
   * NOTE: At some point we may want to generalize the 
   * background velocity assignment. This is mainly motivated by 
   * single cloud acceleration setups, but could be used in 
   * simulations of cluster radio galaxies. */
  EXPAND(halo[VX1] = 0;,halo[VX2] = 0;,
         halo[VX3] = g_inputParam[PAR_HVBG]*ini_code[PAR_HVBG];);
#if USE_FOUR_VELOCITY == YES
  vel = VMAG(x1, x2, x3, halo[VX1], halo[VX2], halo[VX3]);
  scrh = Vel2Lorentz(vel);
  EXPAND(halo[VX1] *= scrh;, halo[VX2] *= scrh;, halo[VX3] *= scrh;);
#endif

  /* Tracer */
  halo[TRC] = 0.0;
#if CLOUDS
  halo[TRC+1] = 0.0;
#endif

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
 * Calculate pixel coordinates in fractal cube for given zone.
 * The origin of the pixel coordinates (here el[]) is at the 
 * corner with the smallest values of x1, x2, and x3 
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

  int get_var[] = {RHO, -1};
#if CLOUD_VELOCITY
#define ID_MAX_NVAR (1 + COMPONENTS + 1)
  int get_var[ID_MAX_NVAR]
  get_var[0] = RHO;
  EXPAND(get_var[1] = VX1;
         get_var[2] = VX2;
         get_var[3] = VX3;);
  get_var[ID_MAX_NVAR-1] = -1;
#endif

  static int once01 = 0;

  /* Read cloud data from external file */
  if (!once01) {
    InputDataSet("./grid_in.out", get_var);
    InputDataRead("./clouds.dbl", CUBE_ENDIANNESS);
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
  double sigma_g2, sigma_turb2;
  double mu, dens;
  double wte, wrad, wrho, wtrb, wrot;

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

  /* The total dense gas velocity dispersion */
  wrho = g_inputParam[PAR_WRHO]*ini_code[PAR_WRHO];
  wrad = g_inputParam[PAR_WRAD]*ini_cgs[PAR_WRAD];
  sigma_g2 = 4*CONST_PI*CONST_G*wrho*pow(wrad, 2)/9.;

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

  /* Turbulent velocity dispersion - not explicitly used */
  //wtrb = g_inputParam[PAR_WTRB]*ini_cgs[PAR_WTRB];
  //sigma_turb2 = wtrb*wtrb;

  /* Mean warm phase temperature. Not explicitly used, 
   * just for reference. Note, we cannot calculate a variable
   * mu because cloud is imcomplete. */
  //mu = MeanMolecularWeight(cloud);
  //wte = sqrt(sigma_g2 - sigma_turb2)*mu*CONST_amu/CONST_kB;


#elif CLOUD_DISTR == CD_TURB_KEPLERIAN_DISC

  /* Must use gravitational potential */
#if BODY_FORCE != POTENTIAL
  fputs("Error: Must use BODY_FORCE = POTENTIAL\n", stderr); 
  fputs("with CLOUD_DISTR = CD_TURB_KEPLERIAN_DISC.\n", stderr); 
  QUIT_PLUTO(1);
#endif


  /* The total dense gas velocity dispersion */
  wrho = g_inputParam[PAR_WRHO]*ini_code[PAR_WRHO];
  wrad = g_inputParam[PAR_WRAD]*ini_cgs[PAR_WRAD];
  sigma_g2 = 4*CONST_PI*CONST_G*wrho*pow(wrad, 2)/9.;

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
  r1 = gr_rad[il]; r2 = gr_rad[il+1];
  frac = (r_cyl - r1)/(r2 - r1);
  y0 = gr_phi[il-1]; y1 = gr_phi[il]; y2 = gr_phi[il+1]; y3 = gr_phi[il+2];
  phi_cyl = CubicCatmullRomInterpolate(y0, y1, y2, y3, frac);

  /* The profile */
  wrot = g_inputParam[PAR_WROT]*g_inputParam[PAR_WROT];
  dens = wrho*exp((-phi + phi_cyl*wrot)/sigma_g2);

  /* Turbulent velocity dispersion - not explicitly used*/
  wtrb = g_inputParam[PAR_WTRB]*ini_cgs[PAR_WTRB];
  sigma_turb2 = wtrb*wtrb;

  /* Mean warm phase temperature. Not explicitly used, 
   * just for reference. Note, we cannot calculate a variable
   * mu because cloud is imcomplete. */
  //mu = MeanMolecularWeight(cloud);
  //wte = sqrt(sigma_g2 - sigma_turb2)*mu*CONST_amu/CONST_kB;


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

  /* We are working in cartesian. Clouds cube is always
   * assumed to be a data file in cartesian */
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
  uthk = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
  inner_circ = rrad/(hemf*(urad + uthk));

  /* Exclude zone outside interior sphere */
  wrad = g_inputParam[PAR_WRAD]*ini_code[PAR_WRAD]
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
  uthk = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
#if NOZZLE_SHAPE == FAN
  inner_circ = rrad/(hemf*(urad + 0.5*uthk));
#else
  inner_circ = rrad/(hemf*(urad + uthk));
#endif

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
  uthk = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
#if NOZZLE_SHAPE == FAN
  inner_circ = rrad/(hemf*(urad + 0.5*uthk));
#else
  inner_circ = rrad/(hemf*(urad + uthk));
#endif

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


#if CLOUDS
/* ************************************************************** */
void CloudVelocity(double* cloud,
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

#if !CLOUD_VELOCITY
  EXPAND(cloud[VX1] = 0;,
         cloud[VX2] = 0;,
         cloud[VX3] = 0;);
#endif

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

  /* Convert coordinates and velocity vectors to polerical */

  EXPAND(xpol1 = POL1(x1, x2, x3),
         xpol2 = POL2(x1, x2, x3),
         xpol3 = POL3(x1, x2, x3));

  v1 = cloud[VX1]; v2 = cloud[VX2]; v3 = cloud[VX3];

  EXPAND(vpol1 = VPOL1(x1, x2, x3, v1, v2, v3),
         vpol2 = VPOL2(x1, x2, x3, v1, v2, v3),
         vpol3 = VPOL3(x1, x2, x3, v1, v2, v3));
  
  /* Apply change to radial component */
  vpol2 += CV_VALUE;

  /* Convert velocity vectors back to the current coordinate system */
  EXPAND(v1 = VPOL_1(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3),
         v2 = VPOL_2(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3),
         v3 = VPOL_3(xpol1, xpol2, xpol3, vpol1, vpol2, vpol3));

  cloud[VX1] = v1; cloud[VX2] = v2; cloud[VX3] = v3;


#elif CLOUD_VEL_DISTR == CV_VIRIAL_FRAC
  /* Not yet programmed */
#endif


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
      CloudVelocity(cloud, x1, x2, x3);

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
  wtemp = TempNrEOS(warm[RHO], warm[PRS], mu)*vn.temp_norm;


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
int RotateGrid2Nozzle(double const r , double const x,
                             double*      rp, double*      xp){
/*
 * This function assumes 2D (cartesian) coor-
 * dinates for a cylindrically symmetric nozzle.
 * First translates a point by -k = r/sin(theta). 
 * Rotates a point from coordinates defined by grid to 
 * coordinates parallel and perpendicular to the theta
 * (OANG parameter) direction. 
 * Then translates point back by k again.
 * The rotation is a 2D rotation by an angle
 * -theta about . Thus, the rotation matrix is defined as
 * 
 *  / rp \       / r \
 * |      | = R |     |
 *  \ xp /       \ y /
 *
 *  where 
 *
 *       / cos(th)   sin(th) \
 *  R = |                     |
 *       \ -sin(th)  cos(th) /
 *
 ************************************************** */

  double th;

  th = g_inputParam[PAR_OANG]*ini_code[PAR_OANG];

  *rp =  cos(th)*r - sin(th)*x;
  *xp =  sin(th)*r + cos(th)*x;

  return 0;
}


/* ************************************************ */
int RotateNozzle2Grid(double const r , double const x,
                             double* rp, double*  xp){
/*
 * This function assumes 2D (cartesian) coor-
 * dinates for a cylindrically symmetric nozzle.
 * Rotates a point from coordinates defined by
 * coordinates parallel and perpendicular to the theta
 * (OANG parameter) direction, to coordinates aligned with grid.
 * The rotation is a 2D rotation by an angle
 * theta around (0,0,0).
 * Thus, the rotation matrix is defined as
 * 
 *  / r \       / rp \
 * |     | = R |      |
 *  \ x /       \ yp /
 *
 *  where 
 *       / cos(th)   sin(th) \
 *  R = |                     |
 *       \ -sin(th)  cos(th) /
 *
 ************************************************** */

  double th;

  th = g_inputParam[PAR_OANG]*ini_code[PAR_OANG];

  *rp =  cos(th)*r + sin(th)*x;
  *xp = -sin(th)*r + cos(th)*x;

  return 0;
}




/* ************************************************ */
int TranslateRotateGrid2Nozzle(double const r , double const x,
                               double*      rp, double*      xp){
/*
 * This function assumes 2D (cartesian) coor-
 * dinates for a cylindrically symmetric nozzle.
 * First translates a point by -k = r/sin(theta). 
 * Rotates a point from coordinates defined by grid to 
 * coordinates parallel and perpendicular to the theta
 * (OANG parameter) direction. 
 * Then translates point back by k again.
 * The rotation is a 2D rotation by an angle
 * -theta about . Thus, the transformation is 
 * 
 *  / rp \       / r     \     / 0 \
 * |      | = R |         | - |    |
 *  \ xp /       \ x + k /     \ k /
 *
 *  where 
 *       / cos(th)   sin(th) \
 *  R = |                     |
 *       \ -sin(th)  cos(th) /
 *
 ************************************************** */

  double th, w, k;

  th = g_inputParam[PAR_OANG]*ini_code[PAR_OANG];
#if NOZZLE_SHAPE == ANN
  w  = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
#elif NOZZLE_SHAPE == FAN
  w  = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
#endif
  k = w*sin(th);

  /* Offset: from (0, 0) to center of rotation*/

  *rp =  cos(th)*r - sin(th)*(x + k);
  *xp =  sin(th)*r + cos(th)*(x + k) - k;

  return 0;
}


/* ************************************************ */
int TranslateRotateNozzle2Grid(double const r , double const x,
                             double*      rp, double*      xp){
/*
 * This function assumes 2D (cartesian) coor-
 * dinates for a cylindrically symmetric nozzle.
 * First translates a point by -k = r/sin(theta).
 * Then rotates a point from coordinates defined by
 * coordinates parallel and perpendicular to the theta
 * (OANG parameter) direction, to coordinates aligned with grid.
 * Then translates point back by k again.
 * The rotation is a 2D rotation by an angle
 * theta around a point -k = r/sin(theta).
 * Thus, the rotation matrix is defined as
 * 
 *  / r \       / rp     \     / 0 \
 * |     | = R |          | - |    |
 *  \ x /       \ xp + k /     \ k /
 *
 *  where 
 *
 *       / cos(th)   sin(th) \
 *  R = |                     |
 *       \ -sin(th)  cos(th) /
 *
 ************************************************** */

  double th, w, k;
  th = g_inputParam[PAR_OANG]*ini_code[PAR_OANG];
#if NOZZLE_SHAPE == ANN
  w  = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
#elif NOZZLE_SHAPE == FAN
  w  = 0.5*g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
#endif
  k = w/tan(th);


  *rp =  cos(th)*r + sin(th)*(x + k);
  *xp = -sin(th)*r + cos(th)*(x + k) - k;

  return 0;
}






/* ************************************************ */
int InNozzleBase(double const x1, double const x2, double const x3){
/*
 * Returns 1 if r is in wind base on x=0 plane, 0 if not. The wind
 * region is the hemisphere defined by r2.
 *
 * A FAIRE:
 * Add special condition for when th < small_value.
 *
 ************************************************** */

  double w, w_c, r0, th, r2, r1, r, z;
  const double small = 1.e-10;

  th = g_inputParam[PAR_OANG]*ini_code[PAR_OANG];
  w = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];

  r = CYL1(x1, x2, x3);
  z = CYL2(x1, x2, x3);

#if NOZZLE_SHAPE == FAN
   
#if GEOMETRY == SPHERICAL
  /* Assumes g_domBeg[IDIR] > w/2, and that w = w_c */
  double k, theta, rsph;
  k = w/(2*tan(th));
  theta = atan(r/(z + k));
  rsph = SPH1(x1, x2, x3);
  return ((fabs(rsph - g_domBeg[IDIR]) < small) && (theta <= th));
#else
  return ((fabs(z) < small) && (r < 0.5*w));
#endif

#elif NOZZLE_SHAPE == ANN
  w_c = w/cos(th);
  r0 = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  r2 = r0 + 0.5*w_c;
  r1 = r0 - 0.5*w_c;

#if GEOMETRY == SPHERICAL
  double k, theta, rsph, thsph, r1p, z1p;
  k = r0/(2*tan(th));
  theta = atan(r/(z + k));

  rsph = SPH1(x1, x2, x3);
  thsph = SPH2(x1, x2, x3);

  RotateGrid2Nozzle(r1, 0, &r1p, &z1p);

  return ((fabs(rsph - g_domBeg[IDIR]) < small) && 
      (thsph >= th + asin(r1p/g_domBeg[IDIR])) && 
      (theta <= th));
#else
  return ((fabs(z) < small) && (r < r2 ) && (r > r1));
#endif

#endif
}



/* ************************************************ */
int InNozzleRegion(double const x1, double const x2, double const x3){
/*
 * Returns 1 if r is in outflow region, 0 if not. The outflow
 * region is the wedge and cap of radius w/2 in the 
 * direction of the outflow. k is the distance from (0,0,0)
 * to the point where the hypothenuse of teh wedge meets the 
 * outflow axis.
 *
 * A FAIRE:
 * Add special condition for when th < small_value.
 *
 ************************************************** */

  double r, r0, r1, r2, th, w, bullet, z;
  double rp, r0p, r1p, r2p, zp, z0p, z1p, z2p;

  th = g_inputParam[PAR_OANG]*ini_code[PAR_OANG];
  w  = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
  
  r = CYL1(x1, x2, x3);
  z = CYL2(x1, x2, x3);


#if NOZZLE_SHAPE == FAN
  double k, rr, d;
  k = w/(2*tan(th));
  rr = sqrt((z + k)*(z + k) + r*r);
  d = 0.5*w/sin(th);
   
#if GEOMETRY == SPHERICAL
  /* Assumes g_domBeg[IDIR] > w/2 */
    return (rr <= g_domBeg[IDIR] + d);
#else
    return (rr <= d);
#endif


#elif NOZZLE_SHAPE == ANN

  r0 = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  r1 = r0 - 0.5*w/cos(th);
  r2 = r0 + 0.5*w/cos(th);

  RotateGrid2Nozzle(r , z, &rp , &zp);
  RotateGrid2Nozzle(r0, 0, &r0p, &z0p);
  RotateGrid2Nozzle(r1, 0, &r1p, &z1p);
  RotateGrid2Nozzle(r2, 0, &r2p, &z2p);

#if GEOMETRY == SPHERICAL
  bullet = sqrt(r1p*r1p + g_domBeg[IDIR]*g_domBeg[IDIR]) + 
    sqrt(w*w/4. - pow((rp - r0p), 2));
#else
  bullet = z2p + sqrt(w*w/4. - pow((rp - r0p), 2));
#endif

  return (zp <= bullet);
#endif
}


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
 * Cylindrically symmetric wind
 * As a function of three coordinates.
 * 
 ************************************************** */
{
  int n = 14;
  double r, z, rp, zp, r0, r0p, z0p, w;

  /* Assumes nozzle axis is along first directional component */
  r = CYL1(x1, x2, x3);
  z = CYL2(x1, x2, x3);

  w  = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
#if NOZZLE_SHAPE == ANN
  r0 = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  RotateGrid2Nozzle(r , z, &rp , &zp);
  RotateGrid2Nozzle(r0, 0, &r0p, &z0p);
#elif NOZZLE_SHAPE == FAN
  rp = r;
  r0p = 0;
#endif


  return 1.0/cosh(pow((rp - r0p)/(w/2.), n));
}



/* ************************************************ */
double Profile_sharp(const double x1, const double x2, const double x3)
/* 
 * Cylindrically symmetric wind
 * r      radial coordinate 
 * x      coordinate perpendicular to disc
 *
 ************************************************** */
{
  double r, z, rp, zp, r0, r0p, z0p, w;

  /* Assumes nozzle axis is along first directional component */
  r = CYL1(x1, x2, x3);
  z = CYL2(x1, x2, x3);

  w  = g_inputParam[PAR_OTHK]*ini_code[PAR_OTHK];
#if NOZZLE_SHAPE == ANN
  r0 = g_inputParam[PAR_ORAD]*ini_code[PAR_ORAD];
  RotateGrid2Nozzle(r , z, &rp , &zp);
  RotateGrid2Nozzle(r0, 0, &r0p, &z0p);
#elif NOZZLE_SHAPE == FAN
  rp = r;
  r0p = 0;
#endif


  return (fabs(rp - r0p) <= w/2.);

}




