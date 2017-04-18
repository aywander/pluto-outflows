#ifndef PLUTO_USR_H
#define PLUTO_USR_H

/* For completeness, parameters which take arbitrary constants, 
 * or just YES or NO are listed here for reference in comments. */

/* Cooling */

/* MU_CALC values. Methods of calculating the mean molecular mass. */
#define MU_CONST        0
#define MU_TABLE        1
#define MU_ANALYTIC     2
#define MU_FRACTIONS    3
#define MU_FNAME "mutable.dat"


/* Inlet */

/* NOZZLE values. 
 * Type of inlet nozzle 
 * This determines e.g. whether UfoPrimitives, 
 * JetPrimitives, etc are called.
 * */
#define NOZZLE_JET      1
#define NOZZLE_UFO      2

/* NOZZLE_FILL values.
 * How nozzle region is set.
 * Either set primitives by overwriting them
 * or dump energy and mass as conservative quantities
 * */
#define NF_PRIMITIVE    1
#define NF_CONSERVATIVE 2


/* Accretion */
/* SIC_METHOD values.
 * See SphereSurfaceIntersectsCell function. */
#define SIC_RADIUS  1
#define SIC_CORNERS 2
#define SIC_HYBRID  3

/* SID_METHOD values.
 * See SphereIntersectsDomain function. */
#define SID_POINTS  1
#define SID_REGIONS 2

/* Sink */
/* SINK_METHOD values */
#define SINK_VACUUM 1
#define SINK_FREEFLOW 2
#define SINK_BONDI 3
#define SINK_FEDERRATH 4

/* Feedback Cycle */
/* Feedback cycle modes */
/* FBC_DEBOOST_MODE values. */
#define FBC_DEBOOST_MODE_0  0
#define FBC_DEBOOST_MODE_1  1
#define FBC_DEBOOST_MODE_2  2
#define FBC_DEBOOST_MODE_3  3

/* Supernovae */
/* SUPERNOVAE implementations */
// TODO: Implement different SN methods

/* Clouds */
/* Use a grid_in.out file to specify dimensions of cube 
 * Input files are rho.dbl, vx1.dbl, vx2.dbl, vx3.dbl, etc*/

/* CLOUD_EXTRACT values. 
 * Method of cloud extraction 
 * DEFAULT is NONE
 * Shape of ELLIPSOID is determined by parameter WROT and WRAD
 * */
// TODO: Test and debug CE_ELLIPSOID and CE_DENS
#define CE_DENS              1
#define CE_ELLIPSOID         2

/* CLOUD_DENSITY values.
 * Cloud distribution 
 * The density profile for the warm phase */
#define CD_HOMOGENEOUS          0
#define CD_KEPLERIAN            1
#define CD_HERNQUIST            2

/* CLOUD_VELOCITY values.
 * Cloud velocity distribution.
 * The velocity cube can be modified in different ways.
 * The units of velocity should be km/s.
 * Use ini_code[WTRB] to normalize */
#define CV_ZERO          0
#define CV_KEPLERIAN     1

/* CLOUD_SCALE values.
 * Decide whether the cloud scale height is given by the velocity
 * dispersion or the radius */
#define CS_SCALE_HEIGHT          0
#define CS_VELOCITY_DISPERSION   1


/* Static Gravity */

/* GRAV_POTENTIAL values. 
 * Gravitational Potential shape 
 * NONE for no potential.
 * */
#define GRAV_HOMOGENEOUS          1
#define GRAV_HERNQUIST            2
#define GRAV_HERNQUIST_NFW        3
#define GRAV_SINGLE_ISOTHERMAL    4
#define GRAV_DOUBLE_ISOTHERMAL    5


/* JD_MODE values.
 * Jet Domain modes: 
 * JD_GRAD is based on pressure gradient (originally in PLUTO) 
 * JD_PRES is based on the absolute value of pressure. 
 *
 * JD_CONST is a constant appearing in the comparison value (see
 * jet_domain.c, GetRightmostIndex. A good value in the case of JD_PRES is 
 * 1.1, a good value in the case of JD_GRAD is 1.e-5, although it could be 
 * higher depending on the background profile and grid resolution. */
#define JD_GRAD 0
#define JD_PRES 1



/* #######################################################
              include definitions_usr.h 
              and macros.h here 
   */
#include "definitions_usr.h"
#include "macros_usr.h"
/* ####################################################### */



/* Defaults for constants  that need a value */

/* Maximum timestep */
#ifndef DTMAX
#define DTMAX NONE
#endif

/* MU_CALC the method of calculating mu */
#ifndef MU_CALC
#define MU_CALC MU_CONST
#endif

/* MU_NORM . A default value for constant mean
 * molecular weight. From latest Mappings model of an ionized ISM
 * this is approximately 0.60364 */
#ifndef MU_NORM
#define MU_NORM 0.60364
#endif


#ifndef NOZZLE
#define NOZZLE NOZZLE_JET
#endif


#ifndef NOZZLE_FILL
#define NOZZLE_FILL NF_PRIMITIVE
#endif


/* Use a hemispherical cap?
 * By default, if INTERNAL_BOUNDARY is off, a hemispherical
 * cap is used, and this cannot be overridden. However, in the
 * case of INTERNAL_BOUNDARY == YES, there is the choice of
 * whether to use a hemispherical cap to buffer the jet inlet,
 * or just use the spherical section above the cone from the
 * NozzleSphere as a buffer. */
#ifndef NOZZLE_CAP
#define NOZZLE_CAP YES
#endif


/* Accretion */
#ifndef ACCRETION
#define ACCRETION NO
#endif

#ifndef ACCRETION_OUTPUT
#define ACCRETION_OUTPUT NO
#endif

#ifndef ACCRETION_OUTPUT_RATE
#define ACCRETION_OUTPUT_RATE 0
#endif

#ifndef SIC_METHOD
#define SIC_METHOD SIC_RADIUS
#endif

#ifndef SID_METHOD
#define SID_METHOD SID_REGIONS
#endif

/* Smoothing parameter of BH potential
 * Eqn (18), Ruffert (1994) */
#ifndef BH_POT_SMOOTH
#define BH_POT_SMOOTH 4.0
#endif

/* Sink */
#ifndef SINK_METHOD
#define SINK_METHOD SINK_FREEFLOW
#endif

/* Measure Bondi Accretion? */
#ifndef MEASURE_BONDI_ACCRETION
#define MEASURE_BONDI_ACCRETION  NO
#endif
#if SINK_METHOD == SINK_BONDI
#define MEASURE_BONDI_ACCRETION  YES
#endif
#if SINK_METHOD == SINK_FEDERRATH
#undef MEASURE_BONDI_ACCRETION
#define MEASURE_BONDI_ACCRETION  NO
#endif




/* Feedback Cycle */

/* Switch */

#ifndef FEEDBACK_CYCLE
#define FEEDBACK_CYCLE FALSE
#endif

/* Feedback cycle modes */

#ifndef FBC_DEBOOST_MODE
#define FBC_DEBOOST_MODE FBC_DEBOOST_MODE_1
#endif


/* Supernovae */

#ifndef NSTARS_MAX
#define NSTARS_MAX 1000
#endif

#ifndef SUPERNOVAE
#define SUPERNOVAE NO
#endif

#if SUPERNOVAE == YES

#ifndef SUPERNOVA_ENERGY
#define SUPERNOVA_ENERGY  1.e51
#endif

#ifndef SUPERNOVA_RCELLS
#define SUPERNOVA_RCELLS  3
#endif

#endif


/* Clouds */

#ifndef CLOUDS
#define CLOUDS NO
#endif

#ifndef CLOUD_DENSITY
#define CLOUD_DENSITY CD_HOMOGENEOUS
#endif

#ifndef CLOUD_SCALE
#define CLOUD_SCALE  CS_SCALE_HEIGHT
#endif

#ifndef CLOUD_VELOCITY
#define CLOUD_VELOCITY NONE
#endif

#ifndef CUBE_ENDIANNESS
#define CUBE_ENDIANNESS "little"
#endif

#ifndef CLOUD_EXTRACT
#define CLOUD_EXTRACT NONE
#endif

/* A factor with which to underpressure clouds w.r.t. ambient medium (<1)*/
#ifndef CLOUD_UNDERPRESSURE
#define CLOUD_UNDERPRESSURE 0.98
#endif

/* A critical temperature for thermal instability (K) */
#ifndef CLOUD_TCRIT
#define CLOUD_TCRIT 3.e4
#endif

#ifndef CLOUD_MUCRIT
#define CLOUD_MUCRIT  0.6212407755077543
#endif

/* Use setup that can place multiple individual clouds into the grid */
#ifndef CLOUDS_MULTI
#define CLOUDS_MULTI NO
#endif

/* Turn off default cloud init if clouds_multi = yes to prevent double initialisation.
   Comment out lines below if both are needed for a set up.
*/
#if CLOUDS_MULTI == YES
#define CLOUDS NO
#endif



/* Parameters that accompany a particular potential */
#ifndef GRAV_POTENTIAL
#define GRAV_POTENTIAL NONE
#endif


/* Whether the gravity potential requires table */
#if (GRAV_POTENTIAL == GRAV_HERNQUIST_NFW) || \
    (GRAV_POTENTIAL == GRAV_SINGLE_ISOTHERMAL) || \
    (GRAV_POTENTIAL == GRAV_DOUBLE_ISOTHERMAL)

#define GRAV_TABLE

/* Default gravity and hot halo filename */
#ifndef GRAV_FNAME
#define GRAV_FNAME            "gravtable.dat"
#endif

#ifndef HOT_FNAME
#define HOT_FNAME             "hottable.dat"
#endif

#endif


#ifndef JD_MODE
#define JD_MODE JD_GRAD
#endif

#if JD_MODE == JD_GRAD
  #ifndef JD_CONST
  #define JD_CONST 1.e-5
  #endif
#elif JD_MODE == JD_GRAD
  #ifndef JD_CONST
  #define JD_CONST 1.1
  #endif
#endif



/* For further refinement modes define these variables in TagCells.cpp */
#define NPLUS_REF_VARS 0
#define REF_VAR2  (TRC)
#define REF_VAR3  ((TRC) + 1)


/* -- GLOBAL VARIABLES -- */

/* Generally, global variables are defined in the file they are 
 * first used, and their extern declaration is in the header of 
 * the file. But in cases where the variable is used very often 
 * or in cases where there isn't a header, e.g. because it's 
 * just a little edit of a source code file, include the extern 
 * declaration here. */

/* Implies Chombo is being used */
#ifdef CH_SPACEDIM
extern int maxLevel;
#endif

/* Extents of input data. Declared in input_data.c */
extern double g_idBoxBeg[3];  /**< Lower limits of the input data domain. */
extern double g_idBoxEnd[3];  /**< Upper limits of the input data domain. */
extern double g_idnx1, g_idnx2, g_idnx3;
extern int g_idnvar, *g_idvarindx;
extern int last_step;

#endif



/* -- PROTOTYPING -- */

/* Generally, prototypes are in the header of the file they are
 * first used. But in cases where the variable is used very often
 * or in cases where there isn't a header, e.g. because it's just
 * a little edit of a source code file, include the prototype here. */