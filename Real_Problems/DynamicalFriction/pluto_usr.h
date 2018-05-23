#ifndef PLUTO_USR_H
#define PLUTO_USR_H

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


/* Clouds */
/* Use a grid_in.out file to specify dimensions of cube
 * Input files are rho.dbl, vx1.dbl, vx2.dbl, vx3.dbl, etc*/


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


/* MU_NORM . A default value for constant mean
 * molecular weight. From latest Mappings model of an ionized ISM
 * this is approximately 0.60364 */
#ifndef MU_NORM
#define MU_NORM 0.60364
#endif


/* Accretion */
#ifndef ACCRETION
#define ACCRETION NO
#endif

#ifndef SIC_METHOD
#define SIC_METHOD SIC_RADIUS
#endif

#ifndef SID_METHOD
#define SID_METHOD SID_REGIONS
#endif

/* Smoothing parameter of BH potential
 * Eqn (18), Ruffert (1994) */
#ifndef BH_SCALE
#define BH_SCALE 4.0
#endif

/* Sink */
#ifndef SINK_METHOD
#define SINK_METHOD SINK_FREEFLOW
#endif


#if SINK_METHOD == SINK_BONDI
#define BONDI_ACCRETION_OUTPUT  YES
#endif




/* Clouds */

#ifndef CLOUDS
#define CLOUDS NO
#endif

#ifndef CUBE_ENDIANNESS
#define CUBE_ENDIANNESS "little"
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

#ifndef CLOUDS_VELOCITY
#define CLOUDS_VELOCITY NO
#endif

/* Turn off default cloud init if clouds_multi = yes to prevent double initialisation.
   Comment out lines below if both are needed for a set up.
*/
#if CLOUDS_MULTI == YES
#define CLOUDS NO
#endif


/* Jet domain */

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



/* Analaysis */

#ifndef ACCRETION_OUTPUT
#define ACCRETION_OUTPUT NO
#endif

#ifndef ACCRETION_OUTPUT_RATE
#define ACCRETION_OUTPUT_RATE 0
#endif

/* Measure Bondi Accretion? */
#ifndef BONDI_ACCRETION_OUTPUT
#define BONDI_ACCRETION_OUTPUT  NO
#endif

#ifndef TURBULENT_BONDI_ACCRETION
#define TURBULENT_BONDI_ACCRETION  NO
#endif

#ifndef CLOUD_OUTPUT
#define CLOUD_OUTPUT NO
#endif

#ifndef CLOUD_OUTPUT_RATE
#define CLOUD_OUTPUT_RATE 0
#endif

/* Debug */

/* Output in userdef_output for debugging of coordiante systems */
#ifndef DEBUG_USING_USERDEF_VARS
#define DEBUG_USING_USERDEF_VARS          NO
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