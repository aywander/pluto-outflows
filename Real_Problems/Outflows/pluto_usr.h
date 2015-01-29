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

#define MUFNAME "mutable.dat"


/* Inlet */

/* NOZZLE values. 
 * Type of inlet nozzle 
 * This determines e.g. whether UfoPrimitives, 
 * JetPrimitives, etc are called.
 *   */
#define NOZZLE_JET      1
#define NOZZLE_UFO      2

/* Clouds */
/* Use a grid_in.out file to specify dimensions of cube 
 * Input files are rho.dbl, vx1.dbl, vx2.dbl, vx3.dbl, etc*/

/* CLOUD_EXTRACT values. 
 * Method of cloud extraction 
 * DEFAULT is NONE
 * BE is Bonnor-Ebert (not programmed yet)
 * Height of HEMISPHERE_RIM is determined by WRAD
 * Shape of ELLIPSOID is determined by parameter WROT and WRAD
 * */
#define CE_SPHERE_RIM        1
#define CE_BE                2
#define CE_DENS              3
#define CE_HEMISPHERE_RIM    4
#define CE_ELLIPSOID         5

/* CLOUD_DISTR values.
 * Cloud distribution 
 * The density profile for the warm phase */
#define CD_FLAT                     0
#define CD_TURB_ISOTH_HYDROSTATIC   1
#define CD_TURB_KEPLERIAN_DISC      2
#define CD_HERNQUIST                3

/* CLOUD_SCALE values.
 * How the scale height of clouds is set.
 * Default is by WRAD input parameter. */
#define CS_WRAD  0
#define CS_WTRB  1

/* CLOUD_VEL_DISTR values.
 * Cloud velocity distribution.
 * The velocity cube can be modified in different ways. */
#define CV_UNIFORM            0
#define CV_RADIAL             1
#define CV_KEPLERIAN_FRAC     2
#define CV_VIRIAL_FRAC        3


/* Static Gravity */

/* GRAV_POTENTIAL values. 
 * Gravitational Potential shape 
 * NONE for no potential.
 *
 * GRAVFNAME <string>
 * GRAVHNAME <string>
 * Gives the filename for the potential or
 * acceleration vector, and its header file.*/
#define HERNQUIST            1
#define SINGLE_ISO           2
#define DOUBLE_ISO           3
#define HQNFW                4
#define HOMOG                5
#define GRAVFNAME            "gravtable.dat"
#define GRAVHNAME            "gravtable.hdr"


/* HOT_DISTR values.
 * Hot phase distribution 
 * The density profile for the hot phase. */
#define HOT_FLAT 0
#define HOT_EXT_DATA 1
#define HOTFNAME "hottable.dat"
#define HOT_HYDROSTATIC_ISOTHERMAL 2
#define HOT_HERNQUIST 3
#define HOT_UNIFORM_FLOW 4  // Not used, but maybe useful?


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


/* CCL_OK {YES|NO}
 * Does CCL library work? E.g., on vayu with Intel compilers, it doesn't. */



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

#ifndef MU_CALC
#define MU_CALC MU_CONST
#endif

/* MU_NORM . A default value for constant mean 
 * molecular weight. From Mappings model of an ionized ISM with 
 * Asplund 2005 abundances, this is approximately 0.6156. */
#ifndef MU_NORM
#define MU_NORM 0.6165
#endif

#ifndef NOZZLE
#define NOZZLE NOZZLE_JET
#endif

/* Clouds */
/* NO for no clouds */
/* A faire
 * - Change this to TRUE or FALSE
 * - cp all files anew from Src, backup existing to *4.0.[ch], and vimdiff to port to 4.1
 * - Make sure there are no print1s anywhere.
 * */
#ifndef CLOUDS
#define CLOUDS NO
#endif

#ifndef CLOUD_VELOCITY
#define CLOUD_VELOCITY NO
#endif

#ifndef CUBE_ENDIANNESS
#define CUBE_ENDIANNESS "little"
#endif

#ifndef CLOUD_EXTRACT
#define CLOUD_EXTRACT NONE
#endif

#ifndef CLOUD_DISTR
#define CLOUD_DISTR CD_FLAT
#endif

#ifndef CLOUD_VEL_DISTR 
#define CLOUD_VEL_DISTR CV_UNIFORM
#endif

/* Each method is asssociated with a value of CV_VALUE 
 * which has a different meaning every time, but is usually a float.  */
#ifndef CV_VALUE
#define CV_VALUE 0.0
#endif

/* A factor with which to underpressure clouds w.r.t. ambient medium */
#ifndef CLOUD_UNDERPRESSURE
#define CLOUD_UNDERPRESSURE 0.98
#endif

/* A critical temperature for thermal instability (K) */
#ifndef CLOUD_TCRIT
#define CLOUD_TCRIT 3.e4
#endif


/* Parameters that accompany a particular potential */
#ifndef GRAV_POTENTIAL
#define GRAV_POTENTIAL NONE
#endif

#if GRAV_POTENTIAL == DOUBLE_ISO
  #define NPOTVAR 5
  #define KAP  0
  #define LAM  1
  #define RDM  2
  #define DHOT 3
  #define THOT 4

#elif GRAV_POTENTIAL == HQNFW
  #define NPOTVAR 6
  #define AHQ  0
  #define DHQ  1
  #define ANFW 2
  #define DNFW 3
  #define NHOT 4
  #define THOT 5

#endif 


/* Whether the gravity potential is of type table*/
#if (GRAV_POTENTIAL == SINGLE_ISO) || \
    (GRAV_POTENTIAL == DOUBLE_ISO) || \
    (GRAV_POTENTIAL == HQNFW)
  #define  GRAV_TABLE 
#endif


#ifndef HOT_DISTR
#define HOT_DISTR HOT_FLAT
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


/* CAn use CCL library? */
#ifndef CCL_OK
#define CCL_OK FALSE
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
//extern int g_idvarindx[ID_MAX_NVAR];

#endif

