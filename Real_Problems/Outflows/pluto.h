/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO main header file.

  Contains basic macro definitions, structure definitions and global 
  variable declarations used by the code.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#ifndef PLUTO_H
#define PLUTO_H

#define PLUTO_VERSION  "4.0"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*! Return the maximum between two numbers. */
#define MAX(a,b)  ( (a) >= (b) ? (a) : (b) ) 

/*! Return the minimum between two numbers. */
#define MIN(a,b)  ( (a) <= (b) ? (a) : (b) ) 
      
/*! Return the number with the smaller absolute value. */
#define ABS_MIN(a,b)  (fabs(a) < fabs(b) ? (a):(b)) 
                         
/*! Return the sign of x. */
#define DSIGN(x)      ( (x) >= 0.0 ? (1.0) : (-1.0))

#define MINMOD(a,b)  ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
#define VAN_LEER(a,b) ((a)*(b) > 0.0 ? 2.0*(a)*(b)/((a)+(b)):0.0)
#define MC(a,b) (MINMOD(0.5*((a)+(b)), 2.0*MINMOD((a),(b))))
#define SWAP_VAR(x) SwapEndian(&x, sizeof(x));

#define YES      1
#define NO       0
#define DEFAULT -1

         /*  ----  GEOMETRY LABELS ( > 0) ----  */

#define CARTESIAN    1
#define AXISYM       1   /*  should be the same as CARTESIAN in order 
                             to avoid complications                    */
#define CYLINDRICAL  2
#define POLAR        3
#define SPHERICAL    4

/*! If this macro is defined, use the geometrical cell-center
    in place of the cell-center in CYLINDRICAL geometry. 
    The cell-center used to be the default in PLUTO 3.1 but lead to 
    spurious features on the symmetry axis.  */
#define NEW_GEOM     1

#define UNIFORM_GRID             1
#define STRETCHED_GRID           2
#define LOGARITHMIC_INC_GRID     3
#define LOGARITHMIC_DEC_GRID     4

    /*  ----  SET LABELS FOR EQUATION OF STATE (EOS)  ----  */

#define IDEAL         1
#define TAUB          3
#define BAROTROPIC    4
#define ISOTHERMAL    5

         /*  ----  TIME_STEPPING LABELS  ----  */

#define EULER                      1
#define HANCOCK                    2
#define CHARACTERISTIC_TRACING     3
#define RK2                        5
#define RK3                        6
#define RK_MIDPOINT                7

#define EXPLICIT             1 /* -- just a number different from 0 !!!  -- */
#define SUPER_TIME_STEPPING  2 /* -- just a number different from EXPLICIT -- */ 
#define RK_CHEBYSHEV         4  

     /*  ----   Output labels:  ---- */

#define DBL_OUTPUT   1
#define FLT_OUTPUT   2
#define VTK_OUTPUT   3
#define DBL_H5_OUTPUT   4
#define FLT_H5_OUTPUT   5
#define TAB_OUTPUT   6
#define PPM_OUTPUT   7
#define PNG_OUTPUT   8

#define VTK_VECTOR  5  /* -- any number but NOT 1  -- */

/*! The maximum number of output file formats is fixed to 11 so that the 
    size of runtime structure (defined below) is 64 bytes. 
    This should prevent, for some compilers, attempts to change the 
    alignment of the structure and therefore troubleshooting when restarting 
    from files written on different architectures.                    */
#define MAX_OUTPUT_TYPES 11 

   /*  ----  SET LABELS FOR COOLING  ----  */

#define POWER_LAW    3
#define MINEq        4
#define SNEq         5
#define TABULATED    6
#define H2_COOL      7

  /*  ----  SET LABELS FOR PHYSICS MODULE  ----  */

#define HD     1
#define RHD    2
#define MHD    3
#define RMHD   5

    /*  ----  SET LABELS FOR DIV.B REMOVAL  ----  
        If you move them to the MHD header, 
        definitions.h (which is included before)
        cannot correctly use them                */
        
#define NONE                   0
#define EIGHT_WAVES            1
#define DIV_CLEANING           2
#define CONSTRAINED_TRANSPORT  3

   /*  ----  SET LABELS FOR BODY_FORCE  ----
      Please do not change them since they are
      used in bitwise operations                */
   
#define VECTOR     4   /* corresponds to  100 in binary  */
#define POTENTIAL  8   /* corresponds to 1000 in binary  */

   /*  ----  SET LABELS FOR BOUNDARIES  ----  */

#define OUTFLOW          1  /* any number except 0 !! */
#define REFLECTIVE       2 
#define AXISYMMETRIC     3
#define EQTSYMMETRIC     4
#define PERIODIC         5
#define SHEARING         6
#define USERDEF          7

#define X1_BEG           101
#define X1_END           102
#define X2_BEG           103
#define X2_END           104
#define X3_BEG           105
#define X3_END           106

    /* ---- LABELS FOR IMAGE SLICING ---- */

#define X12_PLANE       3  
#define X13_PLANE       5
#define X23_PLANE       6

/*! \name Bit flag labels.
    The following macros define the bits that can be turned on or off
    in an unsigned char (= 1 byte = 8 bits) variable. 
    Different bit flags allow to enable or disable certain actions in 
    a given cell at different points in the code, see also flag.c.
    The 3D unsigned char \c ***flag array is used for bookeeping, in each zone 
    (i,j,k), which bits are actually switched on or off.
    A simple bitwise operation is used to enable a flag, e.g., 
    <tt> flag[k][j][i] |= FLAG_XXX </tt>.
    For instance, by turning the ::FLAG_HLL bit on, we have
    <tt> flag = 00000100 </tt>, while by also enabling the ::FLAG_SPLIT_CELL
    one has <tt> flag = 00010100 </tt> and so on.
*/
/**@{ */
#define FLAG_MINMOD      1  /**< Reconstruct using MINMOD limiter. */
#define FLAG_FLAT        2  /**< Reconstruct using FLAT limiter.   */
#define FLAG_HLL         4  /**< Use HLL Riemann solver. */
#define FLAG_ENTROPY     8  /**< Update pressure using entropy equation. */
#define FLAG_SPLIT_CELL  16 /**< Zone is covered by a finer level (AMR only). */
#define FLAG_INTERNAL_BOUNDARY   32  /**< Zone belongs to an internal boundary
                                          region and should be excluded from 
                                          being updated in time              */
#define FLAG_BIT7        64
#define FLAG_BIT8       128  
/**@} */

#define PRS_FAIL  1
#define ENG_FAIL  2
#define RHO_FAIL  4

#define IDIR    0     /*   This sequence (0,1,2) should */
#define JDIR    1     /*   never be changed             */
#define KDIR    2     /*                                */
#define ALL_DIR -1

/* -- location of a variable inside the cell -- */

#define CENTER  0
#define X1FACE  1
#define X2FACE  2
#define X3FACE  3

#define CELL_CENTER    50  /* really needed ? */
#define FACE_CENTER    51
#define EDGE_CENTER    52

        /*  ----  SET LABELS FOR INTERPOLATION  ----   */

#define FLAT              1
#define LINEAR            2
#define CENO3             3
#define PARABOLIC         4
#define LINEAR_MULTID     5
#define MP5               6
#define LimO3             7
#define WENO3             8

#define WENO3_FD             103
#define WENO5_FD             105
#define WENOZ_FD             106
#define WENO7_FD             107
#define MP5_FD               125
#define PPM_FD               140
#define LIMO3_FD             300

#define ONED   1
#define MULTID 3

/* ---- limiter labels ---- */

#define FLAT_LIM          1
#define MINMOD_LIM        2 
#define VANALBADA_LIM     3
#define UMIST_LIM         4
#define VANLEER_LIM       5
#define MC_LIM            6
#define FOURTH_ORDER_LIM  7

/*

#define FLAT_LIMIT(dv, dp, dm) dv = 0.0;

#define MINMOD_LIMIT(dv, dp, dm)  dv = ((dm)*(dp) > 0.0 ? (fabs(dm) < fabs(dp) ? (dm):(dp)):0.0)

#define VANALBADA_LIMIT(dv, dp, dm)\
  if (dp*dm > 0.0) { \
    double _dpp= dp*dp, _dmm = dm_dm; \
    dv = dp*(_dmm + 1.e-18) + dm*(_dpp + 1.e-18))/(_dpp + _dmm + 1.e-18); \
  }else dv = 0.0;

#define VANLEER_LIMIT(dv, dp, dm)  if (dp*dm > 0.0) dv = 2.0*dp*dm/(dp+dm);\
                                 else             dv = 0.0;

#define UMIST_LIMIT(dv,dp,dm)\
  if (dp*dm > 0.0){ \
  double _ddp = 0.25*(dp + 3.0*dm), _ddm = 0.25*(dm + 3.0*dp); \
  double _d2  = 2.0*(fabs(dp) < fabs(dm) ? dp:dm);  \
  _d2 = (fabs(_d2) < fabs(_ddp) ? _d2:_ddp); \
  dv  = (fabs(_d2) < fabs(_ddm) ? _d2:_ddm);} else dv = 0.0;

#define MC_LIMIT(dv, dp, dm, dc)  \
  if (dp*dm > 0.0){  \
    double _d2 = 2.0*(fabs(dp) < fabs(dm)) ? dp:dm); \
    dv = 2.0*(fabs(dvp) < fabs(dvm) ? dvp:dvm); \
  }else dv = 0.0;
*/
 

/* ###############################################################

                 PROBLEM-DEPENDENT DECLARATIONS

   ############################################################### */

#include "definitions.h"

#ifdef PARALLEL
 #include <al.h>
#endif

/* -- include mpi.h for parallel Chombo, in order to use
      MPI_Abort function in the QUIT_PLUTO macro -- */
 
#ifdef CH_MPI
 #include <mpi.h>
#endif


/* ###################################################################

        SWITCHES USED FOR DEBUGGING PURPOSES ONLY

   ################################################################### */
 
/* -- CHECK_EIGENVECTORS: used in eigenv.c in HD/, MHD/, RHD/
      to check orthogonality and the correctness through 
      the relation the A = L*\Lambda*R  -- */

#define CHECK_EIGENVECTORS      NO

/* -- CHECK_ROE_MATRIX: used in HD/roe.c, MHD/roe.c to 
      check that the characteristic decomp. reproduces
      Roe matrix -- */

#define CHECK_ROE_MATRIX        NO

/* -- CHECK_CONSERVATIVE_VAR: used in RHD/mappers.c to 
      check that conservative vars are physical -- */

#define CHECK_CONSERVATIVE_VAR  NO

/* -- CHECK_DIVB_CONDITION: used in MHD/CT/ct.c 
      to check if div.B = 0 -- */

#ifndef CHECK_DIVB_CONDITION
 #define CHECK_DIVB_CONDITION   NO
#endif

/* ###############################################################

                 Define precision (useless now)
   
   ############################################################### */

typedef double real;

/* ###############################################################
    
     EXIT MACRO. For Chombo it is defined elsewhere.

   ############################################################### */

#ifdef PARALLEL
/* #define QUIT_PLUTO(e_code)   {MPI_Finalize(); exit(e_code);}  */
 #define QUIT_PLUTO(e_code)   \
        {MPI_Abort(MPI_COMM_WORLD, e_code);MPI_Finalize(); exit(e_code);}
#elif (defined CH_MPI)
 #define QUIT_PLUTO(e_code)   \
        {MPI_Abort(MPI_COMM_WORLD, e_code); exit(e_code);}
#else
 #define QUIT_PLUTO(e_code)   exit(e_code);
#endif
        
/* ##################################################################

           MACROS TO EXPAND DIMENSION-DEPENDENT LINES

   #################################################################  */

/*! \def EXPAND(a,b,c)
    Allows to write dimension-independent code by expanding as many arguments
    as the value of COMPONENTS. 
    The result is that only the first argument will be expanded in 1D, 
    the first two arguments in 2D and all of them in 3D. */

/*! \def D_EXPAND(a,b,c)
    Similar to the EXPAND() macro but the expansion depends on DIMENSIONS.  */

/*! \def SELECT(a,b,c)
    Expand only the 1st, 2nd or 3rd argument based on the value of
    COMPONENTS.                                                       */

/*! \def D_SELECT(a,b,c)
    Expand only the 1st, 2nd or 3rd argument based on the value of
    DIMENSIONS.                                                       */

#if COMPONENTS == 1
 #define EXPAND(a,b,c) a
 #define SELECT(a,b,c) a
#endif

#if COMPONENTS == 2
 #define EXPAND(a,b,c) a b
 #define SELECT(a,b,c) b
#endif

#if COMPONENTS == 3
 #define EXPAND(a,b,c) a b c
 #define SELECT(a,b,c) c
#endif

#if DIMENSIONS == 1
 #define D_EXPAND(a,b,c)  a
 #define D_SELECT(a,b,c)  a
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 0
 #define KOFFSET 0
#endif

#if DIMENSIONS == 2
 #define D_EXPAND(a,b,c) a b
 #define D_SELECT(a,b,c) b
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 1
 #define KOFFSET 0
#endif

#if DIMENSIONS == 3
 #define D_EXPAND(a,b,c) a b c
 #define D_SELECT(a,b,c) c
 #define IOFFSET 1   /* ** so we know when to loop in boundary zones ** */
 #define JOFFSET 1
 #define KOFFSET 1
#endif

#if WARNING_MESSAGES == YES
 #define WARNING(a)  a
#else
 #define WARNING(a)
#endif

/* ################################################################# 

                       SPATIAL LOOP MACROS        

   ################################################################# */

#define VAR_LOOP(n)   for ((n) = NVAR; (n)--;    )
#define DIM_LOOP(d)   for ((d) = 0; (d) < DIMENSIONS; (d)++)

/*! \name Spatial loop macros.
    The following macros provide a compact way to perform 1D or multi-D
    loops in selected regions of the (local) computational domain.
    "BEG" and "END" refers to the leftmost or rightmost boundary region in 
    the corresponding direction I, J or K.
    "DOM" is used for looping inside the computational domain (boundaries
    excluded) while "TOT" is used for looping inside the ghost zones.
*/
/**@{ */
#define IBEG_LOOP(i)  for ((i) = IBEG; (i)--;    )
#define JBEG_LOOP(j)  for ((j) = JBEG; (j)--;    )
#define KBEG_LOOP(k)  for ((k) = KBEG; (k)--;    )

#define IEND_LOOP(i)  for ((i) = IEND + 1; (i) < NX1_TOT; (i)++)
#define JEND_LOOP(j)  for ((j) = JEND + 1; (j) < NX2_TOT; (j)++)
#define KEND_LOOP(k)  for ((k) = KEND + 1; (k) < NX3_TOT; (k)++)

#define IDOM_LOOP(i)  for ((i) = IBEG; (i) <= IEND; (i)++)
#define JDOM_LOOP(j)  for ((j) = JBEG; (j) <= JEND; (j)++)
#define KDOM_LOOP(k)  for ((k) = KBEG; (k) <= KEND; (k)++)

#define ITOT_LOOP(i)  for ((i) = 0; (i) < NX1_TOT; (i)++)
#define JTOT_LOOP(j)  for ((j) = 0; (j) < NX2_TOT; (j)++)
#define KTOT_LOOP(k)  for ((k) = 0; (k) < NX3_TOT; (k)++)

#define DOM_LOOP(k,j,i) KDOM_LOOP(k) JDOM_LOOP(j) IDOM_LOOP(i)

#define TOT_LOOP(k,j,i) KTOT_LOOP(k) JTOT_LOOP(j) ITOT_LOOP(i)

#define X1_BEG_LOOP(k,j,i) KTOT_LOOP(k) JTOT_LOOP(j) IBEG_LOOP(i)
#define X2_BEG_LOOP(k,j,i) KTOT_LOOP(k) JBEG_LOOP(j) ITOT_LOOP(i)
#define X3_BEG_LOOP(k,j,i) KBEG_LOOP(k) JTOT_LOOP(j) ITOT_LOOP(i)

#define X1_END_LOOP(k,j,i) KTOT_LOOP(k) JTOT_LOOP(j) IEND_LOOP(i)
#define X2_END_LOOP(k,j,i) KTOT_LOOP(k) JEND_LOOP(j) ITOT_LOOP(i)
#define X3_END_LOOP(k,j,i) KEND_LOOP(k) JTOT_LOOP(j) ITOT_LOOP(i)

#define TRANSVERSE_LOOP(indx, in, i,j,k) \
 if (g_dir == IDIR) {i = &in; j = &indx.t1; k = &indx.t2;} \
 if (g_dir == JDIR) {j = &in; i = &indx.t1; k = &indx.t2;} \
 if (g_dir == KDIR) {k = &in; i = &indx.t1; j = &indx.t2;} \
 g_i = i; g_j = j; g_k = k;\
 for (indx.t2 = indx.t2_beg; indx.t2 <= indx.t2_end; indx.t2++) \
 for (indx.t1 = indx.t1_beg; indx.t1 <= indx.t1_end; indx.t1++)
/**@} */

/* ********************************************************************* */
/*! The BOX_LOOP() macro implements a loop over (i,j,k) in a rectangular 
    portion of the domain with indices defined by the RBox structure B. 
    The loop increments (di,dj,dk) are members of the structure which 
    are here initialized to either 1 or -1 depending on whether the 
    lower corner index lies below or above the upper index 
    (e.g. B->ib <= B->ie or not). 
   ********************************************************************* */
#define BOX_LOOP(B,k,j,i) \
 for ((B)->dk = ((k=(B)->kb) <= (B)->ke ? 1:-1); k != (B)->ke+(B)->dk; k += (B)->dk)\
 for ((B)->dj = ((j=(B)->jb) <= (B)->je ? 1:-1); j != (B)->je+(B)->dj; j += (B)->dj)\
 for ((B)->di = ((i=(B)->ib) <= (B)->ie ? 1:-1); i != (B)->ie+(B)->di; i += (B)->di)


/*
#define DIR_LOOP(grid,in,k,j,i) \
    if (g_dir == IDIR) {in = &i; j = *g_j; k = *g_k;} \
    if (g_dir == JDIR) {in = &j; i = *g_i; k = *g_k;} \
    if (g_dir == KDIR) {in = &k; i = *g_i; j = *g_j;} \
    for (*in = grid[g_dir].lbeg - 2; *in <= grid[g_dir].lend + 2; (*in)++)
*/

/*  not used for now... will serve in EMF_BOUNDARY...
#define ISTAG_LOOP(i)  for ((i) = IBEG - 1; (i) <= IEND; (i)++)
#define JSTAG_LOOP(j)  for ((j) = JBEG - 1; (j) <= JEND; (j)++)
#define KSTAG_LOOP(k)  for ((k) = KBEG - 1; (k) <= KEND; (k)++)
*/


/* ################################################################# 

                   VARIABLE LOOP MACROS        

   ################################################################# */

 /* ------------------------------------------------------------
     The LRVAR_LOOP and LRFLX_LOOP will speed up computation 
     of the left/right states at a given interface by excluding 
     the normal component of magnetic field, which is usually  
     assigned at the end of the reconstruction step.
    ------------------------------------------------------------ */
/*      
#ifdef STAGGERED_MHD
 #define LR_VAR_LOOP(i) for (i = NVAR; i--; i -= (i == BXn+1)) 
 #define LR_FLX_LOOP(i) for (i = NFLX; i--; i -= (i == BXn+1)) 
#else
 #define LR_VAR_LOOP(i) for (i = NVAR; i--; ) 
 #define LR_FLX_LOOP(i) for (i = NFLX; i--; ) 
#endif
*/
/* ################################################################# 

                  TIME INTEGRATOR DEFINITIONS

   ################################################################# */

#if (TIME_STEPPING == HANCOCK) || (TIME_STEPPING == CHARACTERISTIC_TRACING)
 #define SINGLE_STEP   1
 #if DIMENSIONAL_SPLITTING == NO
  #define CTU      1    /* -- Corner Transport Upwind method of Colella -- */
 #endif
#endif

/* ------------------------------------------------------------
    the GET_MAX_DT switch determines how the time step is
    computed. For pure advection, setting it to YES will force
    PLUTO to compute dt in the old way, i.e., by taking the
    maximum.
    When set to NO, the time step is computed only during the
    predictor step and, for UNSPLIT RK schemes, it will be
    calculated as the average over dimensions resulting in 
    slightly larger time increments.
   ------------------------------------------------------------ */

#if ((TIME_STEPPING == RK2) || (TIME_STEPPING == RK3) || \
     (TIME_STEPPING == RK_MIDPOINT)) && DIMENSIONAL_SPLITTING == NO

 #define GET_MAX_DT    NO
#else
 #define GET_MAX_DT    YES
#endif

/* ############################################################### 
  
          Check whether a separate source step has to be done

   ############################################################### */

#if (COOLING != NO) 
 #define INCLUDE_SPLIT_SOURCE   YES
#else
 #define INCLUDE_SPLIT_SOURCE   NO
#endif


/* #################################################################

    set to NO all undefined identifier that may appear in the code.
    (compile with -Wundef).
    Diffusion operators: PARABOLIC_FLUX is the bitwise OR
    combination of all operators, each being either one of 
    NO, EXPLICIT (1st bit), STS (2nd bit). 
    It can take the following values

      00   --> no diffusion operator is being used
      01   --> there's at least one explicit diffusion operator and
               no sts.
               
      10   --> there's at least one sts diffusion operator and
               no explicit one.
      11   --> mixed: there is at least one explicit and sts operator
     
   ################################################################# */

#ifndef RESISTIVE_MHD 
 #define RESISTIVE_MHD   NO
#endif

#ifndef THERMAL_CONDUCTION
 #define THERMAL_CONDUCTION NO
#endif

#ifndef VISCOSITY
 #define VISCOSITY NO
#endif

#ifndef INCLUDE_PARTICLES
 #define INCLUDE_PARTICLES NO
#endif

#ifndef UPDATE_VECTOR_POTENTIAL
 #define UPDATE_VECTOR_POTENTIAL  NO
#endif

#ifndef BACKGROUND_FIELD
 #define BACKGROUND_FIELD NO
#endif

#ifndef USE_FOUR_VELOCITY
 #define USE_FOUR_VELOCITY NO
#endif

#ifndef CHAR_LIMITING
 #define CHAR_LIMITING  NO
#endif

#ifndef ARTIFICIAL_VISCOSITY 
 #define ARTIFICIAL_VISCOSITY  NO
#endif

#ifndef ROTATING_FRAME
 #define ROTATING_FRAME NO
#endif

#ifndef INTERNAL_BOUNDARY
 #define INTERNAL_BOUNDARY NO
#endif

#define PARABOLIC_FLUX (RESISTIVE_MHD|THERMAL_CONDUCTION|VISCOSITY)

/* ################################################################# 

                       USEFUL CONSTANTS        

   ################################################################# */

#define ACCURACY    1.e-14          /* ZERO ACCURACY */

/*! \name Constant physical value macros.
     The following set of macros express some useful physical constants
     in c.g.s units (erg, cm and sec). Values have been taken from
     http://physic.nist.gov/cuu/Constants/index.html
*/
/**@{ */
#define CONST_PI      3.14159265358979   /**<  \f$ \pi \f$.               */
#define CONST_amu     1.66053886e-24     /**<  Atomic mass unit.          */
#define CONST_mp      1.67262171e-24     /**<  Proton mass.               */
#define CONST_mn      1.67492728e-24     /**<  Neutron mass.              */
#define CONST_me      9.1093826e-28      /**<  Electron mass.             */
#define CONST_mH      1.6733e-24         /**<  Hydrogen atom mass.        */
#define CONST_kB      1.3806505e-16      /**<  Boltzmann constant.        */
#define CONST_sigma   5.67051e-5         /**<  Stephan Boltmann constant. */
#define CONST_sigmaT  6.6524e-25         /**<  Thomson Cross section.    */
#define CONST_NA      6.0221367e23       /**<  Avogadro Contant.          */
#define CONST_c       2.99792458e10      /**<  Speed of Light.            */
#define CONST_Msun    2.e33              /**<  Solar Mass.                */
#define CONST_Rsun    6.96e10            /**<  Solar Radius.              */
#define CONST_Mearth  5.9736e27          /**<  Earth Mass.                */
#define CONST_Rearth  6.378136e8         /**<  Earth Radius.              */
#define CONST_G       6.6726e-8          /**<  Gravitational Constant.    */
#define CONST_h       6.62606876e-27     /**<  Planck Constant.           */
#define CONST_pc      3.0856775807e18    /**<  Parsec.                    */
#define CONST_ly      0.9461e18          /**<  Light year.                */
#define CONST_au      1.49597892e13      /**<  Astronomical unit.         */
#define CONST_eV      1.602176463158e-12 /**<  Electron Volt in erg.      */
/**@} */

/* ############################################################### 
 
                     Define Structures here     

   ############################################################### */

typedef struct CMD_LINE {
  int restart;       
  int h5restart;       
  int nrestart;      
  int makegrid;      
  int write;         
  int maxsteps;
  int jet;  /* -- follow jet evolution in a given direction -- */
  int parallel_dim[3];
  int nproc[3];  /* -- user supplied number of processors -- */
  int show_dec; /* -- show domain decomposition ? -- */
  int xres; /* -- change the resolution via command line -- */
  /* AYW -- 2012-06-19 00:41 JST */
  float maxtime;
  int fill[15]; /* useless, it makes the struct a power of 2 */ 
  //int fill; /* useless, it makes the struct a power of 2 */ 
  /* --AYW */
} Cmd_Line;

/* ********************************************************************* */
/*! The Data structure contains the main solution 3D arrays used by 
    the code. 
   ********************************************************************* */
typedef struct DATA{
  double ****Uc;
  double ****Vc;  /**< The main four-index data array used for cell-centered
                       primitive variables. The index order is Vc[nv][k][j][i],
                       where nv gives the variable index, k,j and i are the
                       locations of the cell in the \f$x_3\f$,
                       \f$x_2\f$ and \f$x_1\f$ direction. */
  double ****Vs;  /**< The main four-index data array used for face-centered
                       staggered magnetic fields. 
                       The index order is Vc[nv][k][j][i],
                       where nv gives the variable index, k,j and i are the
                       locations of the cell in the \f$x_3\f$,
                       \f$x_2\f$ and \f$x_1\f$ direction. */
  double ****Vuser; /**< Array storing user-defined supplementary variables 
                         written to disk. */ 
  double ***Ax1; /**< The vector potential component in the \f$x_1\f$ direction.*/
  double ***Ax2; /**< The vector potential component in the \f$x_2\f$ direction.*/
  double ***Ax3; /**< The vector potential component in the \f$x_3\f$ direction.*/
  unsigned char ***flag; /**< Pointer to a 3D array setting useful integration
                              flags that are retrieved during integration. */
  int fill; /* useless, it makes the struct a power of 2 */ 
} Data;
   
/* ********************************************************************* */
/*! This structure contains one-dimensional vectors of conserved variables,   
    primitive variables, fluxes and so on, used during the 
    reconstruct-Solve-Average strategy. 
    It is a frequently passed to the Riemann solver routines, source and 
    flux functions, etc.
   ********************************************************************* */
typedef struct STATE_1D{
  double **v;    /**< Cell-centered primitive varables at the base time level,
                      v[i] = \f$ \vec{V}^n_i \f$ . */
  double **vL;   /**< Primitive variables to the left of the interface, 
                       \f${\rm vL[i]} \equiv \vec{V}_{i,+} = 
                           \vec{V}^L_{i+\HALF} \f$. */
  double **vR;   /**< Primitive variables to the right of the interface, 
                       \f$\mathrm{vR[i]} \equiv \vec{V}^R_{i+\HALF} \f$. */
  double **vm;   /**< prim vars at i-1/2 edge, vm[i] = vR(i-1/2)     */
  double **vp;   /**< prim vars at i+1/2 edge, vp[i] = vL(i+1/2)     */

  double **uL;   /**< same as vL, in conservative vars */
  double **uR;   /**< same as vR, in conservative vars */
  double **um;   /**< same as vm, in conservative vars */
  double **up;   /**< same as vp, in conservative vars */

  double **flux;    /**< upwind flux computed with the Riemann solver */
  double **par_flx; /**< parabolic fluxes (viscosity + thermal cond. +
                          resistivity) */
  double **pnt_flx;
  double **dff_flx;

  double ***Lp, ***Rp; /**< Left and right primitive eigenvectors */
  double **lambda;     /**< Characteristic speed associated to Lp and Rp */
  double *lmax;   /**< Define the maximum k-characteristic speed over the domain */
  double *a2;     /**< Sound speed squared */ 
  double *h;      /**< Enthalpy. */
  double **vt;
  double **src;     
  double **par_src; /**< Geometrical source terms for parabolic operators in
                         curvilinear coords */
  double **vh;      /**< Primitive    state at n+1/2 (only for one step method) */
  double **rhs;     /**< Conservative right hand side */
  double *press;    /**< Upwind pressure term computed with the Riemann solver */
  double *bn;       /**< Face magentic field, bn = bx(i+1/2) */
  double *SL;       /**< Leftmost  velocity in the Riemann fan at i+1/2 */
  double *SR;       /**< Rightmost velocity in the Riemann fan at i+1/2 */
  double q1;            /**< Useless, just to make struct a power of 2 */
  unsigned char uc1;    /**< Useless, just to make struct a power of 2 */
  unsigned char *flag;
} State_1D;

/* ********************************************************************* */
/*! The PLUTO Grid structure contains information pertaining to the 
    computational mesh in a specific 1D coordinate direction. 
    Since PLUTO assumes a logically rectangular system of coordinates, 
    the whole computational domain is obtained as the cartesian product
    of 2 or 3 grid structures.\n

    In parallel, each processor owns a different portion of the domain
    and the grid structures will be different.
    For this reason, in the following member description, we use the 
    word "global" or "local" to refer the the whole computational domain
    or to the sub-domain owned by a single processor.
 
    Similarly, variables ending with a "glob" suffix are intended to be 
    global, i.e., they refer to the whole computational stencil and 
    not to the local processor sub-domain.
   ********************************************************************* */
typedef struct GRID{
  double xi, xf;        /**< Leftmost and rightmost point in the local domain. */
  double *x, *x_glob;   /**< Cell geometrical central points. */
  double *xr, *xr_glob; /**< Cell right interface. */
  double *xl, *xl_glob; /**< Cell left interface. */
  double *dx, *dx_glob; /**< Cell size.  */ 
  double *xgc;          /**< Cell volumetric centroid 
                             (!= x when geometry != CARTESIAN).  */
  double *dV;           /**< Cell volume.  */
  double *A;            /**< Right interface area, A[i] = \f$A_{i+\HALF}\f$. */
  double *r_1;          /**< Geometrical factor 1/r.  */
  double *ct;           /**< Geometrical factor cot(theta).  */
  double *inv_dx;       /**<      */
  double *inv_dxi;      /**< inverse of the distance between the center of 
                             two cells, inv_dxi = \f$\DS \frac{2}{\Delta x_i +
                             \Delta x_{i+1}}\f$.     */
  double dl_min;      /**<  minimum cell length (e.g. min[dr, r*dth,
                           r*sin(th)*dphi] (GLOBAL DOMAIN).  */

  double *dfg;  /**< Coefficient weight for computing forward differences. */
  double *df2g; /**< Coefficient weight for computing central differences. */
  double *dfR;  /**< Interpolation weight for right interface state. */
  double *dfL;  /**< Interpolation weight for left  interface state. */
 
  int np_tot_glob; /**< Total number of points in the global domain 
                        (boundaries included). */
  int np_int_glob; /**< Total number of points in the global domain 
                        (boundaries excluded). */
  int np_tot;      /**< Total number of points in the local domain 
                        (boundaries included). */
  int np_int;      /**< Total number of points in the local domain 
                        (boundaries excluded). */
  int nghost;      /**< Number of ghost zones. */
  int lbound;      /**< When different from zero, it specifies the boundary
                        condition to be applied at leftmost grid side where  
                        the physical boundary is located.
                        Otherwise, it equals zero if the current 
                        processor does not touch the leftmost physical boundary. 
                        This evantuality (lbound = 0) is possible only
                        in PARALLEL mode.  */ 
  int rbound;      /**< Same as lbound, but for the right edge of the grid. */
  int gbeg;        /**< Global start index for the global array. */
  int gend;        /**< Global end   index for the global array. */
  int beg;         /**< Global start index for the local array. */
  int end;         /**< Global end   index for the local array. */
  int lbeg;        /**< Local start  index for the local array. */
  int lend;        /**< Local end    index for the local array. */
  int uniform;     /* = 1 when the grid is cartesian AND uniform everywhere  */
  int nproc;       /**< number of processors for this grid. */
  int rank_coord;  /**< Parallel coordinate in a Cartesian topology. */
  int level;       /**< The current refinement level (chombo only). */
  char fill[84];   /* useless, just to make the structure size a power of 2 */
} Grid;

/* ********************************************************************* */
/*! The Time_Step structure contains essential information for 
    determining the time step.
   ********************************************************************* */

typedef struct TIME_STEP{
  double *cmax;     /**< Maximum signal velocity for hyperbolic eqns. */
  double inv_dta;   /**< Inverse of advection (hyperbolic) time step, 
                         \f$ \lambda/\Delta l\f$.*/
  double inv_dtp;   /**< Inverse of diffusion (parabolic)  time step 
                         \f$ \eta/\Delta l^2\f$. */
  double dt_cool;   /**< Cooling time step. */
  double cfl;       /**< Courant number for advection. */
  double cfl_par;   /**< Courant number for diffusion (STS only). */
  double rmax_par;  
  int    Nsts;      /**< Maximum number of substeps used in STS. */
  int    Nrkc;      /**< Maximum number of substeps used in RKC. */
  char  fill[24];   /* useless, just to make the structure size a power of 2 */
} Time_Step;


/* ********************************************************************* */
/*! The Output structure contains essential information for I/O.
   ********************************************************************* */
typedef struct OUTPUT{
  int    type;       /**< output format (DBL, FLT, ...) - one per output */
  int    nvar;       /**< tot. # of vars that can be written - same for all   */
  int    user_outs;  
  int    nfile;      /**< current number being saved - one per output */
  int    dn;         /**< step increment between outputs    - one per output */
  int    *stag_var;  /**< centered or staggered variable    - same for all   */
  int    *dump_var;  /**< select vars being written         - one per output */
  char   mode[32];   /**< single or multiple files          - one per output */
  char   **var_name; /**< variable names                    - same for all   */
  char   ext[8];     /**< output extension                                   */
  double dt;         /**< time increment between outputs    - one per output */
  double dclock;     /**< time increment in clock hours     - one per output */
  double ***V[64];   /**< pointer to arrays being written   - same for all  */
  char   fill[168];  /**< useless, just to make the structure size a power of 2 */
} Output;

typedef struct INPUT{
  int    npoint[3];           /**< Global number of zones in the interior. */
  int    lft_bound_side[3];       /* left  boundary type */
  int    rgt_bound_side[3];      /* right boundary type */
  int    grid_is_uniform[3];    /* = 1 when grid is uniform, 0 otherwise */
  int    npatch[5];               /* number of grid patches  */
  int    patch_npoint[5][16];     /* number of points per patch */
  int    patch_type[5][16];             
  int    log_freq;              /* log frequency */
  int    user_var;              /* number of additional user-variables being
                                   held in memory and written to disk */
  int    anl_dn;                /*  number of step increment for ANALYSIS */
  char   solv_type[64];
  char   user_var_name[128][128];
  Output output[MAX_OUTPUT_TYPES];  
  double patch_left_node[5][16];  /*  self-expl. */
  double  xbeg[3];
  double  xend[3];
  double  cfl;                    /* hyperbolic cfl number */
  double  cfl_max_var;
  double  cfl_par;           /* (STS) parabolic  cfl number */
  double  rmax_par;          /* (STS) max ratio between current time
                                step and parabolic time step */
  double  tstop;
  double  first_dt;
  double  anl_dt;          /* time step increment for ANALYSIS */
  double  aux[32];         /* we keep aux inside this structure, 
                              since in parallel execution it has
                              to be comunicated to all processors  */
} Input;

typedef struct RUNTIME{
  int nstep;
  int nfile[MAX_OUTPUT_TYPES];
  double t;
  double dt;
} Runtime;

typedef struct RGB{
  unsigned char r, g, b;
} RGB;

typedef struct IMAGE{
  int    nrow, ncol;    /* -- image rows and columns -- */
  int    slice_plane;   /* -- one of X12_PLANE, X13_PLANE, X23_PLANE -- */
  int    logscale;      /* -- YES/NO for log scale -- */
  char   *colormap;     /* -- colormap name -- */
  char   basename[32];  /* -- image base name (no extensions) -- */
  unsigned char r[256], g[256], b[256]; /* -- colortable saved here -- */
  RGB    **rgb;         /* -- rgb array containing image values -- */
  double max;           /* -- max image value -- */
  double min;           /* -- min image value -- */
  double slice_coord;   /* -- slice coord orthogonal to slice_plane -- */
} Image;

typedef struct FLOAT_VECT{
  float v1, v2, v3;
} Float_Vect;

typedef struct INDEX{
  int ntot, beg, end;
  int t1, t1_beg, t1_end;
  int t2, t2_beg, t2_end;
  char fill[28]; /* useless, just to make the structure size a power of 2 */
} Index;

/* ********************************************************************* */
/*! The RBox (= Rectangular Box) defines a rectangular portion of the 
    domain in terms of the grid indices [ib,jb,kb] corresponding to 
    the lower corner and [ie,je,ke] corresponding to the upper corner. 
    The integer vpos specifies the variable location with respect to 
    the grid (e.g. center/staggered). 

    \note The lower and upper grid indices may also be reversed 
          (e.g. box->ib > box->ie). In this case the macro ::BOX_LOOP
          automatically reset the directional increment (box->di) to -1.
   ********************************************************************* */
typedef struct RBOX{
  int ib; /**< Lower corner index in the x1 direction. */
  int ie; /**< Upper corner index in the x1 direction. */
  int jb; /**< Lower corner index in the x2 direction. */
  int je; /**< Upper corner index in the x2 direction. */
  int kb; /**< Lower corner index in the x3 direction. */
  int ke; /**< Upper corner index in the x3 direction. */
  int di; /**< Directional increment (+1 or -1) when looping over the 1st 
               dimension of the box. Automatically set by the ::BOX_LOOP macro. */
  int dj; /**< Directional increment (+1 or -1) when looping over the 2nd 
               dimension of the box. Automatically set by the ::BOX_LOOP macro. */
  int dk; /**< Directional increment (+1 or -1) when looping over the 3rd 
               dimension of the box. Automatically set by the ::BOX_LOOP macro. */
  int vpos; /**< Location of the variable inside the cell. */
} RBox;

/* ####################################################### 

            Recurrent function types 
 
   ####################################################### */

typedef void Riemann_Solver (const State_1D *, int, int, double *, Grid *);
typedef void Limiter        (double *, double *, double *, int, int, Grid *);
typedef double Reconstruct  (double *, double, int);
typedef double ****Data_Arr;

/* -- when using Finite Difference Schemes, the "Riemann Solver" 
      function computes the fluxes with high order interpolants. -- */

/* ----------------------------------------------------------
             Include module header files
   ---------------------------------------------------------- */

#if   PHYSICS == HD

 #include "HD/mod_defs.h"

#elif PHYSICS == RHD

 #include "RHD/mod_defs.h"

#elif PHYSICS == MHD

 #include "MHD/mod_defs.h"
 #ifdef SHEARINGBOX
  #include "MHD/ShearingBox/shearingbox.h"
 #endif

#elif PHYSICS == RMHD

 #include "RMHD/mod_defs.h"

#endif

#if COOLING == MINEq
 #include "Cooling/MINEq/cooling.h"
#elif COOLING == SNEq 
 #include "Cooling/SNEq/cooling.h"
#elif COOLING == POWER_LAW
 #include "Cooling/Power_Law/cooling.h"
#elif COOLING == TABULATED
 #include "Cooling/Tab/cooling.h"
#elif COOLING == H2_COOL
 #include "Cooling/H2_COOL/cooling.h"
#endif

#if THERMAL_CONDUCTION != NO
 #include "Thermal_Conduction/tc.h"
#endif

#if VISCOSITY != NO
 #include "Viscosity/viscosity.h"
#endif

#if INCLUDE_PARTICLES == YES
 #include "Particles/particles.h"
#endif

#ifdef FARGO
 #include "Fargo/fargo.h"
#endif

/* *****************************************************
   *****************************************************

                  GLOBAL VARIABLES

   *****************************************************
   ***************************************************** */

 extern int SZ;
 extern int SZ_stagx;
 extern int SZ_stagy;
 extern int SZ_stagz;
 extern int SZ_char;
 extern int SZ_float;
 extern int SZ_Float_Vect;
 extern int SZ_rgb;
 extern int SZ_short;
 extern int prank;

extern long int IBEG, IEND, JBEG, JEND, KBEG, KEND;
extern long int NX1, NX2, NX3;
extern long int NX1_TOT, NX2_TOT, NX3_TOT;

extern long int NMAX_POINT;

extern int VXn, VXt, VXb;
extern int MXn, MXt, MXb;
extern int BXn, BXt, BXb;

extern int *g_i, *g_j, *g_k;

extern int g_dir;
extern int g_maxRiemannIter;
extern long int g_usedMem;
extern long int g_stepNumber;
extern int      g_intStage;

extern double g_unitDensity, g_unitLength, g_unitVelocity;

extern double g_maxCoolingRate, g_minCoolingTemp;

extern double g_smallDensity, g_smallPressure;

extern double g_time, g_dt;
extern double g_maxMach;
#if ROTATING_FRAME
 extern double g_OmegaZ;
#endif

extern double g_domBeg[3], g_domEnd[3];

extern double g_inputParam[32];
#if EOS != ISOTHERMAL
 extern double g_gamma;
#elif EOS == ISOTHERMAL
 extern double g_isoSoundSpeed;
#endif

#ifdef CH_SPACEDIM
 extern double ch_max, ch_max_loc, coeff_dl_min;
 extern int use_glm;
 extern int LEVEL;
#endif

extern FILE *pluto_log_file;

/* ----------------------------------------------------------
    define NX_MAX, NY_MAX, NZ_MAX as the maximum 
    size on a patch.
    For the static grid of PLUTO, they simply coincide
    with NX1_TOT, NX2_TOT, NX3_TOT. Otherwise
    the grid is taken to be the maximum, according
    to maxGridSize. 
    This avoids allocating and freeing memory all the times.
   ---------------------------------------------------------- */

#ifdef CH_SPACEDIM
 #define NX_MAX NMAX_POINT
 #define NY_MAX NMAX_POINT
 #define NZ_MAX NMAX_POINT
#else
 #define NX_MAX NX1_TOT
 #define NY_MAX NX2_TOT
 #define NZ_MAX NX3_TOT
#endif 

/* ------------------------------------------------------------------
    Define total number of variables to be integrated;
    this includes:

     NFLX     = Number of equations defining the system of 
                conservation laws. For example, for the 
                HD module, it consists of density, momentum and energy.
                It is defined in the physics module header file mod_defs.h.

     NTRACER  = Number of user-defined tracers; defined in the problem
                directory header file definitions.h

     NIONS = Number of chemical species; it is defined in 
                the cooling modules cooling.h, if present.

  nv = 0...NFLX - 1; NFLX...NFLX+NIONS-1; TRC...TRC + NTRACER-1; ENTR; 
       <----------> <------------------> <------------------>  
          NFLX              NIONS               NTRACER              

                    <---------------------------------------------->
                                    NSCL

   ------------------------------------------------------------------ */

#ifndef ENTROPY_SWITCH
 #define ENTROPY_SWITCH  NO 
#endif 

#ifndef NIONS
 #define NIONS 0
#endif

#define TRC   (NFLX + NIONS)
#define TR   TRC
#define NSCL (NTRACER + NIONS + ENTROPY_SWITCH)

#if ENTROPY_SWITCH == YES
 #define ENTR  (TRC + NTRACER)
#else
 #if EOS != ISOTHERMAL && EOS != BAROTROPIC
  #define ENTR (ENG)
 #endif
#endif

#define NVAR (NFLX + NSCL)

/*! Define the conversion constant between dimensionless 
    temperature prs/rho and physical temperature in Kelvin,
    T = (prs/rho)*KELVIN*mu                                   */
#define KELVIN (g_unitVelocity*g_unitVelocity*CONST_amu/CONST_kB) 

#include "prototypes.h"
#endif /* PLUTO_H */
