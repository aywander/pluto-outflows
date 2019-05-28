#define  PHYSICS                 HD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              POTENTIAL
#define  FORCED_TURB             NO
#define  COOLING                 NO
#define  RECONSTRUCTION          PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     4

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          SELECTIVE

/* -- user-defined parameters (labels) -- */

#define  PAR_MACH             0
#define  PAR_DENS             1
#define  PAR_DRATIO           2
#define  PAR_RCORE            3

/* [Beg] user-defined constants (do not change this line) */

#define  MU_NORM                 0.60364
#define  UNIT_DENSITY            1.
#define  UNIT_LENGTH             1.
#define  UNIT_VELOCITY           1.

#define  ARTIFICIAL_VISC         NO
//#define  ASSIGN_VECTOR_POTENTIAL YES
#define  CHAR_LIMITING           YES
//#define  CT_EN_CORRECTION        YES
#define  ID_NZ_MAX               4
#define  INITIAL_SMOOTHING       NO
#define  INTERNAL_BOUNDARY       YES
#define  LIMITER                 MC_LIM
#define  RECONSTRUCT_4VEL        NO
#define  SHOCK_FLATTENING        MULTID
#define  PPM_ORDER               4
//#define  UPDATE_VECTOR_POTENTIAL YES
#define  WARNING_MESSAGES        NO

//#define  FORCED_TURB_ENERGY      2.e-3
//#define  FORCED_TURB_DECAY       0.5
//#define  FORCED_TURB_KMIN        (2 * CONST_PI)
//#define  FORCED_TURB_KMAX        (6 * CONST_PI)
//#define  FORCED_TURB_WEIGHT      0.3

#define  SHOW_TIME_STEPS         YES
#define  SHOW_TIMING             YES

/* [End] user-defined constants (do not change this line) */

