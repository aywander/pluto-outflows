#define  PHYSICS                 HD
#define  DIMENSIONS              1
#define  COMPONENTS              1
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  FORCED_TURB             NO
#define  COOLING                 TABULATED
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     8

/* -- physics dependent declarations -- */

#define    EOS                     IDEAL
#define    ENTROPY_SWITCH          NO
#define    THERMAL_CONDUCTION      NO
#define    VISCOSITY               NO
#define    ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  PAR_MACH           0
#define  PAR_SSPD           1
#define  PAR_DENS           2
#define  PAR_TEMP           3
#define  PAR_BPRP           4
#define  PAR_BPAR           5
#define  PAR_SLOC           6
#define  PAR_TMIN           7

/* [Beg] user-defined constants (do not change this line) */

// CIE OS 2013
#define  MU_NORM                 0.618

// NEQ 2015
//#define  MU_NORM                 0.60364

#define  UNIT_DENSITY            (CONST_amu * MU_NORM)

// v300-t1e6-cie
#define  UNIT_LENGTH             (69.1284408863 * 1.e3 * CONST_pc)
#define  UNIT_VELOCITY           (208.3426991365 * 1.e5)

// v200-t1e6-cie
//#define  UNIT_LENGTH             (2.0762906886 * 1.e3 * CONST_pc)
//#define  UNIT_VELOCITY           (98.1195372780 * 1.e5)

// v140-t3e5-cie
//#define  UNIT_LENGTH             (0.0639909744 * 1.e3 * CONST_pc)
//#define  UNIT_VELOCITY           (39.7094589211 * 1.e5)

// v200-t1e6-neq
//#define  UNIT_LENGTH             (11.7126525214 * 1.e3 * CONST_pc)
//#define  UNIT_VELOCITY           (98.1195372780 * 1.e5)

// v140-t3e5-neq
//#define  UNIT_LENGTH             (0.4847402503 * 1.e3 * CONST_pc)
//#define  UNIT_VELOCITY           (39.7094589211 * 1.e5)

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
