#define  PHYSICS                 HD
#define  DIMENSIONS              1
#define  COMPONENTS              1
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 TABULATED
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     8
#define  USER_DEF_CONSTANTS      4

/* -- physics dependent declarations -- */

#define    EOS                     IDEAL
#define    ENTROPY_SWITCH          NO
#define    THERMAL_CONDUCTION      NO
#define    VISCOSITY               NO
#define    ROTATING_FRAME          NO

/* -- pointers to user-def parameters -- */

#define  PAR_MACH           0
#define  PAR_SSPD           1
#define  PAR_DENS           2
#define  PAR_TEMP           3
#define  PAR_BPRP           4
#define  PAR_BPAR           5
#define  PAR_SLOC           6
#define  PAR_TMIN           7

/* -- user-defined symbolic constants -- */

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

/* -- supplementary constants (user editable) -- */

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      NO
#define  PRINT_TO_FILE         NO
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      ONED
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         YES
#define  LIMITER               DEFAULT
