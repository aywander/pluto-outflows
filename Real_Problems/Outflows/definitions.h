#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 TABULATED
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 2
#define  USER_DEF_PARAMETERS     17
#define  USER_DEF_CONSTANTS      4

/* -- physics dependent declarations -- */

#define    EOS                     IDEAL
#define    ENTROPY_SWITCH          YES
#define    THERMAL_CONDUCTION      NO
#define    VISCOSITY               NO
#define    ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  PAR_OPOW               0
#define  PAR_OANG               1
#define  PAR_OSPD               2
#define  PAR_OMDT               3
#define  PAR_ORAD               4
#define  PAR_OTHK               5
#define  PAR_ODIR               6
#define  PAR_HRHO               7
#define  PAR_HTE                8
#define  PAR_HVBG               9
#define  PAR_HRAD               10
#define  PAR_WRHO               11
#define  PAR_WTRB               12
#define  PAR_WRAD               13
#define  PAR_WROT               14
#define  PAR_WSMF               15
#define  PAR_LEV1               16

/* -- user-defined symbolic constants -- */

#define MU_NORM 0.6165
#define UNIT_DENSITY      (CONST_amu)*(MU_NORM)
#define UNIT_LENGTH       1.e3*(CONST_pc)
#define UNIT_VELOCITY     (UNIT_LENGTH)/(1.e3*(CONST_ly)/(CONST_c))

/* -- supplementary constants (user editable) -- */ 
#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      NO
#define  PRINT_TO_FILE         YES
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      MULTID
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         YES
#define  LIMITER               DEFAULT

