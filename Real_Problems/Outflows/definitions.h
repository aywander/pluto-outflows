#define  PHYSICS                 RHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 TABULATED
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 2
#define  USER_DEF_PARAMETERS     19
#define  USER_DEF_CONSTANTS      4

/* -- physics dependent declarations -- */

#define    EOS                     TAUB
#define    ENTROPY_SWITCH          YES
#define    USE_FOUR_VELOCITY       YES

/* -- user-defined parameters (labels) -- */

#define  PAR_OPOW                0
#define  PAR_OSPD                1
#define  PAR_OMDT                2
#define  PAR_OANG                3
#define  PAR_ORAD                4
#define  PAR_ODBH                5
#define  PAR_ODIR                6
#define  PAR_OOMG                7
#define  PAR_OPHI                8
#define  PAR_HRHO                9
#define  PAR_HTMP               10
#define  PAR_HVBG               11
#define  PAR_HRAD               12
#define  PAR_WRHO               13
#define  PAR_WTRB               14
#define  PAR_WRAD               15
#define  PAR_WROT               16
#define  PAR_WSMF               17
#define  PAR_LEV1               18

/* -- user-defined symbolic constants -- */

#define MU_NORM           0.6165
#define UNIT_DENSITY      (CONST_amu)*(MU_NORM)
#define UNIT_LENGTH       1.e3*(CONST_pc)
#define UNIT_VELOCITY     CONST_c

/* -- supplementary constants (user editable) -- */ 
#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         YES
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      MULTID
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         YES
#define  LIMITER               DEFAULT

