#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              VECTOR
#define  COOLING                 TABULATED
#define  RECONSTRUCTION          PARABOLIC
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 2
#define  USER_DEF_PARAMETERS     39
#define  USER_DEF_CONSTANTS      4

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          SELECTIVE
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO
#define  RECONSTRUCT_4VEL        NO

/* -- user-defined parameters (labels) -- */

#define  PAR_OPOW                0
#define  PAR_OSPD                1
#define  PAR_OMDT                2
#define  PAR_OANG                3
#define  PAR_ORAD                4
#define  PAR_ODBH                5
#define  PAR_OSPH                6
#define  PAR_ODIR                7
#define  PAR_OOMG                8
#define  PAR_OPHI                9
#define  PAR_ARAD                10
#define  PAR_AMBH                11
#define  PAR_AEFF                12
#define  PAR_ASNK                13
#define  PAR_HRHO                14
#define  PAR_HTMP                15
#define  PAR_HVX1                16
#define  PAR_HVX2                17
#define  PAR_HVX3                18
#define  PAR_HVRD                19
#define  PAR_HRAD                20
#define  PAR_WRHO                21
#define  PAR_WTRB                22
#define  PAR_WRAD                23
#define  PAR_WROT                24
#define  PAR_WX1L                25
#define  PAR_WX1H                26
#define  PAR_WX2L                27
#define  PAR_WX2H                28
#define  PAR_WX3L                29
#define  PAR_WX3H                30
#define  PAR_WVRD                31
#define  PAR_WVPL                32
#define  PAR_WVPP                33
#define  PAR_WVAN                34
#define  PAR_SGAV                35
#define  PAR_NCLD                36
#define  PAR_LOMX                37
#define  PAR_LCMX                38

/* [Beg] user-defined constants (do not change this line) */

#define  MU_NORM                 0.60364
#define  UNIT_DENSITY            ((CONST_amu) * (MU_NORM))
#define  UNIT_LENGTH             ((CONST_pc) * 1.e3)
#define  UNIT_VELOCITY           CONST_c

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING      NO
#define  WARNING_MESSAGES       NO
#define  PRINT_TO_FILE          YES
#define  INTERNAL_BOUNDARY      YES
#define  SHOCK_FLATTENING       MULTID
#define  ARTIFICIAL_VISC        NO
#define  LIMITER                MC_LIM
#define  SHOW_TIME_STEPS        YES
#define  CHAR_LIMITING          NO
#define  PPM_ORDER              5
