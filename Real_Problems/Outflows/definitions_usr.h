/* -- Own definitions choices -- */

#undef DTMAX

#define NOZZLE               NOZZLE_UFO
#define NOZZLE_SHAPE         ANN

#define CLOUDS               SINGLE_CUBE
#define CLOUD_VELOCITY       NO
#define CUBE_BYTESWAP        NO
#define TCRIT                3.0e4
#define CLOUD_EXTRACT        CE_ELLIPSOID
#define CLOUD_DISTR          TURB_KEPLERIAN_DISC

#define GRAV_POTENTIAL       HQNFW

#define HOT_DISTR            HOT_EXT_DATA

#define JD_MODE              JD_PRES
#define JD_CONST             1.1

#define CCL_OK               NO

#define MU_CALC              MU_CONST
