/* -- Own definitions choices -- */

#define NOZZLE               NOZZLE_UFO

#define CLOUDS               YES
#define CLOUD_VELOCITY       NO
#define CLOUD_TCRIT          3.0e4
#define CLOUD_EXTRACT        CE_ELLIPSOID
#define CLOUD_DISTR          TURB_KEPLERIAN_DISC

#define GRAV_POTENTIAL       HQNFW

#define HOT_DISTR            HOT_EXT_DATA

#define JD_MODE              JD_GRAD

#define MU_CALC              MU_CONST

#define CCL_OK               TRUE

/* -- user-defined symbolic constants -- */
/* These are just here as a backup */
#if !(defined UNIT_DENSITY) || !(defined UNIT_DENSITY) || !(defined UNIT_DENSITY)
#define MU_NORM           0.6165
#define UNIT_DENSITY      (CONST_amu)*(MU_NORM)
#define UNIT_LENGTH       1.e3*(CONST_pc)
#define UNIT_VELOCITY     (UNIT_LENGTH)/(1.e3*(CONST_ly)/(CONST_c))
#endif
