/* -- Own definitions choices -- */

#define NOZZLE               NOZZLE_JET

#define CLOUDS               YES
#define CLOUD_VELOCITY       NO
#define CLOUD_TCRIT          3.0e4
#define CLOUD_EXTRACT        CE_HEMISPHERE_RIM
#define CLOUD_DISTR          CD_FLAT

#define GRAV_POTENTIAL       NONE

#define HOT_DISTR            HOT_FLAT

#define JD_MODE              JD_GRAD

#define MU_CALC              MU_CONST

#define CCL_OK               TRUE

/* -- user-defined symbolic constants -- */
#if !(defined UNIT_DENSITY) || !(defined UNIT_DENSITY) || !(defined UNIT_DENSITY)
#define MU_NORM           0.6165
#define UNIT_DENSITY      (CONST_amu)*(MU_NORM)
#define UNIT_LENGTH       1.e3*(CONST_pc)
#define UNIT_VELOCITY     (UNIT_LENGTH)/(1.e3*(CONST_ly)/(CONST_c))
#endif
