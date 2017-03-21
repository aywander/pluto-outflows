/* -- Own definitions choices -- */

#define NOZZLE                      NOZZLE_UFO
#define NOZZLE_CAP                  YES

#define ACCRETION                   YES
#define ACCRETION_OUTPUT            YES
#define ACCRETION_OUTPUT_RATE       0.15318627450980393
#define SIC_METHOD                  SIC_HYBRID
#define SID_METHOD                  SID_REGIONS
#define SINK_METHOD                 SINK_FREEFLOW
#define MEASURE_BONDI_ACCRETION     YES
#define FEEDBACK_CYCLE              YES
#define FBC_DEBOOST_MODE            FBC_DEBOOST_MODE_3

#define SUPERNOVAE                  NO
#define SUPERNOVA_ENERGY            1.e51

#define GRAV_POTENTIAL              GRAV_HERNQUIST

#define CLOUDS                      YES
#define CLOUD_REPEAT	              NO
#define CLOUDS_MULTI   	            NO

#define CLOUD_DENSITY               CD_KEPLERIAN
//#define CLOUD_VELOCITY              CV_KEPLERIAN
//#define CLOUD_SCALE                 CS_VELOCITY_DISPERSION

#define MU_CALC                     MU_ANALYTIC

/* --- Not usually changed ---- */
#define CLOUD_TCRIT                 3.0e4
#define CLOUD_MUCRIT                0.6212407755077543
//#define CLOUD_MUCRIT                1.24
#define CLOUD_UNDERPRESSURE         0.98
#define CLOUD_EXTRACT               NONE
#define JD_MODE                     JD_GRAD
#define BH_POT_SMOOTH               4.0


