/* -- Own definitions choices -- */

#define NOZZLE                      NOZZLE_UFO
#define NOZZLE_FILL                 NF_PRIMITIVE
#define NOZZLE_CAP                  YES

#define ACCRETION                   NO
#define ACCRETION_OUTPUT            NO
#define ACCRETION_OUTPUT_RATE       0.15318627450980393
#define SIC_METHOD                  SIC_HYBRID
#define SID_METHOD                  SID_REGIONS
#define SINK_METHOD                 SINK_FEDERRATH
#define MEASURE_BONDI_ACCRETION     NO
#define FEEDBACK_CYCLE              NO
#define FBC_DEBOOST_MODE            FBC_DEBOOST_MODE_3

#define GRAV_POTENTIAL              GRAV_DOUBLE_ISOTHERMAL

#define CLOUDS                      YES
#define CLOUD_REPEAT	            NO
#define CLOUDS_MULTI   	            NO

#define CLOUD_DENSITY               CD_KEPLERIAN
#define CLOUD_VELOCITY              CD_KEPLERIAN
#define CLOUD_SCALE                 CS_VELOCITY_DISPERSION
#define CLOUD_EXTRACT               NONE


#define MU_CALC                     MU_ANALYTIC
/* --- Not usually changed ---- */
#define CLOUD_TCRIT                 3.0e4
#define CLOUD_MUCRIT                0.6212407755077543
#define JD_MODE                     JD_GRAD
#define BH_POT_SMOOTH               4.0


