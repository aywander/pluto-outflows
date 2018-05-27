/* -- Own definitions choices -- */

#define NOZZLE                            NOZZLE_UFO
#define NOZZLE_FILL                       NF_PRIMITIVE
#define NOZZLE_CAP                        YES
#define NOZZLE_DT                         1.e-5

#define ACCRETION                         NO
#define SINK_METHOD                       SINK_FREEFLOW

#define FEEDBACK_CYCLE                    NO
#define FBC_DEBOOST                       NO
#define FBC_DEBOOST_MODE                  FBC_DEBOOST_MODE_3

#define GRAV_POTENTIAL                    NONE

#define CLOUDS                            NO
#define CLOUD_REPEAT	                  NO
#define CLOUDS_MULTI   	                  NO

#define CLOUD_DENSITY                     NONE
#define CLOUD_VELOCITY                    NONE
#define CLOUD_SCALE                       CS_VELOCITY_DISPERSION
#define CLOUD_EXTRACT_ELLIPSOID           YES
#define CLOUD_EXTRACT_CENTRAL_BUFFER      YES

#define ACCRETION_OUTPUT                  NO
#define ACCRETION_OUTPUT_RATE             1.0
#define BONDI_ACCRETION_OUTPUT            NO
#define TURBULENT_BONDI_ACCRETION         NO
#define OUTFLOW_OUTPUT                    NO
#define OUTFLOW_OUTPUT_RATE               1.0
#define CLOUD_OUTPUT                      NO
#define CLOUD_OUTPUT_RATE                 1.0

#define COORDINATE_SYSTEM_DEBUG           FALSE


/* --- Not usually changed ---- */
#define MU_CALC                           MU_ANALYTIC
#define CLOUD_TCRIT                       3.0e4
#define CLOUD_MUCRIT                      0.6212407755077543
#define JD_MODE                           JD_GRAD
#define BH_POT_SMOOTH                     4.0


