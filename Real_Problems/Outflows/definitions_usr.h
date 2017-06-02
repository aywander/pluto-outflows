/* -- Own definitions choices -- */

#define NOZZLE                            NOZZLE_UFO
#define NOZZLE_FILL                       NF_PRIMITIVE
#define NOZZLE_CAP                        YES

#define ACCRETION                         YES
#define ACCRETION_OUTPUT                  YES
#define ACCRETION_OUTPUT_RATE             0.1
#define SIC_METHOD                        SIC_HYBRID
#define SID_METHOD                        SID_REGIONS
#define SINK_METHOD                       SINK_FREEFLOW
#define MEASURE_BONDI_ACCRETION           YES
#define TURBULENT_BONDI_ACCRETION         YES
#define FEEDBACK_CYCLE                    NO
#define FBC_DEBOOST_MODE                  FBC_DEBOOST_MODE_3

#define GRAV_POTENTIAL                    GRAV_DOUBLE_ISOTHERMAL

#define CLOUDS                            YES
#define CLOUD_REPEAT	                  NO
#define CLOUDS_MULTI   	                  NO

#define CLOUD_DENSITY                     CD_KEPLERIAN
#define CLOUD_VELOCITY                    NONE
#define CLOUD_SCALE                       CS_VELOCITY_DISPERSION
#define CLOUD_EXTRACT_ELLIPSOID           FALSE
#define CLOUD_EXTRACT_CENTRAL_BUFFER      FALSE

#define COORDINATE_SYSTEM_DEBUG           FALSE

/* --- Not usually changed ---- */
#define MU_CALC                           MU_ANALYTIC
#define CLOUD_TCRIT                       3.0e4
#define CLOUD_MUCRIT                      0.6212407755077543
#define JD_MODE                           JD_GRAD
#define BH_POT_SMOOTH                     4.0


