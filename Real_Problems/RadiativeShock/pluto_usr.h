/* COOLING_OFIT <int>
 * COOLING_NFIT <int>
 * Advanced power-law cooling parameters for direct integration of internal energy equation. 
 * By using POWER_LAW cooling and setting these we can use
 *  - piece-wise power-law fits to a cooling curve: COOLING_OFIT = 1
 *    which is a fit to log Lambda(log T/Tcrit)
 *  - piece-wise polynomial fits to a cooling curve: COOLING_OFIT > 1
 *    which is a fit to 1./Lambda(log T/Tcrit)
 * The two cases are treated separately in cooling.c, PowerLawCooling. 
 * COOLING_PL_NFIT is the number of piecewise segments. */


/* Include user definitions here (last) */
#include "definitions_usr.h"

