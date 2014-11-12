#include "pluto.h"
#include "pluto_usr.h"
#include "nrEOS.h"

/* ************************************************ 
 * Quick routines to get a specific variable
 * To get enthalpy, entropy, and sound speeds, use functions
 * Enthalpy, Entropy, and SoundSpeed2
 * Only Nonrelativistic EOS
 ************************************************** */


/* ************************************************ */
double PresNrEOS(const double dens, const double temp, const double mu) {
/*!
 * Return pressure from density and temperature
 * in code units
 *
 ************************************************** */

    return dens*temp/mu;
}



/* ************************************************ */
double TempNrEOS(const double dens, const double pres, const double mu) {
/* !
 * Return temperature from density and pressure
 * in code units
 *
 ************************************************** */

    return pres/dens*mu;
}


/* ************************************************ */
double DensNrEOS(const double pres, const double temp, const double mu) {
/* !
 * Return density from pressure and temperature
 * in code units
 *
 ************************************************** */

    return pres/temp*mu;
}




