#include "pluto.h"
#include "pluto_usr.h"
#include "idealEOS.h"

/* ************************************************ 
  Quick routines to get a specific variable
  To get enthalpy, entropy, and sound speeds, use 
  functions Enthalpy, Entropy, and SoundSpeed2.
  Only Nonrelativistic EOS.
      The most general ways to define eos s are 
  through the enthalpy (which is computed from 
  primitives) and then calculating entropy, 
  soundspeed. Temperature is ill-defined for 
  relativistic fluids.
 ************************************************** */


/* ************************************************ */
double PresIdealEOS(const double dens, const double temp, const double mu) {
/*!
 * Return pressure from density and temperature
 * in code units
 *
 ************************************************** */

    return dens * temp / mu;
}



/* ************************************************ */
double TempIdealEOS(const double dens, const double pres, const double mu) {
/* !
 * Return temperature from density and pressure
 * in code units
 *
 ************************************************** */

    return pres / dens * mu;
}


/* ************************************************ */
double DensIdealEOS(const double pres, const double temp, const double mu) {
/* !
 * Return density from pressure and temperature
 * in code units
 *
 ************************************************** */

    return pres / temp * mu;
}



/* ************************************************ */
double SoundSpeed2IdealEOS(const double dens, const double pres) {
/* !
 * Return Adiabatic Sound speed  from pressure and temperature
 * in code units
 *
 ************************************************** */

    return g_gamma * pres / dens;
}


