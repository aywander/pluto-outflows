//
// Created by Alexander Y. Wagner on 2018/03/23.
//

#ifndef PLUTO_FRICTION_H
#define PLUTO_FRICTION_H

#include "pluto.h"

/* Structure for dynamical friction. All are kept in code units. */
typedef struct {
    double fx1;               // total gas-dynamical friction in x1 dir
    double fx2;               // total gas-dynamical friction in x2 dir
    double fx3;               // total gas-dynamical friction in x3 dir

    int nr;                   // number of r bins
    double rmin;              // minimum r
    double rmax;              // maximum r
    double *r;                // radius array
    double *dr;               // radius array bin widths
    double dy;                // radius array bin widths (in log space)
    double *frx1;             // force as a function of radius
    double *frx2;             // force as a function of radius
    double *frx3;             // force as a function of radius

    int nth;                  // number of theta bins
    double thmin;             // minimum theta
    double thmax;             // maximum theta
    double *th;               // theta array
    double dth;               // theta array
    double *fthx1;            // force as a function of theta
    double *fthx2;            // force as a function of theta
    double *fthx3;            // force as a function of theta

} GasdynamicalFriction;

extern GasdynamicalFriction gf;

void SetFriction(const Grid *grid);

void FrictionForce(const Data *d, Grid *grid);

void FrictionForceRadius(const Data *d, Grid *grid);

void FrictionForcePolar(const Data *d, Grid *grid);

void FrictionForceSlab(const Data *d, Grid *grid);

void FrictionForceOutput();

void FrictionForceRadiusOutput();

void FrictionForcePolarOutput();

void FrictionForceSlabOutput();

void OptimalRadialBinning(const Grid *grid, double *rmin, double *rmax, int *nr);

void OptimalPolarBinning(const Grid *grid, double *thmin, double *thmax, int *nth);

void MakeLogRadiusArray(double *r, double *dr, double *dy, double rmin, double rmax, int nr);

void MakeUniformPolarArray(double *th, double *dth, double thmin, double thmax, int nth);

#endif //PLUTO_FRICTION_H

