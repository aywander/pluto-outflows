//
// Created by Alexander Y. Wagner on 4/6/16.
//

#ifndef PLUTO_ACCRETION_H
#define PLUTO_ACCRETION_H

/* Need the following for the reference nozzle and outflow states */
#include "outflow.h"
#include "nozzle.h"

/* Structure for accretion physics. All are kept in code units. */
typedef struct {
    double rad;               // Radius of surface through which accretion rate is measured
    double mbh;               // BH mass
    double edd;               // Eddington power
    double eff;               // Fraction of accretion rate going into outflow
    double mld;               // Mass loading efficiency
    double snk;               // Sink radius
    double area;              // Area of surface
    double accr_rate;         // Mass accretion rate, random spherical sampling
    double accr_rate_bondi;   // Bondi accretion rate
    double deboost;           // Deboost factor for outflow
    Nozzle nzi;               // Initial nozzle parameters
    OutflowState osi;         // Initial outflow state parameters

} Accretion;

extern Accretion ac;


void SetAccretionPhysics();

int InSinkRegion(const double x1, const double x2, const double x3);

void SphericalFreeflow(double *prims, double ****VC, const double *x1, const double *x2, const double *x3,
                       const int k, const int j, const int i);

void SphericalFreeflowInternalBoundary(double ****Vc, int i, int j, int k, const double *x1, const double *x2,
                                       const double *x3, double *result);


double BondiAccretionRate(const double mbh, const double rho_far, const double snd_far);

double BondiAccretionRateLocal(const double mbh, const double rho_acc, const double snd_acc, const double snd_far);

double BondiLambda();

double BondiRadius(double m, double *v);

void BondiFlowInternalBoundary(const double x1, const double x2, const double x3, double *result);


void VacuumInternalBoundary(double *result);

double EddingtonLuminosity(const double mbh);

double BondiAccretion(const Data *d, Grid *grid, const double radius);

void SphericalAccretion(const Data *d, Grid *grid);

double SphericalSampledAccretion(const Data *d, Grid *grid, const double radius);

double SphericalSelectedAccretion(const Data *d, Grid *grid, const double radius);

void BondiAccretionOutput();

void SphericalAccretionOutput();

#endif //PLUTO_ACCRETION_H
