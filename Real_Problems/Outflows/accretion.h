//
// Created by Alexander Y. Wagner on 4/6/16.
//

#ifndef PLUTO_ACCRETION_H
#define PLUTO_ACCRETION_H

#include "init_tools.h"

/* Structure for accretion physics. All are kept in code units. */
typedef struct {
    double rad;               // Radius of surface through which accretion rate is measured
    double mbh;               // BH mass
    double edd;               // Eddington power
    double eff;               // Fraction of accretion rate going into outflow
    double snk;               // Sink radius
    double accr_rate;         // Accretion rate (global)
    double area;              // Area of surface
    double rho_acc;           // denisty far away
    double prs_acc;           // pressure far away
    double snd_acc;           // sound speed far away
    double accr_rate_bondi;   // sound speed far away
} Accretion;

extern Accretion ac;


void SetAccretionPhysics();

int InSinkRegion(const double x1, const double x2, const double x3);

void SphericalFreeflow(double *prims, double ****VC, const double *x1, const double *x2, const double *x3,
                       const int k, const int j, const int i);

void SphericalFreeflowInternalBoundary(const double ****Vc, int i, int j, int k, const double *x1, const double *x2,
                                       const double *x3, double *result);

double BondiAccretionRate(const double mbh, const double rho_far, const double snd_far);

double BondiAccretionRateLocal(const double mbh, const double rho_acc, const double snd_acc, const double snd_far);

double BondiLambda();

void BondiFlowInternalBoundary(const double x1, const double x2, const double x3, double *result);


void VacuumInternalBoundary(double *result);



double EddingtonLuminosity(const double mbh);

void SphericalAccretion(const Data *d, Grid *grid);

void SphericalAccretionOutput();

#endif //PLUTO_ACCRETION_H
