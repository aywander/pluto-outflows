//
// Created by Alexander Y. Wagner on 2017/06/29.
//

#ifndef PLUTO_SFR_H
#define PLUTO_SFR_H

void FederrathAccretion(Data *d, const Grid *grid);

double JeansResolvedDensity(const double *prim);

double FederrathSinkInternalBoundary(const double ****Vc, int i, int j, int k, const double *x1, const double *x2,
                                     const double *x3, const double ***vol, double *result);

double VirialParameter(const double * prim, const double mass,
                       const double x1, const double x2, const double x3);

double GravitationallyBound(const double *prim, const double mass, const double vol,
                            const double x1, const double x2, const double x3);

double StarFormationRateDensity(Data *d, Grid *grid, int i, int j, int k);

double StarFormationEfficiency(Data *d, Grid *grid, int i, int j, int k);

double StarFormationCriteria(Data *d, Grid *grid, int i, int j, int k);

//#include "accretion.c"

#endif //PLUTO_SFR_H
