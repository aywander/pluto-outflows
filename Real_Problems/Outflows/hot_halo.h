//
// Created by Alexander Y. Wagner on 4/7/16.
//

#ifndef PLUTO_HOT_HALO_H
#define PLUTO_HOT_HALO_H

#include "init_tools.h"

void HaloOuterBoundary(const int side, const Data *d, int i, int j, int k, Grid *grid, int *touch);

void HotHaloPrimitives(double *halo, const double x1, const double x2, const double x3);

int InFlankRegion(const double x1, const double x2, const double x3);

#endif //PLUTO_HOT_HALO_H
