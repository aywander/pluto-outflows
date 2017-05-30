//
// Created by Alexander Y. Wagner on 4/6/16.
//

#ifndef PLUTO_CLOUDS_H
#define PLUTO_CLOUDS_H

int CloudCubePixel(int *el, const double x1, const double x2, const double x3);

void ReadFractalData();

void GetFractalData(double *cloud, const double x1, const double x2, const double x3);

void CloudDensity(double *cloud, const double x1, const double x2, const double x3);

int CloudExtract(double *cloud, const double *halo, const int *pixel,
                 const double x1, const double x2, const double x3);

double CloudExtractEllipsoid(double fdratio, const double x1, const double x2, const double x3);

double CloudExtractCentralBuffer(double fdratio, const double x1, const double x2, const double x3);

void CloudVelocity(double *cloud, double *halo, const double x1, const double x2, const double x3);

int CloudPrimitives(double *cloud, const double x1, const double x2, const double x3);

int WarmTcrit(double *const warm);

#include "init_tools.h"

#endif //PLUTO_CLOUDS_H
