//
// Created by Alexander Y. Wagner on 4/6/16.
//

#ifndef PLUTO_GRID_GEOMETRY_H
#define PLUTO_GRID_GEOMETRY_H

int RotateGrid2Nozzle(const double cx1, const double cx2, const double cx3,
                      double *cx1p, double *cx2p, double *cx3p);

int RotateNozzle2Grid(const double cx1, const double cx2, const double cx3,
                      double *cx1p, double *cx2p, double *cx3p);

int SphereSurfaceIntersectsNCells(const double dx1, const double dx2, const double dx3, const double r);

int SphereSurfaceIntersectsCellByRadius(const double x1, const double x2, const double x3,
                                        const double dV1, const double dV2, const double dV3,
                                        const double r);

int SphereSurfaceIntersectsCellByCorners(const double x1, const double x2, const double x3,
                                         const double dx1, const double dx2, const double dx3,
                                         const double r);

int SphereSurfaceIntersectsDomain(struct GRID *grid, double r);

int SphereIntersectsDomain(struct GRID *grid, const double r);

int BoxIntersectsDomain(struct GRID *grid,
                        const double x1i, const double x1f,
                        const double x2i, const double x2f,
                        const double x3i, const double x3f);

int PointInDomain(struct GRID *grid, const double x1, const double x2, const double x3);

double FindDxMax(const Grid *grid);

#include "init_tools.h"

#endif //PLUTO_GRID_GEOMETRY_H
