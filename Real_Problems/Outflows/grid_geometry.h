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

int SphereSurfaceIntersectsCellByRadius(const double x1, const double x2, const double x3, const double vol,
                                        const double r);

int SphereSurfaceIntersectsCellByCorners(const double x1, const double x2, const double x3,
                                         const double dx1, const double dx2, const double dx3,
                                         const double r);

int SphereSurfaceIntersectsDomain(const Grid *grid, double r);

int SphereIntersectsDomain(const Grid *grid, const double r);

int BoxIntersectsDomain(const Grid *grid,
                        const double x1i, const double x1f,
                        const double x2i, const double x2f,
                        const double x3i, const double x3f);

int PointInDomain(const Grid *grid, const double x1, const double x2, const double x3);

double FindDxMax(const Grid *grid);

double ElevateCellVolume(Grid *grid, int i, int j, int k);

double ElevateVolume(double vol);

#endif //PLUTO_GRID_GEOMETRY_H
