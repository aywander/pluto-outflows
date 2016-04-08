//
// Created by Alexander Y. Wagner on 4/7/16.
//

#ifndef PLUTO_OUTFLOW_H
#define PLUTO_OUTFLOW_H

/* Structure for nozzle physics. All are kept in code units. */
typedef struct {
    double ang;               // Half opening angle.
    double rad;               // Radius of cone.
    double dir;               // Angle direction of cone.
    double dbh;               // Location of rotation axis (BH) relative to (0,0,0).
                              // Can be positive or negative.
    double cbh;               // Height of cone from cap up to rotation axis.
                              // Mathematically always positive.
    double omg;               // Precession angular velocity.
    double phi;               // Precession starting angle.
    double sph;               // Radius of spherical region if INTERNAL_BOUNDARY = YES
    double orig;              // Height of beginning of domain as measured
                              // from the origin along flowaxis.
                              // = g_domBeg[FLOWAXIS(IDIR, JDIR, KDIR)];.
                              // Can be positive or negative.
    double area;              // Area through which flux is calculated.
    double cone_height;       // Height of cone. Always positive.
    double cone_apex;         // Location of apex relative to (0,0,0)
                              // for a cone whose apex is on the flow axis.
                              // Can be used for a cone after RotateGrid2Nozzle transform.
                              // Can be positive or negative.
    int isfan;                // Is the nozzle a fan (conical) or bullet-shaped (parallel)
    int fill[15];             // Useless, just to make the struct a power of 2.
} Nozzle;

extern Nozzle nz;

void SetNozzleGeometry();

void OutflowPrimitives(double *out_primitives, const double x1, const double x2, const double x3,
                       const double accr_rate);

void OutflowVelocity(double *out_primitives, double speed, const double x1, const double x2, const double x3);

void JetPrimitives(double *jet_primitives, const double x1, const double x2, const double x3, const double accr_rate);

void UfoPrimitives(double *ufo_primitives, const double x1, const double x2, const double x3, const double accr_rate);

int InNozzleCap(const double x1, const double x2, const double x3);

int InNozzleRegion(const double x1, const double x2, const double x3);

double Profile(const double x1, const double x2, const double x3);

#include "init_tools.h"

#endif //PLUTO_OUTFLOW_H
