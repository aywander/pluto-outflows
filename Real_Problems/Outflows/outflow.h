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
    double vol;               // Volume into which energy etc is dumped.
    double cone_height;       // Height of cone. Always positive.
    double cone_apex;         // Location of apex relative to (0,0,0)
                              // for a cone whose apex is on the flow axis.
                              // Can be used for a cone after RotateGrid2Nozzle transform.
                              // Can be positive or negative.
    int is_fan;               // Is the nozzle a fan (conical) or bullet-shaped (parallel)
    int is_two_sided;         // Is the nozzle two-sided?
} Nozzle;

extern Nozzle nz;

/* Structure for outflow state parameters. All are kept in code units. */
typedef struct {
    double pow;               // Power of outflow as given by input parameter
    double mdt;               // Mass related parameter of outflow as given by input parameter
    double spd;               // Speed related parameter of outflow as given by input parameter
    double prs;               // Pressure as calculated from pow, mdt, and speed
    double rho;               // Density as calculated from pow, mdt, and speed
    double eth;               // Enthalpy inj. rate as calculated from pow, mdt, and speed
    double kin;               // Initial kinetic energy inj. rate as calculated from pow, mdt, and speed
} OutflowState;

extern OutflowState os;

void SetNozzleGeometry(Nozzle *noz);

void OutflowPrimitives(double *out_primitives, const double x1, const double x2, const double x3);

void OutflowVelocity(double *out_primitives, double speed, const double x1, const double x2, const double x3);

void SetOutflowState(OutflowState *ofs);

void SetJetState(OutflowState *ofs);

void SetUfoState(OutflowState *ofs);

int InNozzleCap(const double x1, const double x2, const double x3);

int InNozzleRegion(const double x1, const double x2, const double x3);

double Profile(const double x1, const double x2, const double x3);

#include "init_tools.h"

#endif //PLUTO_OUTFLOW_H
