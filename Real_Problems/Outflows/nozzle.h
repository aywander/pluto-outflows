//
// Created by Alexander Y. Wagner on 2019-01-02.
//

#ifndef PLUTO_NOZZLE_H
#define PLUTO_NOZZLE_H

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
    double vol_counted;       // Is the nozzle halved (for cylindrical symmetry?)
    double cone_height;       // Height of cone. Always positive.
    double cone_apex;         // Location of apex relative to (0,0,0)
                              // for a cone whose apex is on the flow axis.
                              // Can be used for a cone after RotateGrid2Nozzle transform.
                              // Can be positive or negative.
    int is_fan;               // Is the nozzle a fan (conical) or bullet-shaped (parallel)
    int is_two_sided;         // Is the nozzle two-sided?
} Nozzle;

extern Nozzle nz;

int InNozzleCap(double x1, double x2, double x3);

int InNozzleRegion(double x1, double x2, double x3);

void SetNozzleGeometry(Nozzle *noz);

void NozzleFill(Data *d, const Grid *grid);

double Profile(double x1, double x2, double x3);

void NozzleVolume(Data *d, Grid *grid);

void InitDomainNozzle(Data *d, Grid *grid);

void ClearNozzleSurrounding(double *cell, const double *halo, double x1, double x2, double x3);

#endif //PLUTO_NOZZLE_H
