//
// Created by Alexander Y. Wagner on 4/6/16.
//

#ifndef PLUTO_CLOUDS_H
#define PLUTO_CLOUDS_H


typedef struct {
    int nrad;
    double radii[3];
    double mdot_out[3];
    double pdot_out[3];
    double edot_out[3];
    double vrho_out[3];
    double mdot_out_dom;
    double pdot_out_dom;
    double edot_out_dom;
    double vrho_out_dom;
    double mass_tr2_dom;
    double mass_rtc_dom;
    double sfr;
    double porosity;
} CloudAnalytics;

extern CloudAnalytics ca;

int CloudCubePixel(int *el, const double x1, const double x2, const double x3);

void NormalizeFractalData(double *cloud, const double x1, const double x2, const double x3);

void CloudDensity(double *cloud, const double x1, const double x2, const double x3);

int CloudExtract(double *cloud, const double *halo, const int *pixel,
                 const double x1, const double x2, const double x3);

double CloudExtractEllipsoid(double fdratio, const double x1, const double x2, const double x3);

double CloudExtractCentralBuffer(double fdratio, const double x1, const double x2, const double x3);

void CloudVelocity(double *cloud, double *halo, const double x1, const double x2, const double x3);

int CloudPrimitives(double *cloud, const double x1, const double x2, const double x3);

int WarmTcrit(double *const warm);

void CloudAnalysis(Data *d, Grid *grid);

void WarmOutflowRates(Data *d, Grid *grid);

void WarmPhaseMass(Data *d, Grid *grid);

void WarmPhaseMassByRhoTcut(Data *d, Grid *grid);

void WarmPhaseStarFormationRate(Data *d, Grid *grid);

void WarmPhaseConditions(double *rho_c, double *te_c, double *tr2_c);

void WarmPhasePorosity(Data *d, Grid *grid);

void CloudOutput();

void InputDataClouds(const Data *d, const Grid *grid);

#endif //PLUTO_CLOUDS_H
