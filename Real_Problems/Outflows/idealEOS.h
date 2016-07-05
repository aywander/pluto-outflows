#ifndef idealEOS_h
#define idealEOS_h

double PresIdealEOS(const double dens, const double temp, const double mu);
double TempIdealEOS(const double dens, const double pres, const double mu);
double DensIdealEOS(const double pres, const double temp, const double mu);
double SoundSpeed2IdealEOS(const double dens, const double pres);

#endif
