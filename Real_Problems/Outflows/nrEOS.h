#ifndef nrEOS_h
#define nrEOS_h

double PresNrEOS(const double dens, const double temp, const double mu);
double TempNrEOS(const double dens, const double pres, const double mu);
double DensNrEOS(const double pres, const double temp, const double mu);

#endif
