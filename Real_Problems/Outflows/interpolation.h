#ifndef interpolation_h
#define interpolation_h

int hunter(const double *arr, const int narr, const double val);

int hunter2(const double *arr, const int narr, const double val);

double LinearInterpolate(double y1, double y2, double fc);

double CosineInterpolate(double y1, double y2, double fc);

double CubicInterpolate(double y0, double y1, 
                        double y2, double y3, 
                        double fc);

double CubicCatmullRomInterpolate(double y0, double y1, 
                                  double y2, double y3, 
                                  double fc);


double HermiteInterpolate(double y0, double y1,
                          double y2, double y3,
                          double fc, double tension, 
                          double bias);

int sgn(double val);

int locate(double *,double, int);

double InterpolationWrapper(const double arg_arr[], const double val_arr[], const int narr, const double arg);

void InterpolateGrid(const Data *data, const Grid *grid, int *vars, double x1, double x2, double x3, double *v);

int UniformSamplingSphericalSurface(const int npoints, const double radius, double *x1, double *x2, double *x3);

void RandomSamplingSphericalSurface(const int npoints, const double radius, double *x1, double *x2, double *x3);

#endif
