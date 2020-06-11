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

double InterpolationWrapper2D(const double arg_i_arr[], const double arg_j_arr[], const double **val_arr,
                              const int ni, const int nj, const double arg_i, const double arg_j);

void InterpolateGrid(const Data *data, const Grid *grid, int *vars, double x1, double x2, double x3, double *v);

double BilinearInterpolate(const double val1, const double val2_i, const double val2_j, const double val2_ij,
                           const double frac_i, const double frac_j);

int UniformSamplingSphericalSurface(const int npoints, const double radius, double *x1, double *x2, double *x3);

void RandomSamplingSphericalSurface(const int npoints, const double radius, double *x1, double *x2, double *x3);

#endif
