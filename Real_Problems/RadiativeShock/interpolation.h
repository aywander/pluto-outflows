#ifndef interpolation_h
#define interpolation_h

int hunter(const double arr[], const int narr, const double val);

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

#endif
