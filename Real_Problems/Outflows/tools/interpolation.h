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

void x3u_3d_extrapol(double ***a, int kb, int i, int j, int k, Grid *grid);

void x2l_3d_extrapol (double ***a, int jb, int i, int j, int k, Grid *grid);

#endif
