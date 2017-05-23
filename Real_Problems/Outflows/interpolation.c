#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pluto.h"
#include "pluto_usr.h"
#include "interpolation.h"



int hunter2(const double *arr, const int narr, const double val) {

    /*
     * This function returns index, il, of an array, arr, such that
     * arr[il] < val < arr[il + 1]
     */

    int il = 0;
    int ih = narr;

    while (il != (ih - 1)){
        int im = (il + ih) / 2;
        if   (val <= arr[im])
            ih = im;
        else
            il = im;
    }
    return il;
}


int hunter(const double *arr, const int narr, const double val) {

    /*
     * This function returns index, il, of an array, arr, such that
     * arr[il] < val < arr[il + 1]
     */

    int il, ir, shift;
    double arrb, arre;
    double vl, vr;
    double ldelta, rdelta, delta, dprod;


    /* Beginning and end values */
    arrb = arr[0];
    arre = arr[narr - 1];

    /* Bounds check */
    if (val < arrb) {
        print1("Error: interpolation.c: hunter: interpolation out of (lower) bounds.\n");
        print1("val :  %e\n", val);
        print1("arrb:  %e\n", arrb);
        exit(1);
    }
    else if (val > arre) {
        print1("Error: interpolation.c: hunter: interpolation out of (upper) bounds.\n");
        print1("val  %e\n", val);
        print1("arre %e\n", arre);
        exit(1);
    }

    /* Initial linear guess, left and right indices */
    il = (int) (narr - 1) * (val - arrb) / (arre - arrb);
    ir = il + 1;

    /* Left and right values from arrays */
    vl = arr[il];
    vr = arr[ir];

    /* Partial and total differences in r */
    ldelta = val - vl;
    rdelta = vr - val;
    delta = vr - vl;

    /* If positive, then we're in the right interval */
    dprod = ldelta * rdelta;

    while (dprod < 0) {
        /* Calculate single cell shift direction */
        shift = sgn(ldelta);

        /* New indices */
        il += shift;
        ir = il + 1;

        /* Left and right values from arrays */
        vl = arr[il];
        vr = arr[ir];

        /* Partial and total differences in r */
        ldelta = val - vl;
        rdelta = vr - val;
        delta = vr - vl;

        /* If positive, then we're in the right interval */
        dprod = ldelta * rdelta;
    }

    return il;

}

double LinearInterpolate(double y1, double y2, double fc) {
    /* Linear interpolator */

    return (y1 * (1 - fc) + y2 * fc);
}


double CosineInterpolate(double y1, double y2, double fc) {
    /* Cosine interpolator */

    double fc2;
    double pi = 3.14159265358979323846;
    fc2 = (1. - cos(fc * pi)) / 2.;
    return (y1 * (1. - fc2) + y2 * fc2);
}


double CubicInterpolate(double y0, double y1, 
                        double y2, double y3, 
                        double fc) {
    /* Cubic interpolator */

    double a0, a1, a2, a3, fc2;

    fc2 = fc * fc;
    a0 = y3 - y2 - y0 + y1;
    a1 = y0 - y1 - a0;
    a2 = y2 - y0;
    a3 = y1;

    return (a0 * fc * fc2 + a1 * fc2 + a2 * fc + a3);
}


double CubicCatmullRomInterpolate(double y0, double y1, 
                                  double y2, double y3, 
                                  double fc) {
    /* Cubic interpolator Catmull-Rom Spline */

    double a0, a1, a2, a3, fc2;

    fc2 = fc * fc;
    a0 = -0.5 * y0 + 1.5 * y1 - 1.5 * y2 + 0.5 * y3;
    a1 = y0 - 2.5 * y1 + 2 * y2 - 0.5 * y3;
    a2 = -0.5 * y0 + 0.5 * y2;
    a3 = y1;

    return (a0 * fc * fc2 + a1 * fc2 + a2 * fc + a3);
}


double HermiteInterpolate(double y0, double y1,
                          double y2, double y3,
                          double fc, double tension, 
                          double bias) {
/* Hermite interpolator
 * Tension: 1 is high, 0 normal, -1 is low
 * Bias: 0 is even,
 *       positive is towards first segment,
 *       negative towards the other
 */


    double m0, m1, fc2, fc3;
    double a0, a1, a2, a3;

    fc2 = fc * fc;
    fc3 = fc2 * fc;
    m0 = (y1 - y0) * (1. + bias) * (1. - tension) / 2.;
    m0 += (y2 - y1) * (1. - bias) * (1. - tension) / 2.;
    m1 = (y2 - y1) * (1. + bias) * (1. - tension) / 2.;
    m1 += (y3 - y2) * (1. - bias) * (1. - tension) / 2.;
    a0 = 2. * fc3 - 3. * fc2 + 1.;
    a1 = fc3 - 2. * fc2 + fc;
    a2 = fc3 - fc2;
    a3 = -2. * fc3 + 3. * fc2;

    return (a0 * y1 + a1 * m0 + a2 * m1 + a3 * y2);
}

int sgn(double val) {
    /* The signum function */

    if (val > 0) return 1;

    else if (val < 0) return -1;

    else return 0;
}



double InterpolationWrapper(const double arg_arr[], const double val_arr[], const int narr, const double arg) {

/* Wrapper around hunter and CubicCatmullRomInterpolate */

    int il;
    double arg1, arg2, frac;
    double val0, val1, val2, val3;

    /* Find cell left index to interpolate at */
    il = hunter(arg_arr, narr, arg);

    /* Linear fractional location of arg in cell */
    arg1 = arg_arr[il];
    arg2 = arg_arr[il + 1];
    frac = (arg - arg1) / (arg2 - arg1);

    /* Interpolation */
    val0 = val_arr[il - 1];
    val1 = val_arr[il];
    val2 = val_arr[il + 1];
    val3 = val_arr[il + 2];
    return CubicCatmullRomInterpolate(val0, val1, val2, val3, frac);

}



int locate(double *xx,double x, int N) {

/*
Locate interval by bisection from an ordered table, based on NR sec. 3.4 (modified).
Returns j such that interval is: xx[j]<x<xx[j+1]. Returns j=0 or j=N-2 for out of bounds data.
*/

    int ju, jl, jm;
    int ascnd, j;

    jl = 0;
    ju = N;
    j = jl;


    ascnd = (xx[N - 1] >= xx[0]);
    while (ju - jl > 1) {
        jm = (ju + jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl = jm;
        else
            ju = jm;
    }
    if (x <= xx[0]) j = 0;
    else if (x >= xx[N - 1]) j = N - 2;
    else j = jl;
    return j;
}


/* ************************************************************* */
void x3u_3d_extrapol(double ***a, int kb, int i, int j, int k, Grid *grid)
/*
 *
 * Quadratic extrapolation in the upper boundary of X3
 *
 *************************************************************** */
{
    double *x1, *x2, *x3;
    double y0, y1, y2, z0, z1, z2;

    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    y0 = a[kb][j][i];
    y1 = a[kb - 1][j][i];
    y2 = a[kb - 2][j][i];

    z0 = x3[kb];
    z1 = x3[kb - 1];
    z2 = x3[kb - 2];

    a[k][j][i] = y0 + (y1 - y0) / (z1 - z0) * (x3[k] - z0) +
                 ((y2 - y1) / (z2 - z1) / (z2 - z0) -
                  (y1 - y0) / (z1 - z0) / (z2 - z0)) *
                 (x3[k] - z0) * (x3[k] - z1);

}



/* ************************************************************* */
void x2l_3d_extrapol (double ***a, int jb, int i, int j, int k, Grid *grid)
/*
 *
 * Quadratic extrapolation in the lower boundary of X2
 *
 *************************************************************** */
{
    double *x1, *x2, *x3;
    double y0, y1, y2, z0, z1, z2;

    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    y0 = a[k][jb][i];
    y1 = a[k][jb + 1][i];
    y2 = a[k][jb + 2][i];

    z0 = x2[jb];
    z1 = x2[jb + 1];
    z2 = x2[jb + 2];

    a[k][j][i] = y0 + (y1 - y0) / (z1 - z0) * (x3[k] - z0) +
                 ((y2 - y1) / (z2 - z1) / (z2 - z0) -
                  (y1 - y0) / (z1 - z0) / (z2 - z0)) *
                 (x3[k] - z0) * (x3[k] - z1);

}


void InterpolateGrid(const Data *data, const Grid *grid, int *vars, double x1, double x2, double x3, double *v){

    /* Get domain data and range */

    D_EXPAND(int gr_nx1 = grid[IDIR].np_tot;,
             int gr_nx2 = grid[JDIR].np_tot;,
             int gr_nx3 = grid[KDIR].np_tot;);
    D_EXPAND(double * gr_x1 = grid[IDIR].x;,
             double * gr_x2 = grid[JDIR].x;,
             double * gr_x3 = grid[KDIR].x;);

    /* Find left indices */
    int il = grid[IDIR].lbeg;
    int jl = grid[JDIR].lbeg;
    int kl = grid[KDIR].lbeg;
    D_EXPAND(il = hunter2(gr_x1, gr_nx1, x1);,
             jl = hunter2(gr_x2, gr_nx2, x2);,
             kl = hunter2(gr_x3, gr_nx3, x3););


    /* Define normalized coordinates between [0,1]: */

    double xx, yy, zz;
    xx = yy = zz = 0.0;

    D_EXPAND(xx = (x1 - gr_x1[il]) / (gr_x1[il + 1] - gr_x1[il]);,
             yy = (x2 - gr_x2[jl]) / (gr_x2[jl + 1] - gr_x2[jl]);,
             zz = (x3 - gr_x3[kl]) / (gr_x3[kl + 1] - gr_x3[kl]););


    /* Perform bi- or tri-linear interpolation. */

    int nv = 0;
    double ***V;

    while (vars[nv] != -1) {

        int inv = vars[nv];
        V = data->Vc[nv];

        D_EXPAND(
        v[inv] = V[kl][jl][il] * (1.0 - xx) * (1.0 - yy) * (1.0 - zz)
                  + V[kl][jl][il + 1] * xx * (1.0 - yy) * (1.0 - zz);,

        v[inv] += V[kl][jl + 1][il] * (1.0 - xx) * yy * (1.0 - zz)
                  + V[kl][jl + 1][il + 1] * xx * yy * (1.0 - zz);,

        v[inv] += V[kl + 1][jl][il] * (1.0 - xx) * (1.0 - yy) * zz
                  + V[kl + 1][jl][il + 1] * xx * (1.0 - yy) * zz
                  + V[kl + 1][jl + 1][il] * (1.0 - xx) * yy * zz
                  + V[kl + 1][jl + 1][il + 1] * xx * yy * zz;
        );

        nv ++;
    }
}

/* ************************************************ */
int UniformSamplingSphericalSurface(const int npoints, const double radius, double *x1, double *x2, double *x3) {
/*!
 * Uniform sampling version of RandomSamplingSphericalSurface
 *
 ************************************************** */

    int i, j, k;

    int lcount = 0;

    double r2, scrh;
    double cx1, cx2, cx3;
    double rad_sc = sqrt(radius);

    int npoints_dir;
    npoints_dir = (int) (sqrt(npoints / CONST_PI * 4.));

    double spacing = 2 * rad_sc / npoints_dir;

    for (int n1 = 0; n1 < npoints_dir; n1++) {
        for (int n2 = 0; n2 < npoints_dir; n2++) {

            double rvar1 = -rad_sc + (0.5 + n1) * spacing;
            double rvar2 = -rad_sc + (0.5 + n2) * spacing;

            // Rejection
            r2 = rvar1 * rvar1 + rvar2 * rvar2;
            if (r2 < rad_sc * rad_sc) {

                // Cartesian coordinates of points on sphere
                scrh = sqrt(rad_sc * rad_sc - rvar1 * rvar1 - rvar2 * rvar2);
                cx1 = 2. * rvar1 * scrh;
                cx2 = 2. * rvar2 * scrh;
                cx3 = rad_sc * rad_sc - 2. * r2;

                D_EXPAND(x1[lcount] = CART_1(cx1, cx2, cx3);,
                         x2[lcount] = CART_2(cx1, cx2, cx3);,
                         x3[lcount] = CART_3(cx1, cx2, cx3););

                lcount++;

                if (lcount == npoints) break;
            }
        }
        if (lcount == npoints) break;

    }

    return lcount;

}

/* ************************************************ */
void RandomSamplingSphericalSurface(const int npoints, const double radius, double *x1, double *x2, double *x3) {
/*!
 * Return random points on surface of sphere.
 * Use method by Marsaglia (1972), see:
 * http://mathworld.wolfram.com/SpherePointPicking.html
 *
 ************************************************** */

    int i, j, k;


    int lcount = 0;

    double r2, scrh;
    double cx1, cx2, cx3;
    double rad_sc = sqrt(radius);

    while (lcount < npoints) {
        double rvar1 = RandomNumber(-rad_sc, rad_sc);
        double rvar2 = RandomNumber(-rad_sc, rad_sc);

        // Rejection
        r2 = rvar1 * rvar1 + rvar2 * rvar2;
        if (r2 < rad_sc * rad_sc) {

            // Cartesian coordinates of points on sphere
            scrh = sqrt(rad_sc * rad_sc - rvar1 * rvar1 - rvar2 * rvar2);
            cx1 = 2. * rvar1 * scrh;
            cx2 = 2. * rvar2 * scrh;
            cx3 = rad_sc * rad_sc - 2. * r2;

            D_EXPAND(x1[lcount] = CART_1(cx1, cx2, cx3);,
                     x2[lcount] = CART_2(cx1, cx2, cx3);,
                     x3[lcount] = CART_3(cx1, cx2, cx3););

            lcount ++;
        }

    }

}