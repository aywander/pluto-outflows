//
// Created by Alexander Y. Wagner on 4/6/16.
//

#include "supernovae.h"
#include "pluto.h"
#include "pluto_usr.h"
#include "grid_geometry.h"
#include "outflow.h"


/* ************************************************ */
int RotateGrid2Nozzle(
        const double cx1, const double cx2, const double cx3,
        double *cx1p, double *cx2p, double *cx3p) {
/*
 * This function assumes 3D (cartesian) coor-
 * dinates. First rotates according to any precession
 * through phi + omg*t. Then translates a point by dbh.
 * Then rotates grid around y-axis by angle dir so that
 * it is parallel to the nozzle cone.
 * Then point is translated back by dbh again.
 *
 *
 * Thus, the transformation matrix is defined as
 *
 *  / cx1p \              / cx1 \
 * |  cx2p  | = iT R T P |  cx2  |
 * |  cx3p  |            |  cx3  |
 *  \  1   /              \  1  /
 *
 *  where
 *
 *  Inverse rotation about Z axis
 *
 *       / cos(pre)   sin(pre)  0  0 \
 *  P = |  -sin(pre)  cos(pre)  0  0  |
 *      |   0            0      1  0  |
 *       \  0            0      0  1 /
 *
 *       where pre = phi + omg*g_time
 *
 *  Inverse rotation about Y axis
 *
 *       /  cos(dir)   0  -sin(dir)  0 \
 *  R = |      0       1      0      0  |
 *      |   sin(dir)   0   cos(dir)  0  |
 *       \     0       0      0      1 /
 *
 *  Translation along z
 *
 *       / 1   0  0   0 \
 *  T = |  0   1  0   0  |
 *      |  0   0  1 -dbh |
 *       \ 0   0  0   1 /
 *
 *       / 1   0  0   0 \
 * iT = |  0   1  0   0  |
 *      |  0   0  1  dbh |
 *       \ 0   0  0   1 /
 *
 ************************************************** */

    /* Precession angle */
    double pre;
    pre = nz.phi + nz.omg * g_time;

    /* The transformations */
    // Sign for nz.dbh translations may be the wrong way around.

    SELECT(*cx1p = cx1;,
           *cx1p = cx1 * cos(nz.dir) - sin(nz.dir) * (cx2 - nz.dbh);,
           *cx1p = (cx1 * cos(pre) - cx2 * sin(pre)) * cos(nz.dir) - sin(nz.dir) * (cx3 - nz.dbh););

    SELECT(*cx2p = cx2;,
           *cx2p = cx1 * sin(nz.dir) + nz.dbh + cos(nz.dir) * (cx2 - nz.dbh);,
           *cx2p = cx2 * cos(pre) + cx1 * sin(pre););

    SELECT(*cx3p = cx3;,
           *cx3p = cx3;,
           *cx3p = (cx1 * cos(pre) - cx2 * sin(pre)) * sin(nz.dir) + nz.dbh + cos(nz.dir) * (cx3 - nz.dbh););

    return 0;
}

/* ************************************************ */
int RotateNozzle2Grid(
        const double cx1, const double cx2, const double cx3,
        double *cx1p, double *cx2p, double *cx3p) {
/*
 * This Function does the inverse transformation of
 * RotateGrid2Nozzle
 *
 ************************************************** */

    /* Precession angle */
    double pre;
    pre = nz.phi + nz.omg * g_time;

    SELECT(*cx1p = cx1;,
           *cx1p = cx1 * cos(nz.dir) + sin(nz.dir) * (cx2 + nz.dbh);,
           *cx1p = cx2 * sin(pre) + cos(pre) * (cx1 * cos(nz.dir) + sin(nz.dir) * (cx3 + nz.dbh)););

    SELECT(*cx2p = cx2;,
           *cx2p = -(cx1 * sin(nz.dir)) - nz.dbh + cos(nz.dir) * (cx2 + nz.dbh);,
           *cx2p = cx2 * cos(pre) - sin(pre) * (cx1 * cos(nz.dir) + sin(nz.dir) * (cx3 + nz.dbh)););

    SELECT(*cx3p = cx3;,
           *cx3p = cx3;,
           *cx3p = -(cx1 * sin(nz.dir)) - nz.dbh + cos(nz.dir) * (cx3 + nz.dbh););

    return 0;
}


/* ************************************************ */
int SphereSurfaceIntersectsNCells(const double dx1, const double dx2, const double dx3, const double r) {
/* Returns the approximate expected number of cells intersected by the
 * SphereIntersectsCell routine in SIC_METHOD modes.
 *
 * Sphere is assumed to be at (0,0,0) and has a radius
 * of rad.
 ************************************************** */

    D_EXPAND(int fx1 = (int) (r / dx1);,
             int fx2 = (int) (r / dx2);,
             int fx3 = (int) (r / dx3););

    /* These expressions derive from the projection of the point, rim, surface, along the
     * first, first and second, first and second and thrid dimensions, respectively.
     * The expression for the 3D case is not exact, and slightly overestimates the number of cells,
     * but is correct within ~1%. */
    int ncells = D_SELECT(2;,
                 4 * (fx1 + fx2 + 1);,
                 (int) (2 * CONST_PI * ((fx1 + 1) * (fx2 + 1) +
                                        (fx1 + 1) * (fx3 + 1) +
                                        (fx2 + 1) * (fx3 + 1))););

#if SIC_METHOD == SIC_CORNERS

    return ncells;

#elif SIC_METHOD == SIC_RADIUS

    /* This is the ratio of surface area probed by the prediction of SphereSurfaceIntersectsNCells
     * which applies for SIC_METHOD == SIC_CORNER, giving a spherical surface of face-connected cells
     * to the true surface area of the sphere at which accretion is measured (4 pi r^2).
     * The SIC_RADIUS method has been verified to return the number of cells with this geometric correction. */

    static double geom_correction;
    geom_correction = D_SELECT(1, CONST_PI / 4, 2. / 3);

    return (int) (ncells * geom_correction);

#elif SIC_METHOD == SIC_HYBRID

    static double geom_correction;
    geom_correction = D_SELECT(1, CONST_PI / 4, 2. / 3);

    /* Return harmonic mean */
    return 2. / (1. / ncells + 1. / (ncells * geom_correction));

#endif

}


/* ************************************************ */
int SphereSurfaceIntersectsCellByRadius(const double x1, const double x2, const double x3,
                                        const double dV1, const double dV2, const double dV3,
                                        const double r) {
/*
 * Returns 1 if sphere surface intersects a cartesian cell.
 * An effective radius range is used.
 * Currently only for Cartesian and cylindrical case,
 * but other geometries are not difficult
 *
 * Sphere is assumed to be at (0, 0, 0) and has a radius
 * of rad. Also assumes that (0, 0, 0) is a corner point.
 ************************************************** */



    /* Get volume of cell */
    double vol = D_EXPAND(dV1, * dV2, * dV3);

    /* Equivalent volume spherical radius of cell */
    double rad = D_SELECT(vol, sqrt(vol / CONST_PI), pow(vol * 3 / (4 * CONST_PI), 1/3.));


    /* Sphere center assumed to be at (0, 0, 0).
     * To generalize, convert current coords to cartesians,
     * and bring sphere to center */

    /* Spherical distance to cell center */
    double sph1 = SPH1(x1, x2, x3);

    /* Cell is in sphere rim if the radial distance to the cell is r +- rad
     * The constants that reduce the radius to obtain roughly the right number of cells
     * were found experimentally. In this method, the weighting of the cells nearer the
     * coordinate axes is inevitably bigger. NOTE: Not sure if constants are necessary */
    return fabs(r - sph1) < rad;
             //* D_SELECT(1., 0.91, 0.79);

}


/* ************************************************ */
int SphereSurfaceIntersectsCellByCorners(const double x1, const double x2, const double x3,
                                         const double dx1, const double dx2, const double dx3,
                                         const double r) {
/*
 * Returns 1 if sphere intersects a cartesian cell local domain.
 * The entire cartesian cell is considered.
 * Currently only for Cartesian and cylindrical case,
 * but other geometries are not difficult
 *
 * Sphere is assumed to be at (0, 0, 0) and has a radius
 * of rad. Also assumes that (0, 0, 0) is a corner point.
 ************************************************** */


    /* Cell corner test */

    /* Cell extents */
    D_EXPAND(double x1i = x1 - dx1 / 2.; double x1f = x1 + dx1 / 2.;,
             double x2i = x2 - dx2 / 2.; double x2f = x2 + dx2 / 2.;,
             double x3i = x3 - dx3 / 2.; double x3f = x3 + dx3 / 2.;);

    /* Spherical radii of all corners */
    double ciii = SPH1(x1i, x2i, x3i);
    double ciif = SPH1(x1i, x2i, x3f);
    double cifi = SPH1(x1i, x2f, x3i);
    double ciff = SPH1(x1i, x2f, x3f);
    double cfii = SPH1(x1f, x2i, x3i);
    double cfif = SPH1(x1f, x2i, x3f);
    double cffi = SPH1(x1f, x2f, x3i);
    double cfff = SPH1(x1f, x2f, x3f);


    /* Minimum and maximum corner radii */
    double rmax = MAX(MAX(MAX(MAX(MAX(MAX(MAX(ciii, ciif), cifi), ciff), cfii), cfif), cffi), cfff);
    double rmin = MIN(MIN(MIN(MIN(MIN(MIN(MIN(ciii, ciif), cifi), ciff), cfii), cfif), cffi), cfff);

    /* Assumes sphere is at (0, 0, 0) */
    return (r > rmin) && (r < rmax);

}

/* ************************************************ */
int SphereSurfaceIntersectsDomain(struct GRID *grid, double r) {
/*
 * Returns 1 if sphere surface intersects local domain.
 * Only for INTERNAL_BOUNDARY YES
 *
 * Currently only for Cartesian and cylindrical case,
 * but other geometries are not difficult
 *
 * grid is the pointer to the grid struct of the local domain.
 *
 * Sphere is assumed to be at (0,0,0) and has a radius
 * of rad.
 ************************************************** */

    if (SphereIntersectsDomain(grid, r)) {

        /* The local domain limits */
        double x1i = grid[IDIR].xi;
        double x1f = grid[IDIR].xf;
        double x2i = grid[JDIR].xi;
        double x2f = grid[JDIR].xf;
        double x3i = grid[KDIR].xi;
        double x3f = grid[KDIR].xf;

        /* Check if at least one of the corners is outside of sphere */
        if (SPH1(x1i, x2i, x3i) > r) return 1;
        if (SPH1(x1i, x2i, x3f) > r) return 1;
        if (SPH1(x1i, x2f, x3i) > r) return 1;
        if (SPH1(x1i, x2f, x3f) > r) return 1;
        if (SPH1(x1f, x2i, x3i) > r) return 1;
        if (SPH1(x1f, x2i, x3f) > r) return 1;
        if (SPH1(x1f, x2f, x3i) > r) return 1;
        if (SPH1(x1f, x2f, x3f) > r) return 1;

        return 0;

    }
    else {
        return 0;
    }
}

/* ************************************************ */
int SphereIntersectsDomain(struct GRID *grid, const double r) {
/*
 * Returns 1 if solid sphere intersects local domain.
 * This considers the entire sphere, surface and body.
 * Only for INTERNAL_BOUNDARY YES
 *
 * Currently only for Cartesian and cylindrical case,
 * but other geometries are not difficult
 *
 * grid is an array of grid structures of the local domain
 *
 * Sphere is assumed to be at (0,0,0) and has a radius
 * of rad.
 ************************************************** */


    /* The local domain limits */
    double x1i = grid[IDIR].xi;
    double x1f = grid[IDIR].xf;
    double x2i = grid[JDIR].xi;
    double x2f = grid[JDIR].xf;
    double x3i = grid[KDIR].xi;
    double x3f = grid[KDIR].xf;

    /* The local domain center */
    double x1c = (x1f + x1i) / 2.;
    double x2c = (x2f + x2i) / 2.;
    double x3c = (x3f + x3i) / 2.;

    /* Location of center of sphere - currently this is hardcoded here */
    D_EXPAND(double s1 = 0;,
             double s2 = 0;,
             double s3 = 0;);


#if SID_METHOD == SID_POINTS

    /* Check if center of domain is in sphere */
    if (SPH1(x1c, x2c, x3c) < r) return 1;

    /* Check if center of sphere is in domain */
    if (D_EXPAND(x1i < s1 && s1 < x1f, &&
                 x2i < s2 && s2 < x2f, &&
                 x3i < s3 && s3 < x3f)) return 1;

    /* Check if at least one of the corners is inside sphere */
    if (SPH1(x1i, x2i, x3i) < r) return 1;
    if (SPH1(x1i, x2i, x3f) < r) return 1;
    if (SPH1(x1i, x2f, x3i) < r) return 1;
    if (SPH1(x1i, x2f, x3f) < r) return 1;
    if (SPH1(x1f, x2i, x3i) < r) return 1;
    if (SPH1(x1f, x2i, x3f) < r) return 1;
    if (SPH1(x1f, x2f, x3i) < r) return 1;
    if (SPH1(x1f, x2f, x3f) < r) return 1;

#if DIMENSIONS > 1
    /* Check if at least one of the edge centers is inside sphere */
    if (SPH1(x1i, x2i, x3c) < r) return 1;
    if (SPH1(x1i, x2f, x3c) < r) return 1;
    if (SPH1(x1f, x2i, x3c) < r) return 1;
    if (SPH1(x1f, x2f, x3c) < r) return 1;
    if (SPH1(x1i, x2c, x3i) < r) return 1;
    if (SPH1(x1i, x2c, x3f) < r) return 1;
    if (SPH1(x1f, x2c, x3i) < r) return 1;
    if (SPH1(x1f, x2c, x3f) < r) return 1;
    if (SPH1(x1c, x2i, x3i) < r) return 1;
    if (SPH1(x1c, x2i, x3f) < r) return 1;
    if (SPH1(x1c, x2f, x3i) < r) return 1;
    if (SPH1(x1c, x2f, x3f) < r) return 1;
#endif

#if DIMENSIONS > 2
    /* Check if at least one of the face centers is inside sphere */
    if (SPH1(x1i, x2c, x3c) < r) return 1;
    if (SPH1(x1f, x2c, x3c) < r) return 1;
    if (SPH1(x1c, x2i, x3c) < r) return 1;
    if (SPH1(x1c, x2f, x3c) < r) return 1;
    if (SPH1(x1c, x2c, x3i) < r) return 1;
    if (SPH1(x1c, x2c, x3f) < r) return 1;
#endif

    /* Finally, if none of the above triggered, do a domain surface loop */

    /* These are the geometrical central points */
    int i, j, k;
    double *x1 = grid[IDIR].x;
    double *x2 = grid[JDIR].x;
    double *x3 = grid[JDIR].x;

    DOM_SURF_LOOP(k, j, i, if (SPH1(x1[i], x2[j], x3[k]) < r) {return 1; } );

    return 0;


#elif SID_METHOD == SID_REGIONS

    /* Check if center of domain is in sphere */
    if (SPH1(x1c, x2c, x3c) <= r) return 1;

    /* Check if center of sphere is in domain */
    if (D_EXPAND(x1i <= s1 && s1 <= x1f, &&
                 x2i <= s2 && s2 <= x2f, &&
                 x3i <= s3 && s3 <= x3f)) return 1;

    /* Roll sphere around faces, edges, and corners of domain,
     * to define limiting distance of center of sphere from
     * surface of domain. */

    /* Define regions and distances */

    /* Face regions */
    D_EXPAND(int rf1i = s1 < x1i, && s2 > x2i && s2 < x2f, && s3 > x3i && s3 < x3f);
    D_EXPAND(int rf1f = s1 > x1f, && s2 > x2i && s2 < x2f, && s3 > x3i && s3 < x3f);
    D_EXPAND(,int rf2i = s2 < x2i, && s1 > x1i && s1 < x1f && s3 > x3i && s3 < x3f);
    D_EXPAND(,int rf2f = s2 > x2f, && s1 > x1i && s1 < x1f && s3 > x3i && s3 < x3f);
    D_EXPAND(,,int rf3i = s3 < x3i && s1 > x1i && s1 < x1f && s2 > x2i && s2 < x2f);
    D_EXPAND(,,int rf3f = s3 > x3f && s1 > x1i && s1 < x1f && s2 > x2i && s2 < x2f);

    /* 1-D distances to domain walls from (s1, s2, s3) */
    D_EXPAND(double f1i = fabs(s1 - x1i); double f1f = fabs(s1 - x1f);,
             double f2i = fabs(s2 - x2i); double f2f = fabs(s2 - x2f);,
             double f3i = fabs(s3 - x3i); double f3f = fabs(s3 - x3f););

    int cond_faces = D_EXPAND((rf1i && f1i < r) || (rf1f && f1f < r), ||
                              (rf2i && f2i < r) || (rf2f && f2f < r), ||
                              (rf3i && f3i < r) || (rf3f && f3f < r));

#if DIMENSIONS > 1

    /* Edge regions */
    D_EXPAND(,int re12ii = s1 < x1i && s2 < x2i, && s3 > x3i && s3 < x3f);
    D_EXPAND(,int re12if = s1 < x1i && s2 > x2f, && s3 > x3i && s3 < x3f);
    D_EXPAND(,int re12fi = s1 > x1f && s2 < x2i, && s3 > x3i && s3 < x3f);
    D_EXPAND(,int re12ff = s1 > x1f && s2 > x2f, && s3 > x3i && s3 < x3f);
    D_EXPAND(,,int re13ii = s1 < x1i && s3 < x3i && s2 > x2i && s2 < x2f);
    D_EXPAND(,,int re13if = s1 < x1i && s3 > x3f && s2 > x2i && s2 < x2f);
    D_EXPAND(,,int re13fi = s1 > x1f && s3 < x3i && s2 > x2i && s2 < x2f);
    D_EXPAND(,,int re13ff = s1 > x1f && s3 > x3f && s2 > x2i && s2 < x2f);
    D_EXPAND(,,int re23ii = s2 < x2i && s3 < x3i && s1 > x1i && s1 < x1f);
    D_EXPAND(,,int re23if = s2 < x2i && s3 > x3f && s1 > x1i && s1 < x1f);
    D_EXPAND(,,int re23fi = s2 > x2f && s3 < x3i && s1 > x1i && s1 < x1f);
    D_EXPAND(,,int re23ff = s2 > x2f && s3 > x3f && s1 > x1i && s1 < x1f);

    /* 2-D distances squared to edges from (s1, s2, s3) */
    D_EXPAND(,double e12ii = f1i * f1i + f2i * f2i;
              double e12if = f1i * f1i + f2f * f2f;
              double e12fi = f1f * f1f + f2i * f2i;
              double e12ff = f1f * f1f + f2f * f2f;,
              double e13ii = f1i * f1i + f3i * f3i;
              double e13if = f1i * f1i + f3f * f3f;
              double e13fi = f1f * f1f + f3i * f3i;
              double e13ff = f1f * f1f + f3f * f3f;
              double e23ii = f2i * f2i + f3i * f3i;
              double e23if = f2i * f2i + f3f * f3f;
              double e23fi = f2f * f2f + f3i * f3i;
              double e23ff = f2f * f2f + f3f * f3f;);

    double r2 = r * r;

    int cond_edges = D_EXPAND(,(re12ii && e12ii < r2) ||
                               (re12if && e12if < r2) ||
                               (re12fi && e12fi < r2) ||
                               (re12ff && e12ff < r2), ||
                               (re13ii && e13ii < r2) ||
                               (re13if && e13if < r2) ||
                               (re13fi && e13fi < r2) ||
                               (re13ff && e13ff < r2) ||
                               (re23ii && e23ii < r2) ||
                               (re23if && e23if < r2) ||
                               (re23fi && e23fi < r2) ||
                               (re23ff && e23ff < r2) );

#endif

#if DIMENSIONS > 2

    /* Corner regions */
    D_EXPAND(,,int rciii = s1 < x1i && s2 < x2i && s3 < x3i);
    D_EXPAND(,,int rciif = s1 < x1i && s2 < x2i && s3 > x3f);
    D_EXPAND(,,int rcifi = s1 < x1i && s2 > x2f && s3 < x3i);
    D_EXPAND(,,int rciff = s1 < x1i && s2 > x2f && s3 > x3f);
    D_EXPAND(,,int rcfii = s1 > x1f && s2 < x2i && s3 < x3i);
    D_EXPAND(,,int rcfif = s1 > x1f && s2 < x2i && s3 > x3f);
    D_EXPAND(,,int rcffi = s1 > x1f && s2 > x2f && s3 < x3i);
    D_EXPAND(,,int rcfff = s1 > x1f && s2 > x2f && s3 > x3f);

    /* 3-D distances squared to corners from (s1, s2, s3) */
    D_EXPAND(,,double ciii = e12ii + f3i * f3i;
               double ciif = e12ii + f3f * f3f;
               double cifi = e12if + f3i * f3i;
               double ciff = e12if + f3f * f3f;
               double cfii = e12fi + f3i * f3i;
               double cfif = e12fi + f3f * f3f;
               double cffi = e12ff + f3i * f3i;
               double cfff = e12ff + f3f * f3f;);

    int cond_corners = D_EXPAND(,,(rciii && ciii < r2) ||
                                  (rciif && ciif < r2) ||
                                  (rciff && ciff < r2) ||
                                  (rcifi && cifi < r2) ||
                                  (rcfii && cfii < r2) ||
                                  (rcfif && cfif < r2) ||
                                  (rcffi && cffi < r2) ||
                                  (rcfff && cfff < r2))

#endif


    if (D_EXPAND(cond_faces, || cond_edges, || cond_corners)){

         return 1;
    }

    else return 0;

#endif
}

/* ************************************************ */
int BoxIntersectsDomain(struct GRID *grid,
                        const double b1i, const double b1f,
                        const double b2i, const double b2f,
                        const double b3i, const double b3f){
/*!
 * Check whether Box is in Domain
 * grid is an array of grid structures
 *
 * WARNING: untested
 *
 ************************************************** */


    /* The local domain limits */
    double x1i = grid[IDIR].xi;
    double x1f = grid[IDIR].xf;
    double x2i = grid[JDIR].xi;
    double x2f = grid[JDIR].xf;
    double x3i = grid[KDIR].xi;
    double x3f = grid[KDIR].xf;

    /* The local domain center */
    double x1c = (x1f + x1i) / 2.;
    double x2c = (x2f + x2i) / 2.;
    double x3c = (x3f + x3i) / 2.;

    /* The box center */
    double b1c = (b1f + b1i) / 2.;
    double b2c = (b2f + b2i) / 2.;
    double b3c = (b3f + b3i) / 2.;

    /* Check whether box is in domain */
    if (D_EXPAND(x1i < b1c && b1c < x1f, &&
                 x2i < b2c && b2c < x2f, &&
                 x3i < b3c && b3c < x3f)) return 1;

    /* Check whether domain is in box */
    if (D_EXPAND(b1i < x1c && x1c < b1f, &&
                 b2i < x2c && x2c < b2f, &&
                 b3i < x3c && x3c < b3f)) return 1;

    /* Compare extents of the two rectangular regions */
    if (SGN(b1i - x1f) == SGN(b1f - x1i) && SGN(b1i - x1f) != 0) return 0;
    if (SGN(b2i - x2f) == SGN(b2f - x2i) && SGN(b2i - x2f) != 0) return 0;
    if (SGN(b3i - x3f) == SGN(b3f - x3i) && SGN(b3i - x3f) != 0) return 0;

    return 1;
}

/* ****************************************************** */
double FindDxMax(const Grid *grid) {
/*!
 *   Find largest cell width in local grid
 *
 * ****************************************************** */

    double dxmax[DIMENSIONS], scrh;
    double dx_max = 0;

    for (int idim = 0; idim < DIMENSIONS; idim++) dxmax[idim] = 0;

    int k, j, i;
    DOM_LOOP(k, j, i) {

                        EXPAND(scrh = Length_1(i, j, k, grid); dxmax[IDIR] = MAX(dxmax[IDIR], scrh);,
                               scrh = Length_2(i, j, k, grid); dxmax[JDIR] = MAX(dxmax[JDIR], scrh);,
                               scrh = Length_3(i, j, k, grid); dxmax[KDIR] = MAX(dxmax[KDIR], scrh);
                        );

                    }

    for (int idim = 1; idim < DIMENSIONS; idim++) dx_max = MAX(dx_max, dxmax[idim]);

    return dx_max;
}