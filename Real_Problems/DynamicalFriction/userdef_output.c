#include "pluto.h"
#include "macros_usr.h"
#include "init_tools.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
    int i, j, k, nv;
    double ***te, ***spd, ***lmd;
    double ***v1, ***v2, ***v3;
    double ***prs, ***rho, ***vx1, ***vx2, ***vx3;
    double mu, vs1, vs2, vs3;
    double v[NVAR], rhs[NVAR];
    double *x1, *x2, *x3;

    /* New variables - names must exist under uservar */
    te = GetUserVar("te");

    /* Change to v instead of u = lorentz v */
    EXPAND(v1 = GetUserVar("v1");,
           v2 = GetUserVar("v2");,
           v3 = GetUserVar("v3"););

    /* State variables */
    rho = d->Vc[RHO];
    prs = d->Vc[PRS];
    EXPAND(vx1 = d->Vc[VX1];,
           vx2 = d->Vc[VX2];,
           vx3 = d->Vc[VX3];);

    /* These are the geometrical central points */
    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    DOM_LOOP(k, j, i) {

                /* Temperature */
                for (nv = 0; nv < NVAR; nv++) v[nv] = d->Vc[nv][k][j][i];
                te[k][j][i] = prs[k][j][i] / rho[k][j][i] * MU_NORM * KELVIN;

                vs1 = vs2 = vs3 = 0;
                EXPAND(vs1 = vx1[k][j][i];,
                       vs2 = vx2[k][j][i];,
                       vs3 = vx3[k][j][i];);

#if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
                /* This is useful in polar and spherical geometries, where the vectors are rotated. */

                EXPAND(double vc1 = VCART1(x1[i], x2[j], x3[k], vs1, vs2, vs3);,
                       double vc2 = VCART2(x1[i], x2[j], x3[k], vs1, vs2, vs3);,
                       double vc3 = VCART3(x1[i], x2[j], x3[k], vs1, vs2, vs3););

                EXPAND(vs1 = vc1;,
                       vs2 = vc2;,
                       vs3 = vc3;);

#endif

                /* vs[123] is now always cartesian */

                /* Normalize to km / s */
                EXPAND(v1[k][j][i] = vs1 * vn.v_norm / 1.e5;,
                       v2[k][j][i] = vs2 * vn.v_norm / 1.e5;,
                       v3[k][j][i] = vs3 * vn.v_norm / 1.e5;);

            }

}

/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





