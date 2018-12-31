#include "pluto.h"
#include "pluto_usr.h"
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
    double ***rho, ***vx1, ***vx2, ***vx3, ***prs;
    double ***fg1, ***fg2, ***fg3;
    double ***v1, ***v2, ***v3;
    double ***te, ***spd;
    double mu, sp1, sp2, sp3;
    double v[NVAR];
    double *x1, *x2, *x3;
    double xc1, xc2, xc3;
    double fg;

    /* New variables - names must exist under uservar */
    fg1 = GetUserVar("fg1");
    fg2 = GetUserVar("fg2");
    fg3 = GetUserVar("fg3");
    te = GetUserVar("te");
    spd = GetUserVar("spd");

#if DEBUG_USING_USERDEF_VARS

    double ***deltax1, ***deltax2, ***deltax3;
    deltax1 = GetUserVar("dx1");
    deltax2 = GetUserVar("dx2");
    deltax3 = GetUserVar("dx3");

    double *dx1, *dx2, *dx3;
    dx1 = grid->dx[IDIR];
    dx2 = grid->dx[JDIR];
    dx3 = grid->dx[KDIR];

    double ***xcart1, ***xcart2, ***xcart3;
    xcart1 = GetUserVar("xc1");
    xcart2 = GetUserVar("xc2");
    xcart3 = GetUserVar("xc3");

    double ***xsph1, ***xsph2, ***xsph3;
    xsph1 = GetUserVar("xs1");
    xsph2 = GetUserVar("xs2");
    xsph3 = GetUserVar("xs3");

    double ***vol;
    vol = GetUserVar("vol");

#endif

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
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    double cell_vol, rad;

    DOM_LOOP(k, j, i) {

                /* Geometrical factors */
                rad = SPH1(x1[i], x2[j], x3[k]);
                cell_vol = grid->dV[k][j][i]

#if DEBUG_USING_USERDEF_VARS

                deltax1[k][j][i] = dx1[i];
                deltax2[k][j][i] = dx2[j];
                deltax3[k][j][i] = dx3[k];

                xc1 = CART1(x1[i], x2[j], x3[k]);
                xc2 = CART2(x1[i], x2[j], x3[k]);
                xc3 = CART3(x1[i], x2[j], x3[k]);

                xcart1[k][j][i] = xc1;
                xcart2[k][j][i] = xc2;
                xcart3[k][j][i] = xc3;

                xsph1[k][j][i] = r;
                xsph2[k][j][i] = SPH2(x1[i], x2[j], x3[k]);
                xsph3[k][j][i] = SPH3(x1[i], x2[j], x3[k]);

                vol[k][j][i] = cell_vol;
#endif

                /* Gravitational force */

                fg = -rho[k][j][i] * cell_vol / (rad * rad);
                fg1[k][j][i] = VSPH2CART1(x1[i], x2[j], x3[k], fg, 0, 0);
                fg2[k][j][i] = VSPH2CART2(x1[i], x2[j], x3[k], fg, 0, 0);
                fg3[k][j][i] = VSPH2CART3(x1[i], x2[j], x3[k], fg, 0, 0);


                /* Temperature */
                for (nv = 0; nv < NVAR; nv++) v[nv] = d->Vc[nv][k][j][i];
                mu = MeanMolecularWeight(v);
                te[k][j][i] = prs[k][j][i] / rho[k][j][i] * mu;


                /* Speed */
                sp1 = sp2 = sp3 = 0;
                EXPAND(sp1 = vx1[k][j][i];,
                       sp2 = vx2[k][j][i];,
                       sp3 = vx3[k][j][i];);
                spd[k][j][i] = VMAG(x1[i], x2[j], x3[k], sp1, sp2, sp3);

#if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
                /* This is useful in polar and spherical geometries, where the vectors are rotated. */

                EXPAND(double csp1 = VCART1(x1[i], x2[j], x3[k], sp1, sp2, sp3);,
                       double csp2 = VCART2(x1[i], x2[j], x3[k], sp1, sp2, sp3);,
                       double csp3 = VCART3(x1[i], x2[j], x3[k], sp1, sp2, sp3););

                EXPAND(sp1 = csp1;,
                       sp2 = csp2;,
                       sp3 = csp3;);

#endif

                EXPAND(v1[k][j][i] = sp1;,
                       v2[k][j][i] = sp2;,
                       v3[k][j][i] = sp3);

            }

}

/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





