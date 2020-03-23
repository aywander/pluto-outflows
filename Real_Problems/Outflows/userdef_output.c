#include "pluto.h"
#include "pluto_usr.h"
#include "idealEOS.h"
#include "init_tools.h"

/* *************************************************************** */
void ComputeUserVar(const Data *d, Grid *grid)
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
    double ***v1, ***v2, ***v3;
    double ***te, ***spd, ***lmd;
    double mu, sp1, sp2, sp3;
    double v[NVAR], rhs[NVAR];
    double *x1, *x2, *x3;
#if (PHYSICS == RHD || PHYSICS == RMHD) && RECONSTRUCT_4VEL == YES
    double vel, speed, lorentz;
#endif


#if COORDINATE_SYSTEM_DEBUG == TRUE
    double ***xs1, ***xs2, ***xs3;
    double ***xc1, ***xc2, ***xc3;
    D_EXPAND(xs1 = GetUserVar("xs1");,
             xs2 = GetUserVar("xs2");,
             xs3 = GetUserVar("xs3"););
    D_EXPAND(xc1 = GetUserVar("xc1");,
             xc2 = GetUserVar("xc2");,
             xc3 = GetUserVar("xc3"););
#endif



    /* New variables - names must exist under uservar */
    te = GetUserVar("te");
    spd = GetUserVar("spd");
#if COOLING != NONE
    lmd = GetUserVar("lmd");
#endif

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
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    DOM_LOOP(k, j, i) {

#if COORDINATE_SYSTEM_DEBUG == TRUE
                D_EXPAND(xc1[k][j][i] = CART1(x1[i], x2[j], x3[k]);,
                         xc2[k][j][i] = CART2(x1[i], x2[j], x3[k]);,
                         xc3[k][j][i] = CART3(x1[i], x2[j], x3[k]););
                D_EXPAND(xs1[k][j][i] = x1[i];,
                         xs2[k][j][i] = x2[j];,
                         xs3[k][j][i] = x3[k];);
#endif

                /* Temperature */
                NVAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];
                mu = MeanMolecularWeight(v);
                te[k][j][i] = TempIdealEOS(rho[k][j][i], prs[k][j][i], mu) * vn.temp_norm;

#if COOLING != NONE
                /* Cooling rate */
                v[RHOE] /= g_gamma - 1.;
                Radiat(v, rhs);
                double lmd_norm = vn.pres_norm / vn.t_norm;
                lmd[k][j][i] = rhs[RHOE] * lmd_norm;
#endif

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

#if (PHYSICS == RHD || PHYSICS == RMHD) && RECONSTRUCT_4VEL == YES
                /* speed at this point is gamma * vel.
                 * Solve for vel. Then get gamma.
                 * Note, c=1 if RECONSTRUCT_4VEL = YES. */

                speed = spd[k][j][i];
                if (speed > 0){
                    vel = speed / sqrt(1 + speed * speed);
                    lorentz = speed / vel;
                }
                else{
                    vel = 0.;
                    lorentz = 1.;
                }

                spd[k][j][i] = vel;
                EXPAND(v1[k][j][i] = sp1 / lorentz;,
                       v2[k][j][i] = sp2 / lorentz;,
                       v3[k][j][i] = sp3 / lorentz;);
#endif

                /* Normalize to km / s */
                EXPAND(v1[k][j][i] = sp1 * vn.v_norm / ini_cgs[PAR_WTRB];,
                       v2[k][j][i] = sp2 * vn.v_norm / ini_cgs[PAR_WTRB];,
                       v3[k][j][i] = sp3 * vn.v_norm / ini_cgs[PAR_WTRB];);

                spd[k][j][i] *= vn.v_norm / ini_cgs[PAR_WTRB];

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

    /* HDF5 output cannot be controlled yet. Everything is output.*/

    /* VTK output */
    EXPAND(SetOutputVar("v1", VTK_OUTPUT, YES);,
           SetOutputVar("v2", VTK_OUTPUT, YES);,
           SetOutputVar("v3", VTK_OUTPUT, YES););
    EXPAND(SetOutputVar("vx1", VTK_OUTPUT, NO);,
           SetOutputVar("vx2", VTK_OUTPUT, NO);,
           SetOutputVar("vx3", VTK_OUTPUT, NO););
    SetOutputVar("prs", VTK_OUTPUT, YES);
    SetOutputVar("tr1", VTK_OUTPUT, YES);
#if CLOUDS == YES
    SetOutputVar("tr2", VTK_OUTPUT, YES);
#endif
    SetOutputVar("te", VTK_OUTPUT, YES);
    SetOutputVar("spd", VTK_OUTPUT, YES);
    SetOutputVar("lmd", VTK_OUTPUT, YES);

    /* FLT output */
    EXPAND(SetOutputVar("v1", FLT_OUTPUT, YES);,
           SetOutputVar("v2", FLT_OUTPUT, YES);,
           SetOutputVar("v3", FLT_OUTPUT, YES););
    EXPAND(SetOutputVar("vx1", FLT_OUTPUT, NO);,
           SetOutputVar("vx2", FLT_OUTPUT, NO);,
           SetOutputVar("vx3", FLT_OUTPUT, NO););
    SetOutputVar("prs", FLT_OUTPUT, YES);
    SetOutputVar("tr1", FLT_OUTPUT, YES);
#if CLOUDS == YES
    SetOutputVar("tr2", FLT_OUTPUT, YES);
#endif
    SetOutputVar("te", FLT_OUTPUT, NO);
    SetOutputVar("spd", FLT_OUTPUT, NO);
    SetOutputVar("lmd", FLT_OUTPUT, NO);


    /* FLT H5 output */
    EXPAND(SetOutputVar("v1", FLT_H5_OUTPUT, YES);,
           SetOutputVar("v2", FLT_H5_OUTPUT, YES);,
           SetOutputVar("v3", FLT_H5_OUTPUT, YES););
    EXPAND(SetOutputVar("vx1", FLT_H5_OUTPUT, NO);,
           SetOutputVar("vx2", FLT_H5_OUTPUT, NO);,
           SetOutputVar("vx3", FLT_H5_OUTPUT, NO););
    SetOutputVar("prs", FLT_H5_OUTPUT, YES);
    SetOutputVar("tr1", FLT_H5_OUTPUT, YES);
#if CLOUDS == YES
    SetOutputVar("tr2", FLT_H5_OUTPUT, YES);
#endif
    SetOutputVar("te", FLT_H5_OUTPUT, YES);
    SetOutputVar("spd", FLT_H5_OUTPUT, YES);
    SetOutputVar("lmd", FLT_H5_OUTPUT, NO);


    /* PNG output */
    EXPAND(SetOutputVar("v1", PNG_OUTPUT, NO);,
           SetOutputVar("v2", PNG_OUTPUT, NO);,
           SetOutputVar("v3", PNG_OUTPUT, NO););
    EXPAND(SetOutputVar("vx1", PNG_OUTPUT, NO);,
           SetOutputVar("vx2", PNG_OUTPUT, NO);,
           SetOutputVar("vx3", PNG_OUTPUT, NO););
    SetOutputVar("prs", PNG_OUTPUT, NO);
    SetOutputVar("tr1", PNG_OUTPUT, NO);
#if CLOUDS == YES
    SetOutputVar("tr2", PNG_OUTPUT, NO);
#endif
    SetOutputVar("te", PNG_OUTPUT, YES);
    SetOutputVar("spd", PNG_OUTPUT, YES);
    SetOutputVar("lmd", PNG_OUTPUT, NO);



    /* density slice */
    image = GetImage("rho");
#if COMPONENTS > 2
    image->slice_plane = X13_PLANE;
    image->slice_coord = 0.0;
#endif
    image->max = 1.e3; image->min = 1.e-2;
    image->logscale = 1;

    /* pressure slice */
    image = GetImage("te");
#if COMPONENTS > 2
    image->slice_plane = X13_PLANE;
    image->slice_coord = 0.0;
#endif
    image->max = 1.e7; image->min = 100.;
    image->logscale = 1;
    image->colormap = "br";


    /* first tracer slice */
#if NTRACER > 0
    image = GetImage("spd");
#if COMPONENTS > 2
    image->slice_plane = X13_PLANE;
    image->slice_coord = 0.0;
#endif
    image->max = 1.e5; image->min = 10;
    image->logscale = 1;
    image->colormap = "bw";
#endif

//#if CLOUDS
//    /* second tracer slice */
//#if NTRACER > 1
//    image = GetImage("tr2");
//#if COMPONENTS > 2
//    image->slice_plane = X13_PLANE;
//    image->slice_coord = 0.0;
//#endif
//    image->max = 1.; image->min = 0.;
//    image->logscale = 0;
//#endif
//#endif

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}


