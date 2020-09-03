#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"


#if BODY_FORCE != NO
/* *************************************************************** */
void ComputeUserVarDphidr(const Data *d, Grid *grid)
/*
 * Calculate radial (cylindrical) acceleration
 *
 ***************************************************************** */
{
    int i, j, k, nv;
    double ***dphidr;
    double gvec[COMPONENTS];
    double v[NVAR];

    dphidr = GetUserVar("dphidr");

    /* These are the geometrical central points */
    double *x1, *x2, *x3;
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

    DOM_LOOP(k, j, i) {
                NVAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];

                BodyForceVector(v, gvec, x1[i], x2[j], x3[k]);  // Acceleration vector (pointing inward)
                dphidr[k][j][i] = -VPOL1(x1[i], x2[j], x3[k], gvec[IDIR], gvec[JDIR], gvec[KDIR]);

            }

}
#endif


/* *************************************************************** */
void ComputeUserVarTe(const Data *d, Grid *grid)
/*
 * Calculate temperature
 *
 ***************************************************************** */
{
    int i, j, k, nv;
    double ***te;
    double v[NVAR];
    double mu;

    te = GetUserVar("te");

    DOM_LOOP(k, j, i) {
                NVAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];
                mu = MeanMolecularWeight(v);
                te[k][j][i] = v[PRS] / v[RHO] * mu * vn.temp_norm;
            }

}


#if COOLING != NONE
/* *************************************************************** */
void ComputeUserVarLmd(const Data *d, Grid *grid)
/*
 * Calculate Cooling rate
 *
 ***************************************************************** */
{
    int i, j, k, nv;
    double ***lmd;
    double v[NVAR], rhs[NVAR];

    lmd = GetUserVar("te");

    DOM_LOOP(k, j, i) {

                NVAR_LOOP(nv) v[nv] = d->Vc[nv][k][j][i];

                v[RHOE] /= g_gamma - 1.;
                Radiat(v, rhs);
                double lmd_norm = vn.pres_norm / vn.t_norm;
                lmd[k][j][i] = rhs[RHOE] * lmd_norm;

            }

}
#endif


/* *************************************************************** */
void ComputeUserVarSpd(const Data *d, Grid *grid)
/*
 * Speed Calculate
 *  - always without Lorentz factor
 *  - normalized to km / s
 *
 ***************************************************************** */
{
    int i, j, k;
    double ***spd;
    double vx1, vx2, vx3;

    spd = GetUserVar("spd");


    DOM_LOOP(k, j, i) {

                vx1 = vx2 = vx3 = 0;
                EXPAND(vx1 = d->Vc[VX1][k][j][i];,
                       vx2 = d->Vc[VX2][k][j][i];,
                       vx3 = d->Vc[VX3][k][j][i];);
                spd[k][j][i] = VMAG(x1[i], x2[j], x3[k], vx1, vx2, vx3);


                /* Normalize to km / s */
                spd[k][j][i] *= vn.v_norm / ini_cgs[PAR_WTRB];

            }

}


/* *************************************************************** */
void ComputeUserVarVel(const Data *d, Grid *grid)
/*
 *  Velocity vectors,
 *  - always in Cartesian
 *  - always without Lorentz factor
 *  - normalized to km / s
 *
 ***************************************************************** */
{
    int i, j, k;
    double spd;
    double ***v1, ***v2, ***v3;
    double vx1, vx2, vx3;

    EXPAND(v1 = GetUserVar("v1");,
           v2 = GetUserVar("v2");,
           v3 = GetUserVar("v3"););

#if GEOMETRY == POLAR || GEOMETRY == SPHERICAL

    /* These are the geometrical central points */
    double *x1, *x2, *x3;
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

#endif

    DOM_LOOP(k, j, i) {


                vx1 = vx2 = vx3 = 0;
                EXPAND(vx1 = d->Vc[VX1][k][j][i];,
                       vx2 = d->Vc[VX2][k][j][i];,
                       vx3 = d->Vc[VX3][k][j][i];);
                spd = VMAG(x1[i], x2[j], x3[k], vx1, vx2, vx3);

#if GEOMETRY == POLAR || GEOMETRY == SPHERICAL

                /* This is useful in polar and spherical geometries, where the vectors are rotated. */

                double cvx1, cvx2, cvx3;
                EXPAND(cvx1 = VCART1(x1[i], x2[j], x3[k], vx1, vx2, vx3);,
                       cvx2 = VCART2(x1[i], x2[j], x3[k], vx1, vx2, vx3);,
                       cvx3 = VCART3(x1[i], x2[j], x3[k], vx1, vx2, vx3););

                EXPAND(vx1 = cvx1;,
                       vx2 = cvx2;,
                       vx3 = cvx3;);

#endif

                /* Normalize to km / s */
//                EXPAND(v1[k][j][i] = vx1 * vn.v_norm / ini_cgs[PAR_WTRB];,
//                       v2[k][j][i] = vx2 * vn.v_norm / ini_cgs[PAR_WTRB];,
//                       v3[k][j][i] = vx3 * vn.v_norm / ini_cgs[PAR_WTRB];);

            }

}


/* *************************************************************** */
void ComputeUserVarCoord(const Data *d, Grid *grid)
/*
 * Calculate Cooling rate
 *
 ***************************************************************** */
{
    int i, j, k;
    double *x1, *x2, *x3;

    double ***xs1, ***xs2, ***xs3;
    double ***xc1, ***xc2, ***xc3;
    D_EXPAND(xs1 = GetUserVar("xs1");,
             xs2 = GetUserVar("xs2");,
             xs3 = GetUserVar("xs3"););
    D_EXPAND(xc1 = GetUserVar("xc1");,
             xc2 = GetUserVar("xc2");,
             xc3 = GetUserVar("xc3"););

    DOM_LOOP(k, j, i) {

                D_EXPAND(xc1[k][j][i] = CART1(x1[i], x2[j], x3[k]);,
                         xc2[k][j][i] = CART2(x1[i], x2[j], x3[k]);,
                         xc3[k][j][i] = CART3(x1[i], x2[j], x3[k]););
                D_EXPAND(xs1[k][j][i] = x1[i];,
                         xs2[k][j][i] = x2[j];,
                         xs3[k][j][i] = x3[k];);
            }

}



/* *************************************************************** */
void ComputeUserVar(const Data *d, Grid *grid)
/*
 *
 *  purpose
 *
 *    define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{

    /* New variables - names must exist under uservar */

    int ip;
    Runtime *runtime;
    runtime = RuntimeGet();

    /* Only do if it exists under uservar */
    for (ip = 0; ip < runtime->user_var; ip++){
        if (strcmp(runtime->user_var_name[ip], "te") == 0) ComputeUserVarTe(d, grid);
#if COOLING != NO
        if (strcmp(runtime->user_var_name[ip], "lmd") == 0) ComputeUserVarLmd(d, grid);
#endif
        if (strcmp(runtime->user_var_name[ip], "spd") == 0) ComputeUserVarSpd(d, grid);
        if (strcmp(runtime->user_var_name[ip], "v1") == 0) ComputeUserVarVel(d, grid);
        if (strcmp(runtime->user_var_name[ip], "xs1") == 0) ComputeUserVarCoord(d, grid);
#if BODY_FORCE != NO
        if (strcmp(runtime->user_var_name[ip], "dphidr") == 0) ComputeUserVarDphidr(d, grid);
#endif
    }

}


/* ************************************************************* */
void SetAllOutputVar(char *var, int yn){
/*
 * Helper function to set all output for a variable to YES or NO.
 * Used in ChangeOutputVar ()
 *
/* ************************************************************* */

    static int output_types_yes[] = {VTK_OUTPUT, FLT_OUTPUT, FLT_H5_OUTPUT, PNG_OUTPUT};
    static int output_types_no[]  = {};
    int iot;

    for (iot = 0; iot < sizeof(output_types_yes) / sizeof(int); iot++) SetOutputVar(var, output_types_yes[iot], yn);
    for (iot = 0; iot < sizeof(output_types_no) / sizeof(int); iot++) SetOutputVar(var, output_types_no[iot], NO);

}

/* ************************************************************* */
void SetAllOutputVarVec(char *var, int yn){
/*
 * Helper function to set all output for a vector variable to YES or NO.
 * Used in ChangeOutputVar (). For vector variables, we disable PNG output.
 *
/* ************************************************************* */

    static int output_types_yes[] = {VTK_OUTPUT, FLT_OUTPUT, FLT_H5_OUTPUT};
    static int output_types_no[]  = {PNG_OUTPUT};
    int iot;

    for (iot = 0; iot < sizeof(output_types_yes) / sizeof(int); iot++) SetOutputVar(var, output_types_yes[iot], yn);
    for (iot = 0; iot < sizeof(output_types_no) / sizeof(int); iot++) SetOutputVar(var, output_types_no[iot], NO);

}

/* ************************************************************* */
void SetImageAttr(Image *image, char *var, double min, double max, int log, char* cmap)
/*
 *  Set default image, with some defaults
 *
 *************************************************************** */
{

    image = GetImage(var);
#if COMPONENTS > 2
    /* Automatically choose which slice to take in 3D
       Default is X13_PLANE, unless it is very thin */
    double width_z = g_domEnd[KDIR] - g_domBeg[KDIR];
    double width_y = g_domEnd[JDIR] - g_domBeg[JDIR];
    double width_x = g_domEnd[IDIR] - g_domBeg[IDIR];
    if      (MIN(width_z / width_x, width_x / width_z) > 0.95) image->slice_plane = X13_PLANE;
    else if (MIN(width_y / width_x, width_x / width_y) > 0.95) image->slice_plane = X12_PLANE;
    else if (MIN(width_z / width_y, width_y / width_z) > 0.95) image->slice_plane = X23_PLANE;
    image->slice_coord = 0.0;
#endif
    image->min = min;
    image->max = max;
    image->logscale = log;
    image->colormap = cmap;

}

/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *  Purpose: Control output, variable by variable
 * 
 *************************************************************** */
{

    int ip;
    Runtime *runtime;
    runtime = RuntimeGet();

    Image *image;

    /* Only do if it exists under uservar */
    for (ip = 0; ip < runtime->user_var; ip++) {

        if (strcmp(runtime->user_var_name[ip], "te") == 0) {

            SetAllOutputVar("te", YES);
            SetImageAttr(image, "te", 30, 3.e5, 1, "parula");

        }

        if (strcmp(runtime->user_var_name[ip], "lmd") == 0) {

            SetAllOutputVar("lmd", YES);
            SetImageAttr(image, "lmd", 1.e23, 1.e25, 1, "jet");

        }

        if (strcmp(runtime->user_var_name[ip], "spd") == 0) {

            SetAllOutputVar("spd", YES);
            SetImageAttr(image, "spd", 10, 1000, 1, "blue");

        }

        if (strcmp(runtime->user_var_name[ip], "v1") == 0) {

            SetAllOutputVarVec("v1", YES);
            SetAllOutputVarVec("v2", YES);
            SetAllOutputVarVec("v3", YES);

            SetAllOutputVarVec("vx1", NO);
            SetAllOutputVarVec("vx2", NO );
            SetAllOutputVarVec("vx3", NO );

        }

        if (strcmp(runtime->user_var_name[ip], "xs1") == 0) {

        }

        if (strcmp(runtime->user_var_name[ip], "dphidr") == 0) {

            SetAllOutputVar("dphidr", YES);
            SetImageAttr(image, "rho", 1.e0, 1.e3, 1, "jet");

        }
    }


    /* density slice */
    SetAllOutputVar("rho", YES);
    SetImageAttr(image, "rho", 1.e-4, 1.e3, 1, "jet");

    /* Cloud Tracer */
    // TODO: This doesn't work - find out why
//#if CLOUDS && NTRACER > 1
//    SetImageAttr(image, "tr2", 0, 1, 0, "bw");
//#endif

#ifdef PARTICLES
    SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
    SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
    SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}


