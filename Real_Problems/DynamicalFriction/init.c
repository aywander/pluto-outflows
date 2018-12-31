/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"
#include "friction.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{

    double mach, rho, prs, dratio, rcore;
    double vsound;
    double vcx1, vcx2, vcx3;
    double vx1, vx2, vx3;
    double cx1, cx2, cx3;

    static int once = 0;

    /* Some things that only need to be done once */
    if (!once) {

        /* Initialize base normalization struct */
        SetBaseNormalization();
        PrintBaseNormalizations();

        /* Set normalization factors for input parameters */
        SetIniNormalization();

        /* Done once now */
        once = 1;
    }


    /* Velocity field in Cartesian coordinates.
     * Due to our chosen normalization, vsound = beta = 1 / √(1 - Mach^2). */

    vcx1 = vcx2 = vcx3 = 0;

    mach = g_inputParam[PAR_MACH];

    /* Choose normalization */
//     vsound = 1. / sqrt(1. + mach * mach);
    vsound = 1.;

    D_SELECT(,
            vcx2 = -vsound * mach;,
            vcx3 = -vsound * mach;);


    /* Velocity field in current coordinates */

    D_EXPAND(cx1 = CART1(x1, x2, x3);,
             cx2 = CART2(x1, x2, x3);,
             cx3 = CART3(x1, x2, x3););

    vx1 = vx2 = vx3 = 0;
    D_EXPAND(vx1 = VCART_1(cx1, cx2, cx3, vcx1, vcx2, vcx3);,
             vx2 = VCART_2(cx1, cx2, cx3, vcx1, vcx2, vcx3);,
             vx3 = VCART_3(cx1, cx2, cx3, vcx1, vcx2, vcx3););


    rho = g_inputParam[PAR_DENS];
    prs = vsound * vsound * rho / g_gamma;

    /* Fill primitive array. */
    v[RHO] = rho;
    v[VX1] = vx1;
    v[VX2] = vx2;
    v[VX3] = vx3;
    v[PRS] = prs;

}


/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

    /* Calculate component-wise friction force */
    SetFriction(grid);

    FrictionForce(d, grid);
    FrictionForceOutput();

    FrictionForceRadius(d, grid);
    FrictionForceRadiusOutput();

//    FrictionForcePolar(d, grid);
//    FrictionForcePolarOutput();

//    FrictionForceSlab(d, grid);
//    FrictionForceSlabOutput();


    /* Output component-wise friction force to file */

}

#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
    int   i, j, k, nv;
    double  *x1, *x2, *x3;
    double v[NVAR];

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];

    if (side == 0) {    /* -- check solution inside domain -- */
        DOM_LOOP(k,j,i){};
    }

    if (side == X1_BEG){  /* -- X1_BEG boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X1FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X2FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X3FACE){
            BOX_LOOP(box,k,j,i){  }
        }
    }

    if (side == X1_END){  /* -- X1_END boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box,k,j,i){

#if GEOMETRY == SPHERICAL

                        /* Free-flow lower boundary */
                        if ((g_time > 1.) && (x2[j] > CONST_PI / 2.)) {
                            for (nv = 0; nv < NVAR; ++nv) {
                                d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
                            }
                        }

                        /* Fixed boundary upper half */
                        else {
                            Init(v, x1[i], x2[j], x3[k]);

                            for (nv = 0; nv < NVAR; ++nv) {
                                d->Vc[nv][k][j][i] = v[nv];
                            }
                        }

#else
                        for (nv = 0; nv < NVAR; ++nv) {
                            d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
                        }
#endif

                    }
        }else if (box->vpos == X1FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X2FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X3FACE){
            BOX_LOOP(box,k,j,i){  }
        }
    }

    if (side == X2_BEG){  /* -- X2_BEG boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X1FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X2FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X3FACE){
            BOX_LOOP(box,k,j,i){  }
        }
    }

    if (side == X2_END){  /* -- X2_END boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X1FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X2FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X3FACE){
            BOX_LOOP(box,k,j,i){  }
        }
    }

    if (side == X3_BEG){  /* -- X3_BEG boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X1FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X2FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X3FACE){
            BOX_LOOP(box,k,j,i){  }
        }
    }

    if (side == X3_END){  /* -- X3_END boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X1FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X2FACE){
            BOX_LOOP(box,k,j,i){  }
        }else if (box->vpos == X3FACE){
            BOX_LOOP(box,k,j,i){  }
        }
    }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{

    /* Black hole is at (0, 0, 0) */

    /* Units of G * Mbh = 1. */

    double r, rc, gr;
    r = SPH1(x1, x2, x3);
    rc = g_inputParam[PAR_RCORE];
    gr = - 1. / (r * r + rc * rc);

    double sx1, sx2, sx3;
    D_EXPAND(sx1 = SPH1(x1, x2, x3);,
             sx2 = SPH2(x1, x2, x3);,
             sx3 = SPH3(x1, x2, x3););

    double g1, g2, g3;
    g1 = g2 = g3 = 0;
    EXPAND(g1 = VSPH_1(sx1, sx2, sx3, gr, 0, 0);,
           g2 = VSPH_2(sx1, sx2, sx3, gr, 0, 0);,
           g3 = VSPH_3(sx1, sx2, sx3, gr, 0, 0););

    g[IDIR] = g1;
    g[JDIR] = g2;
    g[KDIR] = g3;

}

/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{

    /* Black hole is at (0, 0, 0) */

    /* Units of G * Mbh = 1. */

    double rc;
    double r, pot;

    rc = g_inputParam[PAR_RCORE];
    r = SPH1(x1, x2, x3);
    pot = - 1. / (r + rc);

    return pot;

}
#endif
