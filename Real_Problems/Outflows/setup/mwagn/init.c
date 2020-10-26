/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.


  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"
#include "read_grav_table.h"
#include "interpolation.h"
#include "accretion.h"
#include "clouds.h"
#include "grid_geometry.h"
#include "hot_halo.h"
#include "outflow.h"
#include "supernovae.h"
#include "nozzle.h"


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

#if PHYSICS == MHD || PHYSICS == RMHD

    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = 0.0;

    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = 0.0;

#endif

}


/* ********************************************************************* */
void InitDomain(Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{

    /* Initialize base normalization struct */
    SetBaseNormalization();
    PrintBaseNormalizations();

    /* Set normalization factors for input parameters */
    SetIniNormalization();

    /* Get static gravity table if necessary */
#if GRAV_POTENTIAL == GRAV_TABLE || GRAV_POTENTIAL == GRAV_2D_TABLE
    ReadGravTable();
#endif

    /* All about the accretion - this, i think, must be done first.
     * Set outflow geometry struct with parameters of cone */
#if ACCRETION == YES
        SetAccretionPhysics();
#endif

    /* All about the nozzle */
#if NOZZLE != NONE

    /* Set outflow geometry struct with parameters of cone */
    SetNozzleGeometry(&nz);

    /* Set struct of outflow parameters and state variables */
    SetOutflowState(&os);

    double dx;
    dx = FLOWAXIS((g_domEnd[IDIR] - g_domBeg[IDIR]) / NX1;,
                  (g_domEnd[JDIR] - g_domBeg[JDIR]) / NX2;,
                  (g_domEnd[KDIR] - g_domBeg[KDIR]) / NX3;);

    /* Print some data */
    double halo_primitives[NVAR], out_primitives[NVAR];
    OutflowPrimitives(out_primitives, ARG_FLOWAXIS(dx, 0));
    HotHaloPrimitives(halo_primitives, ARG_FLOWAXIS(dx, 0));
    PrintInitData01(out_primitives, halo_primitives);

#endif

#if SUPERNOVAE == YES

    /* Set Supernova physics */
        SetSupernovaePhysics();

#endif


#if COOLING
    g_minCoolingTemp = 3.e2;
    g_maxCoolingRate = 0.4;
#endif



    /* Initialize a hot phase ambient medium by default
     * - this is the bare minimum of a simulation, so it is performed first. */
    InitDomainHotHalo(d, grid);

    /* Input data for clouds initialization */
#if CLOUDS != NO

    /* TODO: Cloud initialization should not take place if this is a restart! */
    /* TODO: Differentiate between cloud initialization and all other things, e.g. cloud analysis, sfr, etc */
    InputDataClouds(d, grid);

#endif

    /* Initialize nozzle. */
#if (NOZZLE != NONE) && (NOZZLE_FILL == NF_PRIMITIVE)
    InitDomainNozzle(d, grid);

    /* Measure Nozzle volume by counting cells. This is the volume into which
     * mass, momentum, energy are dumped. We therefore use DOM_LOOP */
#elif NOZZLE && (NOZZLE_FILL == NF_CONSERVATIVE)
    NozzleVolume(d, grid);

#endif

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

    /* Modules -- Conservative source terms */

    // TODO: Maybe even move this to boundaries
    // TODO: Create a buffer data array in which changes are stored, whilst passing original d into functions

    /* Conservative NOZZLE_FILL Dump */
    // TODO: Move to separate function in outflows.c
    // NOTE: When modularizing, the nozzle filling should maybe be done in either rk_step or ctu_step, or in a rhs routine.

#if (NOZZLE != NONE) && (NOZZLE_FILL == NF_CONSERVATIVE)
    NozzleFill(d, grid);

#endif


    /* Star-formation */
//    WarmPhaseStarFormationRate(d, grid);

    /* Supernova feedback */


    /* Analysis */

    /* Warm phase statistics */
#if CLOUD_OUTPUT
    CloudAnalysis(d, grid);
#endif

    /* Accretion */
#if ACCRETION == YES
    SphericalAccretion(d, grid);

#endif


    /* Analysis output */

#if ACCRETION

#if ACCRETION_OUTPUT == YES
    SphericalAccretionOutput();

#endif

#if BONDI_ACCRETION_OUTPUT == YES
    BondiAccretionOutput();

#endif

#endif // if ACCRETION == YES

#if OUTFLOW_OUTPUT == YES
    OutflowOutput();

#endif

#if CLOUD_OUTPUT == YES
    CloudOutput();

#endif

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

    // TODO: After all the modules are isolated, move the below variables into the respective module sections
    int i, j, k, nv;
    double *x1, *x2, *x3;
    double vc;
    double out_primitives[NVAR], halo_primitives[NVAR], result[NVAR], mirror[NVAR];
    double * prim_buffer = mirror;
    static int touch = 0;
#if ACCRETION == YES
    double ****Vc_new;
#endif
#if CLOUD_EXTRACT_CENTRAL_BUFFER == YES
    static int once_clear_nozzle_surrounding = 0;
#endif

    /* These are the geometrical central points */
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];
    x3 = grid->x[KDIR];

#if INTERNAL_BOUNDARY == YES
    if (side == 0) {    /* -- check solution inside domain -- */


#if SUPERNOVAE

        /* Determine number and locations of SN to explode for this timestep */
        StrewSupernovae();

        if (sn.num_dt > 0) {

            /* Find maximum cell size in this domain to determine injection radius. */
            /* This might be better done in initalize.c, and saved in the grid struct.
             * Although it may change every timestep if AMR is used? */
            // TODO: put this into initialize (overwrite initialize.c)
            double dx_max = FindDxMax(grid);
            double sn_inj = SUPERNOVA_RCELLS * dx_max;

            TOT_LOOP(k, j, i) {

                        /* Inject supernova energy cellwise */
                        NVAR_LOOP(nv) result[nv] = d->Vc[nv][k][j][i];
                        InjectSupernovae(result, x1[i], x2[j], x3[k], sn_inj);
                        NVAR_LOOP(nv) d->Vc[nv][k][j][i] = result[nv];
                    }
        }

#endif   // If SUPERNOVAE



        /* Accretion and Nozzle */

            // TODO: Separate ACCRETION from NOZZLE
            // TODO: Create a AGN switch, and change NOZZLE to AGN_NOZZLE
            // TODO: Change ACCRETION -> AGN_ACCRETION
#if ACCRETION == YES

            // TODO: Perhaps make this available more generally
            /* Create buffer array for internal sink region solution */
            Vc_new = ARRAY_4D(NVAR, NX3_TOT, NX2_TOT, NX1_TOT, double);

            /* Update outflow state according to accretion */
#if FEEDBACK_CYCLE == YES

            int is_on = os.is_on;

            /* Set outflow stat according to measured accretion rate */
            SetOutflowState(&os);

            /* Call SetNozzleGeometry to adjust nozzle radius */
            SetNozzleGeometry(&nz);

            /* If the nozzle switches on, drastically reduce timestep */
#if NOZZLE_FILL == NF_PRIMITIVE
            if (is_on < os.is_on) g_dt = FBC_NOZZLE_DT;
//            if (is_on < os.is_on) g_dt *= pow(10, 0.5 * log10(first_dt / g_dt)))
#endif

#endif  // FEEDBACK_CYCLE

#endif  // ACCRETION


#if CLOUD_EXTRACT_CENTRAL_BUFFER == YES
        TOT_LOOP(k, j, i) {

                    /* Do this before anything else on the nozzle - this clears the surrounding
                     * and is useful when switching on the jet midway through a simulation after a restart.
                     * This bit will also occur before the first time step of a fresh simulation, but
                     * the effect is equivalent to CloudExtractCentralBuffer, so it has no effect in that case.
                     * (We don't have information here whether this is a restart or not.)
                     */
                    if (once_clear_nozzle_surrounding < 2) {

                        HotHaloPrimitives(halo_primitives, x1[i], x2[j], x3[k]);
                        VAR_LOOP(nv) result[nv] = d->Vc[nv][k][j][i];

                        ClearNozzleSurrounding(result, halo_primitives, x1[i], x2[j], x3[k]);
                        VAR_LOOP(nv) d->Vc[nv][k][j][i] = result[nv];

                    }
                }

        once_clear_nozzle_surrounding ++;
#endif


        TOT_LOOP(k, j, i) {

#if (NOZZLE == NO) || (NOZZLE_FILL == NF_CONSERVATIVE)
                        if (0) {}
#else
                        if (InNozzleRegion(x1[i], x2[j], x3[k])) {

                            OutflowPrimitives(out_primitives, x1[i], x2[j], x3[k]);

                            NVAR_LOOP(nv) {
                                vc = d->Vc[nv][k][j][i];
                                d->Vc[nv][k][j][i] = vc + (out_primitives[nv] - vc) *
                                                          Profile(x1[i], x2[j], x3[k]);
                            }

                            d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

                        } // InNozzleRegion
#endif   // If NOZZLE


#if ACCRETION == YES
                            else if (InSinkRegion(x1[i], x2[j], x3[k])) {

#if SINK_METHOD == SINK_FREEFLOW

                                /* Apply free flow internal boundary conditions */
                                SphericalFreeflowInternalBoundary(d->Vc, i, j, k, x1, x2, x3, result);

#elif SINK_METHOD == SINK_VACUUM

                                /* Apply Vacuum internal boundary conditions */
                                VacuumInternalBoundary(result);

#elif SINK_METHOD == SINK_BONDI

                                // TODO: Change this to: extract mass according to Bondi accretion rate.
                                /* Apply Bondi flow solution */
                                BondiFlowInternalBoundary(x1[i], x2[j], x3[k], result);

#elif SINK_METHOD == SINK_FEDERRATH

                                /* Remove mass according to Federrath's sink particle method */
                                FederrathSinkInternalBoundary(d->Vc, i, j, k, x1, x2, x3, grid->dV, result);

// TODO: Create HOT_HALO sink method:just use hot halo fixed boundaries.

#endif  // SINK_METHOD

                                /* Copy results to temporary solution array */
                                NVAR_LOOP(nv) Vc_new[nv][k][j][i] = result[nv];

                                d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

                            } // InSinkRegion

#else  // If no ACCRETION
                        else if (InFlankRegion(x1[i], x2[j], x3[k])) {

                            HotHaloPrimitives(halo_primitives, x1[i], x2[j], x3[k]);

                            NVAR_LOOP(nv) d->Vc[nv][k][j][i] = halo_primitives[nv];
                            d->flag[k][j][i] |= FLAG_INTERNAL_BOUNDARY;

                        } // InFlankRegion
#endif // If ACCRETION


                    } // TOT_LOOP

#if ACCRETION == YES

            /* Copy solution over to Vc array */
            TOT_LOOP(k, j, i) {

                        if (InSinkRegion(x1[i], x2[j], x3[k])) {
                            NVAR_LOOP(nv) d->Vc[nv][k][j][i] = Vc_new[nv][k][j][i];
                        }

                    } // Update TOT_LOOP

            FreeArray4D((void *) Vc_new);

#endif   // If Accretion


    } // side == 0

#endif   // If INTERNAL_BOUNDARY


    if (side == FLOWAXIS(X2_BEG, X3_BEG, X1_BEG)) {
        if (box->vpos == CENTER) {
            BOX_LOOP(box, k, j, i) {

                        HaloOuterBoundary(side, d, i, j, k, grid, &touch);

                    }
        } else if (box->vpos == X1FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X2FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X3FACE) {
            BOX_LOOP(box, k, j, i) { }
        }
    }

    if (side == X1_END) {  /* -- X1_END boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box, k, j, i) {

                        HaloOuterBoundary(side, d, i, j, k, grid, &touch);

                    }
        } else if (box->vpos == X1FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X2FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X3FACE) {
            BOX_LOOP(box, k, j, i) { }
        }
    }


    if (side == FLOWAXIS(X3_BEG, X1_BEG, X2_BEG)) {
        if (box->vpos == CENTER) {
            BOX_LOOP(box, k, j, i) {

                        HaloOuterBoundary(side, d, i, j, k, grid, &touch);

                    }
        } else if (box->vpos == X1FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X2FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X3FACE) {
            BOX_LOOP(box, k, j, i) { }
        }
    }


    if (side == X2_END) {  /* -- X2_END boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box, k, j, i) {

                        HaloOuterBoundary(side, d, i, j, k, grid, &touch);

                    }
        } else if (box->vpos == X1FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X2FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X3FACE) {
            BOX_LOOP(box, k, j, i) { }
        }
    }


    /* This side is where the nozzle is located and the outflow emerges,
     * unless internal boundary is set. */
    if (side == FLOWAXIS(X1_BEG, X2_BEG, X3_BEG)) {
        if (box->vpos == CENTER) {
            BOX_LOOP(box, k, j, i) {

                        /* If we're doing a full galaxy simulation,
                           as opposed to a half-galaxy simulation,
                           treat just like other _BEG boundaries    */

                        if (nz.is_two_sided && GEOMETRY != SPHERICAL) {
                            HaloOuterBoundary(side, d, i, j, k, grid, &touch);

                        } // full galaxy

                        else { // half galaxy

                            /* Reflective boundary */
                            NVAR_LOOP(nv) {
                                mirror[nv] = d->Vc[nv]
                                FLOWAXIS([k][j][2 * IBEG - i - 1],
                                         [k][2 * JBEG - j - 1][i],
                                         [2 * KBEG - k - 1][j][i]);
                            }

                            mirror[FLOWAXIS(VX1, VX2, VX3)] *= -1.0;

#if (NOZZLE == NO) || (NOZZLE_FILL == NF_CONSERVATIVE)
                            NVAR_LOOP(nv) d->Vc[nv][k][j][i] = mirror[nv];
#else

                            if (InNozzleRegion(x1[i], x2[j], x3[k])) {

                                OutflowPrimitives(out_primitives, x1[i], x2[j], x3[k]);

                                NVAR_LOOP(nv) {
                                    d->Vc[nv][k][j][i] = mirror[nv] + (out_primitives[nv] - mirror[nv]) *
                                                                      Profile(x1[i], x2[j], x3[k]);
                                }

                            }
                            else {
                                NVAR_LOOP(nv) {
                                    d->Vc[nv][k][j][i] = mirror[nv];
                                }
                            }
#endif

                        } // half galaxy

                    } // BOX_LOOP
        } else if (box->vpos == X1FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X2FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X3FACE) {
            BOX_LOOP(box, k, j, i) { }
        }
    }

    if (side == X3_END) {  /* -- X3_END boundary -- */
        if (box->vpos == CENTER) {
            BOX_LOOP(box, k, j, i) {

                        HaloOuterBoundary(side, d, i, j, k, grid, &touch);

                    }
        } else if (box->vpos == X1FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X2FACE) {
            BOX_LOOP(box, k, j, i) { }
        } else if (box->vpos == X3FACE) {
            BOX_LOOP(box, k, j, i) { }
        }
    }
}


// TODO: Use vector only for relativistic simulations

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

    /* Common variables */
    double r, gr;
    gr = 0;
#if GRAV_POTENTIAL == GRAV_2D_TABLE
    double fz, sz, gz;
#endif

    /* Types of galactic gravitational potential */
#if GRAV_POTENTIAL == GRAV_HERNQUIST
    /* Hernquist potential */

    double a, rho0;
    r = SPH1(x1, x2, x3);
    a = g_inputParam[PAR_HRAD] * ini_code[PAR_HRAD];
    rho0 = g_inputParam[PAR_HRHO] * ini_code[PAR_HRHO];

    gr = -2. * CONST_PI * CONST_G / vn.newton_norm * rho0 * a / pow((1. + r / a), 2);


    /* Gravity Table */
#elif GRAV_POTENTIAL == GRAV_TABLE || GRAV_POTENTIAL == GRAV_2D_TABLE

#if GRAV_POTENTIAL == GRAV_2D_TABLE
    double fr;
    fr = MIN(MAX(gr_r[0], fabs(POL1(x1, x2, x3))), gr_r[gr_nr-1]);
    fz = MIN(MAX(gr_z[0], fabs(POL3(x1, x2, x3))), gr_z[gr_nz-1]);
    sz = SGN(POL3(x1, x2, x3))

    gr = InterpolationWrapper2D(gr_r, gr_z, gr_acc_r, gr_nr, gr_nz, fr, fz);
    gz = sz * InterpolationWrapper2D(gr_r, gr_z, gr_acc_z, gr_nr, gr_nz, fr, fz);
#else
    r = MIN(MAX(gr_r[0], SPH1(x1, x2, x3)), gr_r[gr_nr-1]);
    gr = InterpolationWrapper(gr_r, gr_acc_r, gr_nr, r);
#endif

#endif /* Types of galactic gravitational potential */

    /* Add potential of point mass at center */
#if CENTRAL_BH == YES

    double mbh;
#if ACCRETION == YES
    mbh = ac.mbh;
#else
    mbh = g_inputParam[PAR_AMBH] * ini_code[PAR_AMBH];
#endif

    /* Softening length */
    double soft;
    soft = MAX(g_inputParam[PAR_ASNK] * ini_code[PAR_ASNK], g_inputParam[PAR_OSPH] * ini_code[PAR_OSPH]);

    r = SPH1(x1, x2, x3);
    gr -= CONST_G * mbh * r / (vn.newton_norm * pow(r * r + soft * soft, 1.5));

#endif


/* Cylindrical potential, if 2D */
/* TODO: Change to #if GRAV_2D_POTENTIAL == TRUE to accommodate 2D analytic potentials */
#if GRAV_POTENTIAL == GRAV_2D_TABLE

    double px1, px2, px3;
    px1 = POL1(x1, x2, x3);
    px2 = POL2(x1, x2, x3);
    px3 = POL3(x1, x2, x3);

    /* Gravity pointing to (0,0,0) - possibly reconsider */
    g[IDIR] = VPOL_1(px1, px2, px3, gr, 0, gz);
    g[JDIR] = VPOL_2(px1, px2, px3, gr, 0, gz);
    g[KDIR] = VPOL_3(px1, px2, px3, gr ,0, gz);

/* Else, spherical acceleration */
#else
    double sx1, sx2, sx3;
    sx1 = SPH1(x1, x2, x3);
    sx2 = SPH2(x1, x2, x3);
    sx3 = SPH3(x1, x2, x3);

    /* Gravity pointing to (0,0,0) - possibly reconsider */
    g[IDIR] = VSPH_1(sx1, sx2, sx3, gr, 0, 0);
    g[JDIR] = VSPH_2(sx1, sx2, sx3, gr, 0, 0);
    g[KDIR] = VSPH_3(sx1, sx2, sx3, gr, 0, 0);
#endif

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

    /* Common variables */
    double r, pot;

    /* Initialize */
    pot = 0;

#if GRAV_POTENTIAL == GRAV_2D_TABLE
    double z;
#endif


    /* Types of galactic gravitational potentials */
#if GRAV_POTENTIAL == GRAV_HERNQUIST

    /* Hernquist potential, given a "galaxy" mass. */
    double a, mg;

    r = SPH1(x1, x2, x3);
    a = g_inputParam[PAR_HRAD] * ini_code[PAR_HRAD];
    mg = g_inputParam[PAR_MGAL] * ini_code[PAR_MGAL];

    pot = - mg * CONST_G / vn.newton_norm / (r + a);

    /* Gravity table */
#elif GRAV_POTENTIAL == GRAV_TABLE || GRAV_POTENTIAL == GRAV_2D_TABLE

    r = MIN(MAX(gr_r[0], SPH1(x1, x2, x3)), gr_r[gr_nr-1]);
#if GRAV_POTENTIAL == GRAV_2D_TABLE
    r = MIN(MAX(gr_r[0], fabs(POL1(x1, x2, x3))), gr_r[gr_nr-1]);
    z = MIN(MAX(gr_z[0], fabs(POL3(x1, x2, x3))), gr_z[gr_nz-1]);
    pot = InterpolationWrapper2D(gr_r, gr_z, gr_phi, gr_nr, gr_nz, r, z);
#else
    pot = InterpolationWrapper(gr_r, gr_phi, gr_nr, r);
#endif

#endif /* Types of galactic gravitational potentials */


    /* Add potential of point mass at center */
#if CENTRAL_BH == YES

    double mbh;
#if ACCRETION == YES
    mbh = ac.mbh;
#else
    mbh = g_inputParam[PAR_AMBH] * ini_code[PAR_AMBH];
#endif

    /* Softening length */
    double soft;
    soft = MAX(g_inputParam[PAR_ASNK] * ini_code[PAR_ASNK], g_inputParam[PAR_OSPH] * ini_code[PAR_OSPH]);

    r = SPH1(x1, x2, x3);
    pot += -CONST_G * mbh / (vn.newton_norm * sqrt(r * r + soft * soft));

#endif

    return pot;

}

#endif

