/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Loop on the computational cells to assign initial conditions.

  This function is called anyway, even if restart from file is enabled.
  This is useful to initialized a number of global variables and/or 
  user-defined parameters by calling Init().

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include "pluto_usr.h"

#if CLOUDS_MULTI == YES
#include "read_hot_table.h"
#include "multicloud_init.h"
#endif


/* ********************************************************************* */
void Startup(Data *d, Grid *G)
/*! 
 *
 *
 *
 *
 *********************************************************************** */
{
    int i, j, k;
    int isub, jsub, ksub, nsub = 5;
    int nv, l_convert;
    static double **ucons, **uprim;
    double x1, x2, x3;
    double x1p, x2p, x3p;
    double x1s, x2s, x3s;
    double dx1, dx2, dx3;
    double us[256], u_av[256], b[3];
    double scrh;
    struct GRID *GX, *GY, *GZ;


#if CLOUDS_MULTI == YES
    double *x1cld, *x2cld, *x3cld, *v1cld, *v2cld, *v3cld;
    int *cldnum;
    char **cldfile_list, temp[100];
    int count, nclouds, nfiles, nc;
    int oncefilelist, oncecloudlist;
    struct InGrid *Gin;
    struct cld_domain cld;
    double x1beg, x1end, x2beg, x2end, x3beg, x3end;
    double x1cl, x1ch, x2cl, x2ch, x3cl, x3ch, delx1cin, delx2cin, delx3cin;
    int ibeg, iend, jbeg, jend, kbeg, kend;

    Gin = (struct InGrid *) malloc(sizeof(struct InGrid));
    //cld=(struct cld_domain) malloc(sizeof(struct cld_domain));
#endif


    GX = G;
    GY = G + 1;
    GZ = G + 2;

    if (uprim == NULL) {
        uprim = ARRAY_2D(NMAX_POINT, NVAR, double);
        ucons = ARRAY_2D(NMAX_POINT, NVAR, double);
    }

/* ----  set labels  ---- */

    EXPAND(MXn = VXn = VX1;,
           MXt = VXt = VX2;,
           MXb = VXb = VX3;)

#if PHYSICS == MHD || PHYSICS == RMHD
    EXPAND(BXn = BX1;  ,
           BXt = BX2;  ,
           BXb = BX3;)
#endif

/* ----------------------------------------------------------
    If the -input-data-file command line argument is given, 
    try to read initial conditions from external files.
   ---------------------------------------------------------- */

    print1("> Assigning initial conditions (Startup) ...\n");

/* --------------------------------------------------------------
                    Assign initial conditions   
   -------------------------------------------------------------- */

    KTOT_LOOP(k) {
        JTOT_LOOP(j) {
            ITOT_LOOP(i) {

#if GEOMETRY == CYLINDRICAL
                x1 = GX->xgc[i]; x1p = GX->xr[i]; dx1 = GX->dx[i];
                x2 = GY->xgc[j]; x2p = GY->xr[j]; dx2 = GY->dx[j];
                x3 = GZ->xgc[k]; x3p = GZ->xr[k]; dx3 = GZ->dx[k];
#else
                x1 = GX->x[i];
                x1p = GX->xr[i];
                dx1 = GX->dx[i];
                x2 = GY->x[j];
                x2p = GY->xr[j];
                dx2 = GY->dx[j];
                x3 = GZ->x[k];
                x3p = GZ->xr[k];
                dx3 = GZ->dx[k];
#endif

                for (nv = NVAR; nv--;) d->Vc[nv][k][j][i] = u_av[nv] = 0.0;

/*  ----------------------------------------------------------------
                Compute volume averages      
    ---------------------------------------------------------------- */

#ifdef PSI_GLM
                u_av[PSI_GLM] = us[PSI_GLM] = 0.0;
#endif

#if INITIAL_SMOOTHING == YES

                for (ksub = 0; ksub < nsub; ksub++){
                for (jsub = 0; jsub < nsub; jsub++){
                for (isub = 0; isub < nsub; isub++){

                  x1s = x1 + (double)(1.0 - nsub + 2.0*isub)/(double)(2.0*nsub)*dx1;
                  x2s = x2 + (double)(1.0 - nsub + 2.0*jsub)/(double)(2.0*nsub)*dx2;
                  x3s = x3 + (double)(1.0 - nsub + 2.0*ksub)/(double)(2.0*nsub)*dx3;

                  Init (us, x1s, x2s, x3s);
                  for (nv = 0; nv < NVAR; nv++) {
                    u_av[nv] += us[nv]/(double)(nsub*nsub*nsub);
                  }
                }}}

#else

                Init(u_av, x1, x2, x3);

#endif

                for (nv = NVAR; nv--;) d->Vc[nv][k][j][i] = u_av[nv];

/* -----------------------------------------------------
        Initialize cell-centered vector potential 
        (only for output purposes)
   ----------------------------------------------------- */

#if PHYSICS == MHD || PHYSICS == RMHD
#if UPDATE_VECTOR_POTENTIAL == YES
                D_EXPAND(                     ,
                  d->Ax3[k][j][i] = u_av[AX3];  ,
                  d->Ax1[k][j][i] = u_av[AX1];
                  d->Ax2[k][j][i] = u_av[AX2];)
#endif
#endif

                /* -------------------------------------------------------------
                    Assign staggered components;
                    If a vector potential is used (ASSIGN_VECTOR_POTENTIAL == YES),
                    use the STAGGERED_INIT routine;
                    otherwise assign staggered components directly from
                    the init.c and ignore the vector potential.

                    NOTE: in N dimensions only N components are assigned
                            through this call.
                   ------------------------------------------------------------- */

#if PHYSICS == MHD || PHYSICS == RMHD
#if ASSIGN_VECTOR_POTENTIAL == YES
                VectorPotentialDiff(b, i, j, k, G);

#ifdef STAGGERED_MHD
                 for (nv = 0; nv < DIMENSIONS; nv++) {
                   d->Vs[nv][k][j][i] = b[nv];
                 }
#else
                 for (nv = 0; nv < DIMENSIONS; nv++) {
                   d->Vc[BX+nv][k][j][i] = b[nv];
                 }
#endif

#else

#ifdef STAGGERED_MHD
                 D_EXPAND(
                   Init (u_av, x1p, x2, x3);
                   d->Vs[BX1s][k][j][i] = u_av[BX1];       ,

                   Init (u_av, x1, x2p, x3);
                   d->Vs[BX2s][k][j][i] = u_av[BX2];       ,

                   Init (u_av, x1, x2, x3p);
                   d->Vs[BX3s][k][j][i] = u_av[BX3];
                 )
#endif
#endif  /* ASSIGN_VECTOR_POTENTIAL */
#endif /* PHYSICS == MHD || PHYSICS == RMHD */

            }
        }
    }

#if CLOUDS_MULTI == YES

    oncefilelist = 0;
    oncecloudlist = 0;

    // Read cloud input data file list //
    if (oncecloudlist == 0) {
        FILE *fp_cldfilelist;
        fp_cldfilelist = fopen("cloud_filelist.dat", "r");
        if (fp_cldfilelist == NULL) {
            print1("Input cloud_filelist.dat not found \n");
            QUIT_PLUTO(1);
        }
        fscanf(fp_cldfilelist, "%d \n", &nfiles);
        cldfile_list = ARRAY_2D(nfiles, 100, char);
        count = 0;

        while (1) {
            fscanf(fp_cldfilelist, "%s \n", temp);
            sprintf(cldfile_list[count], "%s", temp);
            count++;
            if (feof(fp_cldfilelist)) break;
        }
        fclose(fp_cldfilelist);
        oncecloudlist = 1;
    }//--nfiles==NULL


//-----Generate cloud positions and store in arrays-----//
    nclouds = g_inputParam[PAR_NCLD];


#ifdef PARALLEL
    if (prank == 0){
        gen_cldlist(nclouds, nfiles);
        MPI_Barrier(MPI_COMM_WORLD);
    }
#else
    gen_cldlist(nclouds, nfiles);
#endif



    // Read cloud position list & number of clouds //
    if (oncefilelist == 0) {
        FILE *fp_cldposlist;
        fp_cldposlist = fopen("cloud_pos_list.dat", "r");
        if (fp_cldposlist == NULL) {
            print1("Input cloud_pos_list.dat not found \n");
            QUIT_PLUTO(1);
        }
//	fscanf(fp_cldposlist,"%d \n",&nclouds);

        x1cld = ARRAY_1D(nclouds, double);
        x2cld = ARRAY_1D(nclouds, double);
        x3cld = ARRAY_1D(nclouds, double);
        cldnum = ARRAY_1D(nclouds, int);
        v1cld = ARRAY_1D(nclouds, double);
        v2cld = ARRAY_1D(nclouds, double);
        v3cld = ARRAY_1D(nclouds, double);

        count = 0;

        while (1) {
            fscanf(fp_cldposlist, "%d %lf %lf %lf %lf %lf %lf \n", &cldnum[count], &x1cld[count], &x2cld[count],
                   &x3cld[count], &v1cld[count], &v2cld[count], &v3cld[count]);
            count++;
            if (feof(fp_cldposlist)) break;
        }
        fclose(fp_cldposlist);
        oncefilelist = 1;
    }//--nclouds==NULL




    //---Define domain boundaries---//
    x1beg = GX->x[0];
    x1end = GX->x[NX1_TOT - 1];
    x2beg = GY->x[0];
    x2end = GY->x[NX2_TOT - 1];
    x3beg = GZ->x[0];
    x3end = GZ->x[NX3_TOT - 1];

    //----Read input grid file----//
    readgridfile(Gin);

    //----Define extent of input cloud----//
    //----Assumed cloud input centred at origin---//
    x1cl = Gin->x1[0];
    x1ch = Gin->x1[Gin->id_nx1 - 1];
    x2cl = Gin->x2[0];
    x2ch = Gin->x2[Gin->id_nx2 - 1];
    x3cl = Gin->x3[0];
    x3ch = Gin->x3[Gin->id_nx3 - 1];

    //---Store half width of domain extent of input cloud grid---//
    delx1cin = (x1ch - x1cl) * 0.5;
    delx2cin = (x2ch - x2cl) * 0.5;
    delx3cin = (x3ch - x3cl) * 0.5;

    print1("Starting multi cloud initialisation \n");
    print1("Extent of input Grid (grid_in.out) \n");
    print1("X1: %lf to %lf \n", x1cl, x1ch);
    print1("X2: %lf to %lf \n", x2cl, x2ch);
    print1("X3: %lf to %lf \n", x3cl, x3ch);

    //----Loop over cloud position list to check for clouds in domain----//
    for (nc = 0; nc < nclouds; nc++) {
        if (((x1cld[nc] - delx1cin >= x1beg) && (x1cld[nc] - delx1cin <= x1end))
            || ((x1cld[nc] + delx1cin >= x1beg) && (x1cld[nc] + delx1cin <= x1end))) {
            if (((x2cld[nc] - delx2cin >= x2beg) && (x2cld[nc] - delx2cin <= x2end))
                || ((x2cld[nc] + delx2cin >= x2beg) && (x2cld[nc] + delx2cin <= x2end))) {
                if (((x3cld[nc] - delx3cin >= x3beg) && (x3cld[nc] - delx3cin <= x3end))
                    || ((x3cld[nc] + delx3cin >= x3beg) && (x3cld[nc] + delx3cin <= x3end))) {

                    /*----Store cloud centres in cld_domain structure----//
                      Stored values are cloud centres relative to centre of
                      input grid, for ease in interpolation -------------*/

                    cld.x1c = x1cld[nc] - (x1ch + x1cl) * 0.5;
                    cld.x2c = x2cld[nc] - (x2ch + x2cl) * 0.5;
                    cld.x3c = x3cld[nc] - (x3ch + x3cl) * 0.5;

                    cld.v1 = v1cld[nc];
                    cld.v2 = v2cld[nc];
                    cld.v3 = v3cld[nc];

                    //---Read in cloud density---//
                    Read_Multicld(cldfile_list[cldnum[nc]]);


                    //----Locate array indices of cloud domain---//
                    ibeg = (x1cld[nc] - delx1cin < x1beg) ? 0 : locate(GX->x, x1cld[nc] - delx1cin, NX1_TOT);
                    iend = (x1cld[nc] + delx1cin > x1end) ? NX1_TOT - 1 : locate(GX->x, x1cld[nc] + delx1cin, NX1_TOT);
                    jbeg = (x2cld[nc] - delx2cin < x2beg) ? 0 : locate(GY->x, x2cld[nc] - delx2cin, NX2_TOT);
                    jend = (x2cld[nc] + delx2cin > x2end) ? NX2_TOT - 1 : locate(GY->x, x2cld[nc] + delx2cin, NX2_TOT);
                    kbeg = (x3cld[nc] - delx3cin < x3beg) ? 0 : locate(GZ->x, x3cld[nc] - delx3cin, NX3_TOT);
                    kend = (x3cld[nc] + delx3cin > x3end) ? NX3_TOT - 1 : locate(GZ->x, x3cld[nc] + delx3cin, NX3_TOT);


                    for (i = ibeg + 1; i <= iend; i++) {
                        for (j = jbeg + 1; j <= jend; j++) {
                            for (k = kbeg + 1; k <= kend; k++) {
                                x1 = GX->x[i];
                                x2 = GY->x[j];
                                x3 = GZ->x[k];

                                //----Initialise cloud cube---//
                                for (nv = NVAR; nv--;) u_av[nv] = 0.0;
                                Init_multiclds(u_av, x1, x2, x3, cld);

                                for (nv = NVAR; nv--;) d->Vc[nv][k][j][i] = u_av[nv];

                            }
                        }
                    }


                }
            }
        }//--if in cld domain
    }//for nc


#endif


#ifdef STAGGERED_MHD

    /* ---------------------------------------------------
         Compute cell-centered magnetic field
         by simple arithmetic averaging. This
         is useful only for saving the first
         output, since the average will be
         re-computed at the beginning of the computation.
       --------------------------------------------------- */

     DOM_LOOP(k,j,i){
       D_EXPAND(
         d->Vc[BX1][k][j][i] = 0.5*(d->Vs[BX1s][k][j][i] + d->Vs[BX1s][k][j][i-1]); ,
         d->Vc[BX2][k][j][i] = 0.5*(d->Vs[BX2s][k][j][i] + d->Vs[BX2s][k][j-1][i]); ,
         d->Vc[BX3][k][j][i] = 0.5*(d->Vs[BX3s][k][j][i] + d->Vs[BX3s][k-1][j][i]);
       )
     }

#endif

/* --------------------------------------------------------------------
         Check if values have physical meaning...
   -------------------------------------------------------------------- */

#if PHYSICS != ADVECTION
    KDOM_LOOP(k) {
        x3 = GZ->x[k];
        JDOM_LOOP(j) {
            x2 = GY->x[j];
            IDOM_LOOP(i) {
                x1 = GX->x[i];

                for (nv = NVAR; nv--;) us[nv] = d->Vc[nv][k][j][i];

                if (us[RHO] <= 0.0) {
                    print("! Startup: density is negative, zone [%f, %f, %f]\n", x1, x2, x3);
                    QUIT_PLUTO(1);
                }
#if HAVE_ENERGY
                if (us[PRS] <= 0.0) {
                    print("! Startup: pressure is negative, zone [%f, %f, %f]\n", x1, x2, x3);
                    QUIT_PLUTO(1);
                }
#endif
#if (PHYSICS == RHD && USE_FOUR_VELOCITY == NO) || PHYSICS == RMHD
                scrh = EXPAND(us[VX1]*us[VX1], + us[VX2]*us[VX2], + us[VX3]*us[VX3]);
                if (scrh >= 1.0){
                  print ("! Startup: total velocity exceeds 1\n");
                  QUIT_PLUTO(1);
                }
#endif
            }
        }
    }
#endif

/* --------------------------------------------------------------------
     Convert primitive variables to conservative ones
   -------------------------------------------------------------------- */

    Boundary(d, -1, G);

/* --------------------------------------------------------------------
                    Convert to conservative
   -------------------------------------------------------------------- */
/*
  KDOM_LOOP(k) {
  JDOM_LOOP(j){
    IDOM_LOOP(i) {
    for (nv = NVAR; nv--;  ) {
      uprim[i][nv] = d->Vc[nv][k][j][i];
    }}
    PrimToCons(uprim, d->Uc[k][j], IBEG, IEND);
  }}
*/

#if CLOUDS_MULTI == YES
    free((struct InGrid *) Gin);
#endif

}



