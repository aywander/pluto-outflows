#include "pluto.h"
#include "pluto_usr.h"
#include "init_tools.h"
#include "multicloud_init.h"
#include "random.h"


#if CLOUDS_MULTI == YES
/* ********************************************************************* */
void Init_multiclds (double *v, double x1, double x2, double x3, struct cld_domain cld)
//---Same as in init.c but called for multi clouds set up---//
{

    int nv;
    double halo_primitives[NVAR], out_primitives[NVAR];
    double cloud_primitives[NVAR];



    /* Initialize nozzle if we're in hemisphere around
     * nozzle inlet region, otherwise halo */


    if (InNozzleRegion(x1, x2, x3) && (g_inputParam[FLAG_JET] > 0)) {

        OutflowPrimitives(out_primitives, x1, x2, x3);
        HotHaloPrimitives(halo_primitives, x1, x2, x3);

        NVAR_LOOP(nv) {
            v[nv] = halo_primitives[nv] +
                    (out_primitives[nv] - halo_primitives[nv]) * Profile_cap(x1, x2, x3);
        }
    }

        /* Initialize halo. Hot, and warm, if included */
    else {

        /* First get primitives array for hot halo */
        HotHaloPrimitives(halo_primitives, x1, x2, x3);


        /* If we're in the domain of the clouds cube */
        if (MultiCloudPrimitives(cloud_primitives, x1, x2, x3, cld)) {
            NVAR_LOOP(nv) v[nv] = cloud_primitives[nv];
        }
            /* If not a cloud pixel then use hot halo primitives*/
        else {
            NVAR_LOOP(nv) v[nv] = halo_primitives[nv];
        }


    }

#if INTERNAL_BOUNDARY == YES
    if (InNozzleSphere(x1, x2, x3) && (g_inputParam[FLAG_JET] > 0)) {
        HotHaloPrimitives(halo_primitives, x1, x2, x3);
        NVAR_LOOP(nv) {
            v[nv] = halo_primitives[nv];
        }
    }
#endif


#if PHYSICS == MHD || PHYSICS == RMHD

    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = 0.0;

    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = 0.0;

#endif
}



/* ************************************************************** */
int MultiCloudPrimitives(double* cloud, 
                    const double x1, const double x2, const double x3, struct cld_domain cld)
/* 
 * Same as in init_tools.c but for multi cloud set up.
 **************************************************************** */
{

    double halo[NVAR], vel[COMPONENTS], scrh;
    int nv, cube_pixel[DIMENSIONS];
    int is_cloud = 0;



    /* Get the fractal factor for this cell. cld.x#c are cloud centres relative to input domain, 
	passed from startup.c  */
    GetFractalData(cloud, x1 - cld.x1c, x2 - cld.x2c, x3 - cld.x3c);

    /* Apodize cloud with a background mean density profile */
    CloudApodize(cloud, x1, x2, x3);



    /* Extract w.r.t. hot halo. Halo primitives are required 
     * for this step. */
    HotHaloPrimitives(halo, x1, x2, x3);

    /* Calculate cloud primitives so that we can get temperature */

    /* Cloud velocity */
    cloud[VX1] = cld.v1;
    cloud[VX2] = cld.v2;
    cloud[VX3] = cld.v3;
    CloudVelocity(cloud, halo, x1, x2, x3);

    /* Cloud pressure
     * Underpressure the clouds slightly, so that they
     * don't emit sound waves.
     *   Note that the -xnjet option is available and that
     * the factor here should be compared the pressure gradient
     * calculated in SET_JET_DOMAIN jet_domain.c
     * */
    cloud[PRS] = halo[PRS] * CLOUD_UNDERPRESSURE;

    /* Tracers */
    cloud[TRC] = 0.0;
    cloud[TRC + 1] = 1.0;

    /* Final test - is cloud pixel thermally stable? */
    is_cloud = WarmTcrit(cloud);



    /* Fill cloud array with halo primitves if not a cloud cell. This is not
     * strictly necessary, since it is done outside CloudPrimitives, but we
     * do it anyway for completeness. */
    if (is_cloud == 0) { NVAR_LOOP(nv) cloud[nv] = halo[nv]; }

    return is_cloud;
}

/* ************************************************************** */
void Read_Multicld(char *fname)
/*!
 * Read in external density from a input filename. 
 *
 **************************************************************** */
{
  int get_var[] = {RHO, -1};

  /* Read cloud data from external file */
    InputDataSet("./grid_in.out", get_var);
    InputDataRead(fname, CUBE_ENDIANNESS);

}





/* ***********************************************************************/
void readgridfile (struct InGrid *Gin)
/* 
 * Read in external input grid from grid_in.out and store coordinate axes &
 * number of elements in structure InGrid
 *
 *********************************************************************** */
{
    int i, ip, nv, success;
    size_t dsize, dcount;
    char *sub_str, sline[256];
    const char delimiters[] = " \t\r\f\n";
    double xl, xr;
    fpos_t file_pos;
    int id_nx1, id_nx2, id_nx3;
    double *temp;
    FILE *fp;

    fp = fopen("grid_in.out", "r");
    if (fp == NULL) {
        printf("! InputDataSet: grid file not found\n");
        exit(1);
    }
/* --------------------------------------------------------------------- */
/*! - Move file pointer until the first line that does not
      begin with a "#".                                                  */
/* --------------------------------------------------------------------- */

    success = 0;
    while (!success) {
        fgetpos(fp, &file_pos);
        fgets(sline, 512, fp);
        if (sline[0] != '#') success = 1;
    }

    fsetpos(fp, &file_pos);

/* --------------------------------------------------------------------- */
/*! - Start reading number of points and grid coordinates. For the
      input x1 direction these are stored inside the module variables
      ::id_nx1 and ::id_x1.                                              */
/* --------------------------------------------------------------------- */

    fscanf(fp, "%d \n", &Gin->id_nx1);

    Gin->x1 = (double *) malloc(sizeof(double) * Gin->id_nx1);

    for (i = 0; i < Gin->id_nx1; i++) {
        fscanf(fp, "%d  %lf %lf\n", &ip, &xl, &xr);
        Gin->x1[i] = 0.5 * (xl + xr);
    }

    fscanf(fp, "%d \n", &Gin->id_nx2);
    Gin->x2 = (double *) malloc(sizeof(double) * Gin->id_nx2);
    for (i = 0; i < Gin->id_nx2; i++) {
        fscanf(fp, "%d  %lf %lf\n", &ip, &xl, &xr);
        Gin->x2[i] = 0.5 * (xl + xr);
    }

    fscanf(fp, "%d \n", &Gin->id_nx3);
    Gin->x3 = (double *) malloc(sizeof(double) * Gin->id_nx3);
    for (i = 0; i < Gin->id_nx3; i++) {
        fscanf(fp, "%d  %lf %lf\n", &ip, &xl, &xr);
        Gin->x3[i] = 0.5 * (xl + xr);
    }
    fclose(fp);

/* -- reset grid with 1 point -- */

    if (Gin->id_nx1 == 1) Gin->x1[0] = 0.0;
    if (Gin->id_nx2 == 1) Gin->x2[0] = 0.0;
    if (Gin->id_nx3 == 1) Gin->x3[0] = 0.0;


}

/* ***************************************************/
  void gen_cldlist(int Nclouds, int Nfiles) {
/*
 * Generate random cloud position. Radius is weighted as
 * n(r) ~ exp(-phi/sigma^2) if external tabulated gravity.
 * Else uniformly chosen.
 *
******************************************************/
    double x1c, x2c, x3c;
    double sigma, v1, v2, v3;
    long int seed1 = -1000, seed2 = -90000, seed3 = -50000;
    long int vseed1 = -5555, vseed2 = -99999, vseed3 = -222222;
    int ll, count;
    FILE *fp_cldlist;


    //---Open cloud position list file---//

    fp_cldlist = fopen("cloud_pos_list.dat", "w");
    if (fp_cldlist == NULL) {
        print("Can't open cloud_pos_list.dat \n");
        QUIT_PLUTO(1);
    }

    if (Nclouds <= 0) {
        print("PAR_NCLD should be >=1 \n");
        QUIT_PLUTO(1);
    }

/*
  #if CLOUD_DENSITY == CD_TURB_ISOTH_HYDROSTATIC

  if (gr_nr <=0) {print("Something wrong in gen_cldlist. gr_nr=%d \n"); QUIT_PLUTO(1);}

  double sigma2, hx, norm, uran;
  double rad, phi, theta, cphi, ctheta;
  int ii_rad, ii_norm, i, flip;
  double *cumul, *integrand;

  cumul=ARRAY_1D(gr_nr,double);
  integrand=ARRAY_1D(gr_nr,double);

  sigma2=g_inputParam[PAR_WTRB]*1.e5/UNIT_VELOCITY;
  sigma2*=sigma2;


  for (i=0;i<gr_nr;i++){
    integrand[i]=exp(-1.0*gr_phi[i]/sigma2);
  }

  //---Define uniform interval. Any difference would do. Choose any cell besides the first two
  //---as the first row is manually entered from the python output---
  hx=gr_r[10]-gr_r[9];


  for (i=1;i<gr_nr;i++){
    cumul[i]=Simpson_ext(hx,integrand,i+1);
  }

  //---Locate end point of distribution of clouds. Normalise integrand---//
  ii_norm=locate(gr_r,g_inputParam[PAR_WRAD],gr_nr);
  norm=cumul[ii_norm];
  for (i=1;i<gr_nr;i++){
  cumul[i]=cumul[i]/norm;
  }
  cumul[0]=0.0;

 #endif

*/
    //---Create radial distribution---//
    count = 0;
    for (ll = 0; ll < Nclouds; ll++) {

/*
  #if CLOUD_DENSITY == CD_TURB_ISOTH_HYDROSTATIC

  uran=ran1(&seed1);
  ii_rad=locate(cumul,uran,gr_nr);
  if (ii_rad > ii_norm) ii_rad=ii_norm;
  rad=gr_r[ii_rad];
  #endif

   #if CLOUD_DENSITY == CD_HOMOGENEOUS
  //---Choose radius uniformly---
  rad=g_inputParam[PAR_WRAD]*ran1(&seed1);
  #endif

  //---Choose theta and phi uniformly---//
 // theta=M_PI/2.0*ran1(&seed2);
 // phi=2.0*M_PI*ran1(&seed3);

   ctheta=ran1(&seed2);
   theta=acos(ctheta);

   cphi=-1.0+2.0*ran1(&seed3);
   if (ran1(&seed4) > 0.5) flip=1; else flip=0;
   phi=acos(cphi)+M_PI*flip;


  //---Transform to Cartesian-----//
  x1c=SPH_1(rad,theta,phi);
  x2c=SPH_2(rad,theta,phi);
  x3c=SPH_3(rad,theta,phi);
 
*/

        /*----DM 10Aug15: Current setting is to uniformly distribute
               clouds within the domain. Apodization with external
           potential shapes the distribution. Above commented
           code generates random radial distribution weighted by potential along r,
           uniform in theta and phi. But multiplying multiple random variates
           skews net distribution. To be changed later on.     */

        x1c = (g_inputParam[PAR_WX1L] + (g_inputParam[PAR_WX1H] - g_inputParam[PAR_WX1L]) * ran1(&seed1));
        x2c = (g_inputParam[PAR_WX2L] + (g_inputParam[PAR_WX2H] - g_inputParam[PAR_WX2L]) * ran1(&seed2));
        x3c = (g_inputParam[PAR_WX3L] + (g_inputParam[PAR_WX3H] - g_inputParam[PAR_WX3L]) * ran1(&seed3));

        sigma = g_inputParam[PAR_SGAV];
        v1 = gasdev(&vseed1) * sigma;
        v2 = gasdev(&vseed2) * sigma;
        v3 = gasdev(&vseed3) * sigma;


        if (count > Nfiles - 1) count = 0;
        fprintf(fp_cldlist, "%d %lf %lf %lf %lf %lf %lf \n", count, x1c, x2c, x3c, v1, v2, v3);
        count++;
    } // ll


    fclose(fp_cldlist);

/*
  #if CLOUD_DENSITY == CD_TURB_ISOTH_HYDROSTATIC
  FreeArray1D(cumul);
  FreeArray1D(integrand);
  #endif
*/

}//end subroutine


#endif //MULTI_CLOUDS


/* ****************************************************************************/
double Simpson_ext(double delta, double *y, int ni) {
/*
 * Integrate a tabulated array, uniformly spaced. For ni > 6 follows the extended 
 * Simpson rule. Lower ni follows different Newton Coyote formulae.
******************************************************************************/
	double sum;
	double f1, f2, f3, f4, f5, f6;
	int i;

	sum = 0.;

	if (ni <= 1) {
		printf("ni should be atleast 2 \n");
		exit(1);
	}

	switch (ni) {

		case 2:
			sum = (y[0] + y[1]) / 2.;
			break;

		case 3:
			f1 = y[0];
			f2 = y[1];
			f3 = y[2];
			sum = 1. / 3. * f1 + 4. / 3. * f2 + 1. / 3. * f3;
			break;

		case 4:
			f1 = y[0];
			f2 = y[1];
			f3 = y[2];
			f4 = y[3];
			sum = 3. / 8. * f1 + 9. / 8. * f2 + 9. / 8. * f3 + 3. / 8. * f4;
			break;

		case 5:
			f1 = y[0];
			f2 = y[1];
			f3 = y[2];
			f4 = y[3];
			f5 = y[4];
			sum = 14. / 45. * f1 + 64. / 45. * f2 + 24. / 45. * f3 + 64. / 45. * f4 + 14. / 45. * f5;
			break;

		case 6:
			f1 = y[0];
			f2 = y[1];
			f3 = y[2];
			f4 = y[3];
			f5 = y[4];
			f6 = y[5];
			sum = 5. / 288. * (19. * f1 + 75. * f2 + 50. * f3 + 50. * f4 + 75. * f5 + 19. * f6);
			break;

		default:
			//----Eq. 4.1.14 from NR
			sum = 3. / 8. * y[0] + 7. / 6. * y[1] + 23. / 24. * y[2] + 23. / 24. * y[ni - 3] +
				  7. / 6. * y[ni - 2] + 3. / 8. * y[ni - 1];
			for (i = 3; i < ni - 3; i++) sum += y[i];
			break;
	}

	return sum * delta;
}

