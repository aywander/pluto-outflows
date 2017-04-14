#include "pluto.h"
#include "pluto_usr.h"
#include "read_grav_table.h"
#include "init_tools.h"

#ifdef GRAV_TABLE

//double *gr_rad, *gr_phi, *gr_vec;
double *gr_rad, *gr_phi;
double *gr_dphidr;
int gr_ndata;

void ReadGravTable() {
  /*
   * This routine reads the data from a gravity file. The data should be a
   * two-column file in code units. Columns are radius and potential or gravitational
   * acceleration, depending on whether GRAV_USE_POTENTIAL is YES or NO.
   *
   * The values of these are filled into the global arrays gr_rad, gr_phi, and gr_dphidr, which
   * are used throughout the code.
   * 
   * Returns 0
   *
   * */

  FILE *f;

  double buf;
  int i;

  /* Open file */
  if ((f = fopen(GRAV_FNAME, "r")) == NULL) {
    print1("Error: ReadGravData: Unable to open file");
    exit(1);
  }

  /* Scan file first to get number of lines*/
  gr_ndata = 0;
  while (fscanf(f, "%le %le %le", &buf, &buf, &buf) != EOF) {
    gr_ndata++;
  }

  /* Allocate memory for potential profile arrays */
  gr_rad = ARRAY_1D(gr_ndata, double);
//#if BODY_FORCE == POTENTIAL
  gr_phi = ARRAY_1D(gr_ndata, double);
  gr_dphidr = ARRAY_1D(gr_ndata, double);
//#else
//  gr_vec = ARRAY_1D(gr_ndata, double);
//#endif

  /* Read data */
  fseek(f, 0, SEEK_SET);
  for (i = 0; i < gr_ndata; ++i) {
    fscanf(f, "%le ", &gr_rad[i]);
//#if BODY_FORCE == POTENTIAL
    fscanf(f, "%le ", &gr_phi[i]);
    fscanf(f, "%le ", &gr_dphidr[i]);
//#else
//    fscanf(f, "%le ", &gr_vec[i]);
//#endif
  }

  /* Clean up */
  fclose(f);

  /* Convert variables into code units */
  for (i = 0; i < gr_ndata; ++i) {
    gr_rad[i] /= vn.l_norm;
    gr_phi[i] /= vn.pot_norm;
    gr_dphidr[i] /= vn.pot_norm / vn.l_norm;
  }


}


#endif
