#include "pluto.h"
#include "pluto_usr.h"
#include "rMuTable.h"

#ifdef MU_TABLE

double *mu_por, *mu_mu;
int mu_ndata;

int rMuData(){
  /*
   * This routine reads the data of p/rho mu from data file. 
   *
   * The data should be a two file in cgs units. 
   * Columns are p/rho*KELVIN (=T/mu), density, and pressure.
   * 
   * The values of these are filled into the global arrays mu_por
   * and mu_mu (por = p over rho) are used throughout the code.
   *
   * Returns 0
   *
   * */

  FILE * f;

  double buf;
  int i;

  /* Open file */
  if ((f = fopen(MUFNAME, "r")) == NULL){
    print1("Error: rMuData: Unable to open file");
    exit(1);
  }

  /* Scan file first to get number of lines*/
  mu_ndata = 0;
  while (fscanf(f, "%le %le", &buf, &buf) != EOF){
    mu_ndata++;
  }

  /* Allocate memory for potential profile arrays */
  mu_por = ARRAY_1D(mu_ndata, double);
  mu_mu = ARRAY_1D(mu_ndata, double);

  /* Read data */
  fseek(f, 0, SEEK_SET);
  for (i=0; i<mu_ndata; ++i){
    fscanf(f, "%le ", &mu_por[i]);
    fscanf(f, "%le ", &mu_mu[i]);
  }

  /* Clean up */
  fclose(f);

  return 0;
}


void readMuTable(){
  /* This routine currently just calls the data reading routine */

  /* Read Mean molecular Mass table */
  if (rMuData() != 0){
    print1("Error: in rMuTable.c: Table read unsuccessful.");
    QUIT_PLUTO(1);
  }

}
