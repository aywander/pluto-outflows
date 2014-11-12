#include "pluto.h"
#include "pluto_usr.h"
#include "rHotTable.h"

#if HOT_DISTR == HOT_EXT_DATA

double *hot_rad, *hot_rho, *hot_prs;
int hot_ndata;

int rHotData(){
  /*
   * This routine reads the data from a hot phase data file. 
   *
   * The data should be a three-column file in code units. 
   * Columns are radius, density, and pressure.
   * 
   * The values of these are filled into the global arrays hot_rad, hot_rho and hot_prs, which
   * are used throughout the code.
   *
   * Returns 0
   *
   * */

  FILE * f;

  double buf;
  int i;

  /* Open file */
  if ((f = fopen(HOTFNAME, "r")) == NULL){
    print1("Error: rHotData: Unable to open file");
    exit(1);
  }

  /* Scan file first to get number of lines*/
  hot_ndata = 0;
  while (fscanf(f, "%le %le %le", &buf, &buf, &buf) != EOF){
    hot_ndata++;
  }

  /* Allocate memory for profile arrays */
  hot_rad = ARRAY_1D(hot_ndata, double);
  hot_rho = ARRAY_1D(hot_ndata, double);
  hot_prs = ARRAY_1D(hot_ndata, double);

  /* Read data */
  fseek(f, 0, SEEK_SET);
  for (i=0; i<hot_ndata; ++i){
    fscanf(f, "%le ", &hot_rad[i]);
    fscanf(f, "%le ", &hot_rho[i]);
    fscanf(f, "%le ", &hot_prs[i]);
  }

  /* Clean up */
  fclose(f);

  return 0;
}

void readHotTable(){
  /* This routine currently just calls the data reading routine */

  /* Read gravity table */
  if (rHotData() != 0){
    print1("Error: in rHotTable.c: Table read unsuccessful.");
    QUIT_PLUTO(1);
  }

}

#endif
