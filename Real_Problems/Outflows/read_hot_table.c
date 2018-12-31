#include "pluto.h"
#include "pluto_usr.h"
#include "read_hot_table.h"
#include "init_tools.h"

#ifdef GRAV_TABLE

double *hot_rad, *hot_rho, *hot_prs;
int hot_ndata;

void ReadHotTable() {
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

    FILE *f;

    double buf;
    int i;

    /* Open file */
    if ((f = fopen(HOT_FNAME, "r")) == NULL) {
        print("Error: ReadHotData: Unable to open file");
        exit(1);
    }

    /* Scan file first to get number of lines*/
    hot_ndata = 0;
    while (fscanf(f, "%le %le %le", &buf, &buf, &buf) != EOF) {
        hot_ndata++;
    }

    /* Allocate memory for profile arrays */
    hot_rad = ARRAY_1D(hot_ndata, double);
    hot_rho = ARRAY_1D(hot_ndata, double);
    hot_prs = ARRAY_1D(hot_ndata, double);

    /* Read data */
    fseek(f, 0, SEEK_SET);
    for (i = 0; i < hot_ndata; ++i) {
        fscanf(f, "%le ", &hot_rad[i]);
        fscanf(f, "%le ", &hot_rho[i]);
        fscanf(f, "%le ", &hot_prs[i]);
    }

    /* Clean up */
    fclose(f);

    /* Convert variables into code units */
    for (i = 0; i < hot_ndata; ++i) {
        hot_rad[i] /= vn.l_norm;
        hot_rho[i] /= vn.dens_norm;
        hot_prs[i] /= vn.pres_norm;
    }

}

#endif
