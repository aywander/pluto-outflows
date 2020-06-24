#include "pluto.h"
#include "pluto_usr.h"
#include "read_grav_table.h"
#include "init_tools.h"

#ifdef GRAV_TABLE

double *gr_r, *gr_z;
#ifdef GRAV_2D_POTENTIAL
double **gr_phi, **gr_acc_r, **gr_acc_z;
#else
double *gr_phi, *gr_acc_r;
#endif

int gr_nr, gr_nz;

/* ************************************************ */
void ReadGravTable() {
/*!
 * This routine reads the data from a gravity file.
 *
 * The values of these are filled into the global arrays gr_r, gr_phi, and gr_acc_r, which
 * are used throughout the code. If GRAV_2D_POTENTIAL is true, then furthermore, the arrays
 * gr_z, gr_acc_z get filled.
 *
 * NOTE:
 *   - GRAV_2D_POTENTIAL, just like GRAV_TABLE, is set automatically with the choice of GRAV_POTENTIAL
 *   - In 2D arrays, the indices are z and r (i.e. r changing fastest).
 *   - If BODY_FORCE & VECTOR, then the gravitational acceleration must be provided as a table.
 *     The gravitational potential, however, must always be provided in tabular form.
 *     This may change in the future.
 *
 ************************************************** */

    FILE *fg, *fr, *fz, *fgr, *fgz;

    double buf;
    int i, j;

    /* Open r-coordinates file */
    if ((fg = fopen(GRAV_R_FNAME, "r")) == NULL) {
        print("Error: ReadGravData: Unable to open r-domain data file");
        exit(1);
    }

    /* Scan file first to get number of r coordinates */
    gr_nr = 0;
    while (fscanf(fg, "%le", &buf) != EOF)
        gr_nr++;
    gr_r = ARRAY_1D(gr_nr, double);

    /* Read r-coordinates */
    fseek(fg, 0, SEEK_SET);
    for (i = 0; i < gr_nr; ++i) {
        fscanf(fg, "%le ", &gr_r[i]);
        gr_r[i] /= vn.l_norm;
    }
    fclose(fg);

    /* If we're using a 2D potential do the same for the z-coordinate */
#ifdef GRAV_2D_POTENTIAL

    /* Open z-coordinates file */
    if ((fg = fopen(GRAV_Z_FNAME, "r")) == NULL) {
        print("Error: ReadGravData: Unable to open z-domain data file");
        exit(1);
    }

    /* Scan file first to get number of z coordinates */
    gr_nz = 0; while (fscanf(fg, "%le", &buf) != EOF) gr_nz++;
    gr_z = ARRAY_1D(gr_nz, double);

    /* Read z-coordinates */
    fseek(fg, 0, SEEK_SET);
    for (j = 0; j < gr_nz; ++j) {
        fscanf(fg, "%le ", &gr_z[j]);
        gr_z[j] /= vn.l_norm;
    }
    fclose(fg);

#endif

    /* Open potential data file */
    if ((fg = fopen(GRAV_PHI_FNAME, "r")) == NULL) {
        print("Error: ReadGravData: Unable to open gravity data file");
        exit(1);
    }

    /* The data is stored so that r changes fastest (second index) */
#ifdef GRAV_2D_POTENTIAL
    gr_phi = ARRAY_2D(gr_nz, gr_nr, double);
    for (j = 0; j < gr_nz; ++j) {
       for (i = 0; i < gr_nr; ++i) {
        fscanf(fg, "%le ", &gr_phi[j][i]);
        gr_phi[j][i] /= vn.pot_norm;
      }
    }
#else
    gr_phi = ARRAY_1D(gr_nr, double);
    for (i = 0; i < gr_nr; ++i) {
        fscanf(fg, "%le ", &gr_phi[i]);
        gr_phi[i] /= vn.pot_norm;
    }
#endif
    fclose(fg);

    /* If using acceleration vector, also read acceleration files */

    gr_acc_r = ARRAY_2D(gr_nz, gr_nr, double);
    gr_acc_z = ARRAY_2D(gr_nz, gr_nr, double);
#if BODY_FORCE & VECTOR
    if ((fgr = fopen(GRAV_ACCR_FNAME, "r")) == NULL) {
        print("Error: ReadGravData: Unable to open r-acceleration data file");
        exit(1);
    }
    if ((fgz = fopen(GRAV_ACCZ_FNAME, "r")) == NULL) {
        print("Error: ReadGravData: Unable to open z-acceleration data file");
        exit(1);
    }
#ifdef GRAV_2D_POTENTIAL
    for (i = 0; i < gr_nr; ++i) {
        for (j = 0; j < gr_nz; ++j) {
            fscanf(fgr, "%le ", &gr_acc_r[i][j]);
            fscanf(fgz, "%le ", &gr_acc_z[i][j]);
            gr_acc_r[i][j] /= vn.acc_norm;
            gr_acc_z[i][j] /= vn.acc_norm;
        };
    }
    fclose(fgz);
#else
    gr_acc_r = ARRAY_1D(gr_nr, double);
    for (i = 0; i < gr_nr; ++i) {
        fscanf(fgr, "%le ", &gr_acc_r[i]);
        gr_acc_r[i] /= vn.acc_norm;
    }
#endif // ifdef GRAV_2D_POTENTIAL
    fclose(fgr);

    /* If we're not using acceleration vector, calculate derivatives from potential here */
#else

#ifdef GRAV_2D_POTENTIAL

    /* For transposing arrays */
    double **gr_acc_z_T, **gr_phi_T;
    gr_acc_z_T = ARRAY_2D(gr_nr, gr_nz, double);
    gr_phi_T = ARRAY_2D(gr_nr, gr_nz, double);
    TransposeArray(gr_phi, gr_phi_T, gr_nz, gr_nr);

    /* Calculate gradients row-wise */
    for (j = 0; j < gr_nz; ++j) GradFromArray(gr_phi[j], gr_r, gr_acc_r[j], gr_nr);
    for (i = 0; i < gr_nr; ++i) GradFromArray(gr_phi_T[i], gr_z, gr_acc_z_T[i], gr_nz);

    /* Transpose r array back again */
    TransposeArray(gr_acc_z_T, gr_acc_z, gr_nr, gr_nz);

    /* Make acceleration vector point inward (dphi / dr = -g) */
    for (j = 0; j < gr_nz; ++j) {
        for (i = 0; i < gr_nr; ++i) {
            gr_acc_r[j][i] = -gr_acc_r[j][i];
            gr_acc_z[j][i] = -gr_acc_z[j][i];
        }
    }


#else
    /* Gradient of potential */
    GradFromArray(gr_phi, gr_r, gr_acc_r, gr_nr);

    /* Make acceleration vector point inward (dphi / dr = -g) */
    for (i = 0; i < gr_nr; ++i) gr_acc_r[i] = -gr_acc_r[i];

#endif

#endif // if BODY_FORCE & VECTOR

}

#endif // ifdef GRAV_TABLE
