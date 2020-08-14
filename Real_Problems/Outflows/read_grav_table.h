#ifndef read_grav_table_h
#define read_grav_table_h

#if GRAV_POTENTIAL == GRAV_TABLE || GRAV_POTENTIAL == GRAV_2D_TABLE
/* Make sure if included elsewhere it is
 * preceded by #include pluto.h */

/* functions */
void ReadGravTable();

/* global variables*/
extern double *gr_r, *gr_z;
extern int gr_nr, gr_nz;

#if GRAV_POTENTIAL == GRAV_2D_TABLE
extern double **gr_phi, **gr_acc_r, **gr_acc_z;;
#else
extern double *gr_phi, *gr_acc_r;
#endif

#endif
#endif
