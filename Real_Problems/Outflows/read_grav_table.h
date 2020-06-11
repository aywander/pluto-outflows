#ifndef read_grav_table_h
#define read_grav_table_h
#ifdef GRAV_TABLE
/* Make sure if included elsewhere it is 
 * preceded by #include pluto.h */

/* functions */
void ReadGravTable();

/* global variables*/
extern double *gr_r, *gr_z;
extern int gr_nr, gr_nz;

#ifdef GRAV_2D_POTENTIAL
extern double **gr_phi, **gr_acc_r, **gr_acc_z;;
#else
extern double *gr_phi, *gr_acc_r;
#endif

#endif
#endif
