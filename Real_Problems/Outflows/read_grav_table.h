#ifndef read_grav_table_h
#define read_grav_table_h
#ifdef GRAV_TABLE
/* Make sure if included elsewhere it is 
 * preceded by #include pluto.h */

/* functions */
void ReadGravTable();

/* global variables*/
//extern double *gr_rad, *gr_phi, *gr_vec, *gr_dphidr;
extern double *gr_rad, *gr_phi, *gr_dphidr;
extern int gr_ndata;

#endif
#endif
