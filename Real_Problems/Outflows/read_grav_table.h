#ifndef read_grav_table_h
#define read_grav_table_h
#ifdef GRAV_TABLE
/* Make sure if included elsewhere it is 
 * preceded by #include pluto.h */

 //------(DM: Feb20, 2015): for GRAV=POTENTIAL read a three column table of r, phi, dphi/dr normalised to PLUTO units before read------

/* functions */
void ReadGravTable();

/* global variables*/
extern double *gr_rad, *gr_phi, *gr_vec, *gr_dphidr;
extern int gr_ndata;

#endif
#endif
