#ifndef rGravTable_h
#define rGravTable_h
#ifdef GRAV_TABLE
/* Make sure if included elsewhere it is 
 * preceded by #include pluto.h */

/* functions */
int rGravData();
int rGravHeader();
void readGravTable();

/* global variables*/
extern double *gr_rad, *gr_phi, *gr_vec;
extern double potential[NPOT];
extern int gr_ndata;

#endif
#endif
