#ifndef read_hot_table_h
#define read_hot_table_h
#ifdef GRAV_TABLE
/* Make sure if included elsewhere it is 
 * preceded by #include pluto.h */

/* functions */
void ReadHotTable();

/* global variables*/
extern double *hot_rad, *hot_rho, *hot_prs;
extern int hot_ndata;

#endif
#endif

