#ifndef rMuTable_h
#define rMuTable_h
#ifdef MU_TABLE
/* Make sure if included elsewhere it is 
 * preceded by #include pluto.h */

/* functions */
int rMuData();
void readMuTable();

/* global variables*/
/* mu_por is p/rho */
extern double *mu_por, *mu_mu;
extern int mu_ndata;

#endif
#endif

