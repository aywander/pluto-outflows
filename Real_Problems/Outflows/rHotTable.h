#if HOT_DISTR == HOT_EXT_DATA
/* Make sure if included elsewhere it is 
 * preceded by #include pluto.h */

/* functions */
int rHotData();
void readHotTable();

/* global variables*/
extern double *hot_rad, *hot_rho, *hot_prs;
extern int hot_ndata;

#endif

